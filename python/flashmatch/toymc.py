import numpy as np
from flashmatch import flashmatch, geoalgo
import sys

np.random.seed(0)


def generate_random_point(det):
    """
    Generates a random position inside the detector.
    """
    xmin, ymin, zmin = det.ActiveVolume().Min()
    xmax, ymax, zmax = det.ActiveVolume().Max()
    return geoalgo.Vector(np.random.random() * (xmax - xmin) + xmin,
                          np.random.random() * (ymax - ymin) + ymin,
                          np.random.random() * (zmax - zmin) + zmin)


def generate_random_time(tmin=-1000, tmax=1000):
    # in micro-seconds
    time = np.random.random() * (tmax - tmin) + tmin
    return time


# dumb function to generate a stupid trajectory
def gen_trajectories(num_tracks, det):
    """
    Generate N random trajectories.

    Arguments
    ---------
    num_tracks: int
    det: instance of flashmatch.DetectorSpecs.GetME()
    """
    res = []
    for i in range(num_tracks):
        trj = geoalgo.Trajectory(2, 3)  # 2 points, 3d
        trj[0] = generate_random_point(det)
        trj[1] = generate_random_point(det)
        res.append(trj)
    return res


# Create a set ot TPC (flashmatch::QCluster) and PMT (flashmatch::Flash_t) interactions
def gen_input(num_tracks, blob, time_algo="random"):
    track_v = gen_trajectories(num_tracks, blob.det)
    tpc_v = flashmatch.QClusterArray_t()
    pmt_v = flashmatch.FlashArray_t()
    time_v = []
    for idx, track in enumerate(track_v):
        # Decide time range to generate, completely made up for now
        # in us
        if time_algo == 'random':
            time = generate_random_time()
        elif time_algo == 'x':
            x_distance = (track[0][0]-blob.det.ActiveVolume().Min()[0])  # cm
            # 3e-3 us to cover 1m
            time = x_distance * 3e-5
        elif time_algo == 'sequential':
            time = idx * 10
        else:
            raise Exception("Time algo not recognized, must be one of ['random', 'x', 'sequential']")
        time_v.append(time)

        # TPC xyz & light info
        qcluster = blob.make_qcluster(track, idx)
        # PMT pe spectrum
        flash = blob.make_flash(qcluster, time, idx)
        # save
        tpc_v.push_back(qcluster)
        pmt_v.push_back(flash)

    return track_v, tpc_v, pmt_v, time_v


class ToyMC:
    def __init__(self, cfg_file, CustomAlgo="LightPath"):
        self.mgr = flashmatch.FlashMatchManager()
        self.det = flashmatch.DetectorSpecs.GetME()
        self.cfg = flashmatch.CreatePSetFromFile(cfg_file)
        # configure
        self.mgr.Configure(self.cfg)

        # Define algorithms
        self.tpc_algo = self.mgr.GetCustomAlgo(CustomAlgo)
        # PhotonLibHypothesis
        self.pmt_algo = self.mgr.GetAlgo(flashmatch.kFlashHypothesis)

    def make_qcluster(self, track, idx):
        qcluster = self.tpc_algo.FlashHypothesis(track)
        qcluster.idx = idx
        return qcluster

    def make_flash(self, qcluster, time, idx):
        flash = flashmatch.Flash_t()
        flash.idx  = idx
        flash.time = time
        self.pmt_algo.FillEstimate(qcluster, flash)
        # very crude error addition (at least we should change it to poisson)
        for idx,npe in enumerate(flash.pe_v):
            flash.pe_err_v[idx]=np.sqrt(npe) # only good for high npe
        return flash

    def add(self, v):
        for elt in v:
            self.mgr.Add(elt)

    def match(self, tpc_v, pmt_v):
        self.mgr.Reset()
        self.add(tpc_v)
        self.add(pmt_v)
        return self.mgr.Match()


def demo(cfg_file, num_tracks):
    blob = ToyMC(cfg_file)
    track_v, tpc_v, pmt_v, time_v = gen_input(num_tracks, blob)
    match_v = blob.match(tpc_v, pmt_v)

    # Print history of the last minimization if stored
    qll=blob.mgr.GetAlgo(flashmatch.kFlashMatch)
    hist_llhd = qll.HistoryLLHD()
    hist_chi2 = qll.HistoryChi2()
    hist_xpos = qll.HistoryX()
    if hist_llhd.size():
        print('QLLMatch history')
        np.savetxt('hist_llhd.csv', np.array(hist_llhd), delimiter=',')
        np.savetxt('hist_chi2.csv', np.array(hist_chi2), delimiter=',')
        np.savetxt('hist_xpos.csv', np.array(hist_xpos), delimiter=',')
        for i in range(hist_xpos.size()):
            print('Step %d X %f Chi2 %f LLHD %f' % (i,hist_xpos[i],hist_chi2[i],hist_llhd[i]))

    print('Number of match result',len(match_v))
    for idx,match in enumerate(match_v):
        print('Match ID',idx)
        print('  TPC/PMT IDs %d/%d Score %f Min-X %f' % (match.tpc_id,match.flash_id,match.score,match.tpc_point.x))

    #tpc_v = [flashmatch.as_ndarray(tpc) for tpc in tpc_v]
    #pmt_v = [flashmatch.as_ndarray(pmt) for pmt in pmt_v]

    tpc_id_v = [tpc.idx for tpc in tpc_v]
    pmt_id_v = [pmt.idx for pmt in pmt_v]
    all_matches = []
    for match in match_v:
        # FIXME is id same as order in tpc_v?
        qcluster = tpc_v[tpc_id_v.index(match.tpc_id)]
        flash = pmt_v[pmt_id_v.index(match.flash_id)]
        store = np.array([[
            match.score,
            match.tpc_point.x,
            match.tpc_point.y,
            match.tpc_point.z,
            track_v[qcluster.idx][0][0],
            track_v[qcluster.idx][0][1],
            track_v[qcluster.idx][0][2],
            track_v[qcluster.idx][1][0],
            track_v[qcluster.idx][1][1],
            track_v[qcluster.idx][1][2],
            len(qcluster),
            int(qcluster.idx == flash.idx),
            qcluster.sum(),
            np.sum(match.hypothesis),
            flash.TotalPE(),
            match.duration
        ]])
        #)]], dtype=[
        #    ('score', 'f4'),
        #    ('tpc_point_x', 'f4'),
        #    ('tpc_point_y', 'f4'),
        #    ('tpc_point_z', 'f4'),
        #    ('start_point_x', 'f4'),
        #    ('start_point_y', 'f4'),
        #    ('start_point_z', 'f4'),
        #    ('end_point_x', 'f4'),
        #    ('end_point_y', 'f4'),
        #    ('end_point_z', 'f4'),
        #    ('qcluster_num_points', 'f4'),
        #    ('matched', 'B'),
        #    ('qcluster_sum', 'f4'),
        #    ('flash_sum', 'f4'),
        #    ('duration', 'f4')
        #])
        print(store)
        print(store.shape)
        all_matches.append(store)
    x = np.concatenate(all_matches, axis=0)
    print(x.shape)
    names = [
        'score',
        'tpc_point_x',
        'tpc_point_y',
        'tpc_point_z',
        'start_point_x',
        'start_point_y',
        'start_point_z',
        'end_point_x',
        'end_point_y',
        'end_point_z',
        'qcluster_num_points',
        'matched',
        'qcluster_sum',
        'hypothesis_sum',
        'flash_sum',
        'duration'
    ]
    np.savetxt('matches.csv', x, delimiter=',', header=','.join(names))


if __name__ == '__main__':
    import os

    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                            'dat', 'flashmatch.cfg')
    num_tracks = 1
    if len(sys.argv) > 1:
        for argv in sys.argv[1:]:
            if argv.isdigit():
                num_tracks = int(argv)
            else:
                cfg_file = sys.argv[1]

    # Generate some random tracks
    demo(cfg_file, num_tracks)
