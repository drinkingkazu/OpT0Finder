import numpy as np
from flashmatch import flashmatch, geoalgo
import sys


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
    flash = blob.match(tpc_v, pmt_v)
    print(len(flash))
    tpc_v = [flashmatch.as_ndarray(tpc) for tpc in tpc_v]
    pmt_v = [flashmatch.as_ndarray(pmt) for pmt in pmt_v]


if __name__ == '__main__':
    import os

    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                            'dat', 'flashmatch.cfg')
    if len(sys.argv) > 1:
        cfg_file = sys.argv[1]

    # Generate 10 random tracks
    demo(cfg_file, 10)
