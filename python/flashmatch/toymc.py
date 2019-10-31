import numpy as np
from flashmatch import flashmatch, geoalgo
import sys

class ToyMC:
    def __init__(self, cfg_file=None):
        self.mgr = flashmatch.FlashMatchManager()
        self.det = flashmatch.DetectorSpecs.GetME()
        self._qcluster_algo  = None
        self._flash_algo     = None
        self._time_algo      = None
        self._track_algo     = None
        self._ly_variation = 0.0
        self._pe_variation = 0.0
        if not cfg_file is None:
            self.configure(cfg_file)

    def configure(self,cfg_file):
        self.cfg = flashmatch.CreatePSetFromFile(cfg_file)
        # configure
        self.mgr.Configure(self.cfg)
        # PhotonLibHypothesis
        self._flash_algo = self.mgr.GetAlgo(flashmatch.kFlashHypothesis)
        #
        # Toy MC configuration
        #
        pset = self.cfg.get('flashmatch::PSet')('ToyMC')
        # LightPath
        self._qcluster_algo = self.mgr.GetCustomAlgo(pset.get('QClusterAlgo'))
        # PE variation on (fake) reconstructed flash
        self._pe_variation = float(pset.get('PEVariation'))
        # Photon estimation variation
        self._ly_variation = float(pset.get('LightYieldVariation'))
        # time algorithm keyword
        self._time_algo = str(pset.get('TimeAlgo'))
        # track algorithm keyword
        self._track_algo = str(pset.get('TrackAlgo'))
        # Number of tracks to generate at once
        self._num_tracks = int(pset.get('NumTracks'))
        # Set seed if there is any specified
        if pset.contains_value('NumpySeed'):
            seed = int(pset.get('NumpySeed'))
            if seed < 0:
                import time
                seed = int(time.time())
            np.random.seed(seed)


    # Create a set ot TPC (flashmatch::QCluster) and PMT (flashmatch::Flash_t) interactions
    def gen_input(self,num_match=None):
        """
        Make N input samples for flash matching
        ---------
        Arguments
          num_match: int, number of tpc-pmt pair samples to be generated (optional, default self._num_tracks)
        -------
        Returns
          a list of geoalgo::Trajectory, flashmatch::QClusterArray_t, and flashmatch::FLashArray_t
        """
        if num_match is None:
            num_match = self._num_tracks

        track_v = self.gen_trajectories(num_match)
        tpc_v = flashmatch.QClusterArray_t()
        pmt_v = flashmatch.FlashArray_t()
        for idx, track in enumerate(track_v):
            # TPC xyz & light info
            qcluster = self.make_qcluster(track)
            # PMT pe spectrum
            flash = self.make_flash(qcluster)
            # set IDs
            qcluster.idx = idx
            flash.idx    = idx
            # save
            tpc_v.push_back(qcluster)
            pmt_v.push_back(flash)
            
        return track_v, tpc_v, pmt_v


    def gen_trajectories(self,num_tracks):
        """
        Generate N random trajectories.
        ---------
        Arguments
          num_tracks: int, number of tpc trajectories to be generated
        -------
        Returns
          a list of geoalgo::Trajectory objects
        """
        xmin, ymin, zmin = self.det.ActiveVolume().Min()
        xmax, ymax, zmax = self.det.ActiveVolume().Max()
        res = []
        for i in range(num_tracks):
            trj = geoalgo.Trajectory(2, 3)  # 2 points, 3d
            if self._track_algo=="random":
                trj[0] = geoalgo.Vector(np.random.random() * (xmax - xmin) + xmin,
                                        np.random.random() * (ymax - ymin) + ymin,
                                        np.random.random() * (zmax - zmin) + zmin)
                trj[1] = geoalgo.Vector(np.random.random() * (xmax - xmin) + xmin,
                                        np.random.random() * (ymax - ymin) + ymin,
                                        np.random.random() * (zmax - zmin) + zmin)
            elif self._track_algo=="top-bottom":
                trj[0] = geoalgo.Vector(np.random.random() * (xmax - xmin) + xmin,
                                    ymax,
                                        np.random.random() * (zmax - zmin) + zmin)
                trj[1] = geoalgo.Vector(np.random.random() * (xmax - xmin) + xmin,
                                        ymin,
                                        np.random.random() * (zmax - zmin) + zmin)
            else:
                raise Exception("Track algo not recognized, must be one of ['random', 'top-bottom']")
            res.append(trj)
        return res


    def make_qcluster(self, track):
        """
        Create a flashmatch::QCluster_t instance from geoalgo::Trajectory
        ---------
        Arguments
          track: geoalgo::Trajectory object
        -------
        Returns
          a flashmatch::QCluster_t instance
        """
        qcluster = self._qcluster_algo.MakeQCluster(track)
        # apply variation if needed
        if self._ly_variation >0.:
            var = abs(np.random.normal(1.0,self._ly_variation,qcluster.size()))
            for idx in range(qcluster.size()): qcluster[idx].q *= var[idx]

        #for idx in range(qcluster.size()): qcluster[idx].q *= 0.001
        return qcluster

    def make_flash(self, qcluster):
        """
        Create a flashmatch::Flash_t instance from flashmatch::QCluster_t
        ---------
        Arguments
          qcluster: a flashmatch::QCluster_t instance
        -------
        Returns
          a flashmatch::Flash_t instance
        """
        flash = flashmatch.Flash_t()
        self._flash_algo.FillEstimate(qcluster, flash)
        # apply variation if needed + compute error
        var = np.ones(shape=(flash.pe_v.size()),dtype=np.float32)
        if self._pe_variation>0.:
            var = abs(np.random.normal(1.0,self._pe_variation,flash.pe_v.size()))
        for idx in range(flash.pe_v.size()):
            estimate = float(int(np.random.poisson(flash.pe_v[idx] * var[idx])))
            flash.pe_v[idx] = estimate
            flash.pe_err_v[idx]=np.sqrt(estimate) # only good for high npe
        # Decide time range to generate, completely made up for now
        # in us
        if self._time_algo == 'random':
            # random -1e3 to 1e3
            flash.time = np.random.random() * (2000.) + 1000.
        elif self._time_algo == 'x':
            x_distance = (track[0][0] - self.det.ActiveVolume().Min()[0])  # cm
            # 3e-3 us to cover 1m
            flash.time = x_distance * 3e-5
        else:
            raise Exception("Time algo not recognized, must be one of ['random', 'x']")
        return flash

    def match(self, tpc_v, pmt_v):
        """
        Run flash matching for the given arrays of flashmatch::Flash_t and flashmatch::QCluster_t
        ---------
        Arguments
          tpc_v: an array of flashmatch::QCluster_t
          pmt_v: an array of flashmatch::Flash_t
        -------
        Returns
          std::vector<flashmatch::FlashMatch_t>
        """
        self.mgr.Reset()
        for tpc in tpc_v:
            self.mgr.Add(tpc)
        for pmt in pmt_v:
            self.mgr.Add(pmt)
        return self.mgr.Match()

    def print_history(self):
        """
        Print history of the last minimization if stored
        """
        qll=self.mgr.GetAlgo(flashmatch.kFlashMatch)
        hist_llhd = qll.HistoryLLHD()
        hist_chi2 = qll.HistoryChi2()
        hist_xpos = qll.HistoryX()
        if hist_llhd.size():
            print('QLLMatch history')
            #np.savetxt('hist_llhd.csv', np.array(hist_llhd), delimiter=',')
            #np.savetxt('hist_chi2.csv', np.array(hist_chi2), delimiter=',')
            #np.savetxt('hist_xpos.csv', np.array(hist_xpos), delimiter=',')
            for i in range(hist_xpos.size()):
                print('Step %d X %f Chi2 %f LLHD %f' % (i,hist_xpos[i],hist_chi2[i],hist_llhd[i]))

def demo(cfg_file,repeat=1,num_tracks=None,out_file=''):
    """
    Run function for ToyMC
    ---------
    Arguments
      cfg_file:   string for the location of a config file
      repeat:     int, number of times to run the toy MC simulation
      out_file:   string for an output analysis csv file path (optional)
      num_tracks: int for number of tracks to be generated (optional) 
    """
    mgr = ToyMC(cfg_file)
    # dump config
    print(mgr.cfg.dump())
    # override number of tracks to simulate
    if num_tracks is not None:
        num_tracks = int(num_tracks)

    np_result = None
    for event in range(repeat):
        print('Event',event,'/',repeat)
        # Generate samples
        track_v, tpc_v, pmt_v = mgr.gen_input(num_tracks)
        # Run matching
        match_v = mgr.match(tpc_v, pmt_v)
    
        #tpc_v = [flashmatch.as_ndarray(tpc) for tpc in tpc_v]
        #pmt_v = [flashmatch.as_ndarray(pmt) for pmt in pmt_v]

        # If no analysis output saving option given, return
        if not out_file:
            print('Number of match result',len(match_v))
            for idx,match in enumerate(match_v):
                print('Match ID',idx)
                tpc = tpc_v[match.tpc_id]
                pmt = pmt_v[match.flash_id]
                true_x = 1.e20
                for pt in tpc:
                    if pt.x < true_x: true_x = pt.x
                true_pe = np.sum(pmt.pe_v)
                msg = '  TPC/PMT IDs %d/%d Score %f Min-X %f PE sum %f ... true Min-X %f true PE sum %f'
                msg = msg % (match.tpc_id,match.flash_id,match.score,match.tpc_point.x,np.sum(match.hypothesis),true_x,true_pe)
                print(msg)
            continue
    
        tpc_id_v = [tpc.idx for tpc in tpc_v]
        pmt_id_v = [pmt.idx for pmt in pmt_v]
        all_matches = []
        for match in match_v:
            # FIXME is id same as order in tpc_v? YES
            qcluster = tpc_v[tpc_id_v.index(match.tpc_id)]
            flash = pmt_v[pmt_id_v.index(match.flash_id)]
            store = np.array([[
                event,
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
            #    ('event', 'i4'),
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
            #print(store)
            #print(store.shape)
            all_matches.append(store)
        np_event = np.concatenate(all_matches, axis=0)
        if np_result is None:
            np_result = np_event
        else:
            np_result = np.concatenate([np_result,np_event],axis=0)
    if not out_file:
        return
    #print(x.shape)
    names = [
        'event',
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
    np.savetxt(out_file, np_result, delimiter=',', header=','.join(names))

if __name__ == '__main__':
    import os

    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                            'dat', 'flashmatch.cfg')
    num_tracks = None
    gap = 0.5
    tol = 1e3
    if len(sys.argv) > 1:
        for argv in sys.argv[1:]:
            if argv.isdigit():
                num_tracks = int(argv)
            else:
                cfg_file = sys.argv[1]

    # run demo
    demo(cfg_file,num_tracks=num_tracks)

