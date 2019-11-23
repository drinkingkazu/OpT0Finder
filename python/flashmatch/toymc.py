import numpy as np
from flashmatch import flashmatch, geoalgo
import sys, ast
from .utils import FlashMatchInput

class ToyMC:
    def __init__(self, cfg_file=None):
        self.det = flashmatch.DetectorSpecs.GetME()
        self._qcluster_algo  = None
        self._flash_algo     = None
        self._time_algo      = None
        self._track_algo     = None
        self._ly_variation = 0.0
        self._pe_variation = 0.0
        self._periodTPC = [-1000,1000]
        self._periodPMT = [-1000,1000]
        if not cfg_file is None:
            self.configure(cfg_file)

    def event_id(self,entry):
        return -1

    def configure(self,cfg_file):
        self.cfg = flashmatch.CreatePSetFromFile(cfg_file)
        # configure
        #self.mgr.Configure(self.cfg)
        # PhotonLibHypothesis
        #self.mgr.GetAlgo(flashmatch.kFlashHypothesis)
        self._flash_algo = flashmatch.FlashHypothesisFactory.get().create("PhotonLibHypothesis","ToyMCHypothesis")
        self._flash_algo.Configure(self.cfg.get('flashmatch::PSet')("PhotonLibHypothesis"))
        #
        # Toy MC configuration
        #
        pset = self.cfg.get('flashmatch::PSet')('ToyMC')
        # LightPath
        #self._qcluster_algo = self.mgr.GetCustomAlgo(pset.get('QClusterAlgo'))
        self._qcluster_algo = flashmatch.CustomAlgoFactory.get().create("LightPath","ToyMCLightPath")
        self._qcluster_algo.Configure(self.cfg.get('flashmatch::PSet')("LightPath"))
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
        # Period in micro-seconds
        self._periodTPC = ast.literal_eval(pset.get('PeriodTPC'))
        self._periodPMT = ast.literal_eval(pset.get('PeriodPMT'))
        # Truncate TPC tracks (readout effect)
        self._truncate_tpc = int(pset.get("TruncateTPC"))
        # Set seed if there is any specified
        if pset.contains_value('NumpySeed'):
            seed = int(pset.get('NumpySeed'))
            if seed < 0:
                import time
                seed = int(time.time())
            np.random.seed(seed)


    # Create a set ot TPC (flashmatch::QCluster) and PMT (flashmatch::Flash_t) interactions
    def make_flashmatch_input(self,num_match=None):
        """
        Make N input samples for flash matching
        ---------
        Arguments
          num_match: int, number of tpc-pmt pair samples to be generated (optional, default self._num_tracks)
        -------
        Returns
          a list of geoalgo::Trajectory, flashmatch::QClusterArray_t, and flashmatch::FLashArray_t
        """
        result = FlashMatchInput()

        if num_match is None:
            num_match = self._num_tracks
        # Generate 3D trajectories inside the detector
        track_v = self.gen_trajectories(num_match)
        # Generate flash time and x shift (for reco x position assuming trigger time)
        xt_v  = self.gen_xt_shift(len(track_v))
        # Define allowed X recording regions
        min_tpcx, max_tpcx = [t * self.det.DriftVelocity() for t in self._periodTPC]
        pmt_v = flashmatch.FlashArray_t()
        tpc_v = flashmatch.QClusterArray_t()
        raw_tpc_v = flashmatch.QClusterArray_t()
        orphan_pmt=[]
        orphan_tpc=[]
        orphan_raw=[]
        match_id = 0
        true_track_v = []
        for idx, track in enumerate(track_v):
            # TPC xyz & light info
            raw_qcluster = self.make_qcluster(track)
            raw_qcluster.idx = idx
            # PMT pe spectrum
            flash = self.make_flash(raw_qcluster)
            flash.idx = idx
            # Apply x shift and set flash time
            ftime,dx = xt_v[idx]
            flash.time = ftime
            flash.true_time = ftime
            qcluster = raw_qcluster + dx
            qcluster.idx = idx
            qcluster.true_time = ftime
            raw_qcluster.true_time = ftime
            # Drop QCluster points that are outside the recording range
            if self._truncate_tpc:
                qcluster.drop(min_tpcx,max_tpcx)
            # check if this is an orphan
            valid_match = qcluster.size() > 0 and np.sum(flash.pe_v) > 0
            if qcluster.size() >0:
                result.qcluster_v.append(qcluster)
                result.raw_qcluster_v.append(raw_qcluster)
            if np.sum(flash.pe_v) > 0:
                result.flash_v.append(flash)
            if valid_match:
                result.true_match.append((idx,idx))

        return result

    def gen_xt_shift(self,n):
        """
        Generate flash timing and corresponding X shift
        ---------
        Arguments
          n: int, number of track/flash (number of flash time to be generated)
        -------
        Returns
          a list of pairs, (flash time, dx to be applied on TPC points)
        """
        time_dx_v = []
        duration = self._periodPMT[1] - self._periodPMT[0]
        for idx in range(n):
            t,x=0.,0.
            if self._time_algo == 'random':
                t = np.random.random() * duration + self._periodPMT[0]
            elif self._time_algo == 'periodic':
                t = (idx + 0.5) * duration / n + self._periodPMT[0]
            elif self._time_algo == 'same':
                t = 0.
            else:
                raise Exception("Time algo not recognized, must be one of ['random', 'periodic']")
            x = t * self.det.DriftVelocity()
            time_dx_v.append((t,x))
        return time_dx_v


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

        return flash
