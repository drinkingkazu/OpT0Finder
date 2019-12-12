import numpy as np
from flashmatch import flashmatch, geoalgo
import sys, ast
from utils import FlashMatchInput
from scipy.spatial.distance import cdist

class ROOTInput:
    def __init__(self, opflashana, particleana, cfg_file=None):
        self.det = flashmatch.DetectorSpecs.GetME()
        self.geoalgo = geoalgo.GeoAlgo()
        self._qcluster_algo  = None
        self._periodTPC = [-1000,1000]
        if not cfg_file is None:
            self.configure(cfg_file)

        from root_numpy import root2array
        self._particles = root2array(particleana, self._tpc_tree_name)
        self._opflash = root2array(opflashana, self._pmt_tree_name)

    def configure(self,cfg_file):
        self.cfg = flashmatch.CreatePSetFromFile(cfg_file)
        #
        # Toy MC configuration
        #
        pset = self.cfg.get('flashmatch::PSet')('ROOTInput')
        # LightPath
        #self._qcluster_algo = self.mgr.GetCustomAlgo(pset.get('QClusterAlgo'))
        self._qcluster_algo = flashmatch.CustomAlgoFactory.get().create("LightPath","ToyMCLightPath")
        self._qcluster_algo.Configure(self.cfg.get('flashmatch::PSet')("LightPath"))
        # Period in micro-seconds
        self._periodTPC = ast.literal_eval(pset.get('PeriodTPC'))
        # Truncate TPC tracks (readout effect)
        self._truncate_tpc_readout = int(pset.get("TruncateTPCReadout"))
        # Truncate TPC tracks with active volume
        self._truncate_tpc_active  = int(pset.get("TruncateTPCActive"))
        # Shift TPC tracks (for MCTrack to look realistic, needs true timings)
        self._shift_tpc = int(pset.get("ShiftXForMC"))
        # Time window to match MCFlash and MCTrack
        self._matching_window = float(pset.get("MatchingWindow"))
        # Whether to exclude flashes too close to each other
        self._exclude_reflashing = int(pset.get("ExcludeReflashing"))
        self._tpc_tree_name = str(pset.get("TPCTreeName"))
        self._pmt_tree_name = str(pset.get("PMTTreeName"))

        self._clustering = int(pset.get("Clustering"))
        self._clustering_threshold = float(pset.get("ClusteringThreshold"))
        self._clustering_time_window = float(pset.get("ClusteringTimeWindow"))
        self._matching_window_opflash = ast.literal_eval(pset.get("MatchingWindowOpflash"))

        # Set seed if there is any specified
        if pset.contains_value('NumpySeed'):
            seed = int(pset.get('NumpySeed'))
            if seed < 0:
                import time
                seed = int(time.time())
            np.random.seed(seed)

    def get_particle(self,event):
        return self._particles[self._particles['event'] == event]

    def get_opflash(self,event):
        return self._opflash[self._opflash['event'] == event]

    def make_qcluster(self,particles,select_pdg=[]):
        p_idx_v = []
        xyzs_v = []
        ts_v = []
        for p_idx, p in enumerate(particles):
            # Only use specified PDG code if pdg_code is available
            if not self._clustering and len(select_pdg) and int(p['pdg_code']) not in select_pdg:
                continue
            xyzs = np.column_stack([p['x_v'],p['y_v'],p['z_v']]).astype(np.float64)
            if xyzs.shape[0] < 2:
                continue
            merged = False
            if self._clustering:
                # Check if there alreay exists a cluster overlapping in time
                for i, xyzs2 in enumerate(xyzs_v):
                    t2_min = np.min(ts_v[i])
                    t2_max = np.max(ts_v[i])
                    t_min = np.min(p['time_v'])
                    t_max = np.max(p['time_v'])
                    d = -1
                    # if t2_min >= t_min and t2_min <= t_max + self._clustering_time_window:
                    #     # Find whether qclusters are actually close in space
                    #     d = np.min(cdist([xyzs2[0]], xyzs))
                    #
                    # elif t_min >= t2_min and t_min <= t2_max + self._clustering_time_window:
                    #     d = np.min(cdist(xyzs2, [xyzs[0]]))
                    #
                    # if d >= 0 and d < self._clustering_threshold:
                    if (t2_min >= t_min and t2_min <= max(t_max, t_min + self._clustering_time_window)) or (t_min >= t2_min and t_min <= max(t2_max, t2_min + self._clustering_time_window)) :
                        # Merge clusters
                        #print('Merging', p['pdg_code'])
                        new_ts = np.hstack([ts_v[i], p['time_v']])
                        idx = np.argsort(new_ts)
                        ts_v[i] = new_ts[idx]
                        xyzs_v[i] = np.vstack([xyzs_v[i], xyzs])[idx]
                        p_idx_v[i] = min(p_idx_v[i], p_idx)  # FIXME do we want this?
                        merged = True
                        break
            if not merged:
                xyzs_v.append(xyzs)
                ts_v.append(p['time_v'])
                p_idx_v.append(p_idx)

        # Now looping over xyzs and making QClusters
        qcluster_v = []
        #print('Making qclusters from ', len(xyzs_v))
        for i in range(len(xyzs_v)):
            traj = flashmatch.as_geoalgo_trajectory(xyzs_v[i])
            # If configured, truncate the physical boundary here
            # (before shifting or readout truncation)
            if self._truncate_tpc_active:
                bbox = self.det.ActiveVolume()
                traj = self.geoalgo.BoxOverlap(bbox,traj)
            # Need at least 2 points to make QCluster
            if traj.size() < 2: continue;
            # Make QCluster
            qcluster = self._qcluster_algo.MakeQCluster(traj)

            qcluster.time_true = np.min(ts_v[i]) * 1.e-3

            # Assign the index number of a particle
            qcluster.idx = p_idx_v[i]
            qcluster_v.append(qcluster)
        return qcluster_v


    def make_flash(self,opflash):
        flash_v = []
        for f_idx, f in enumerate(opflash):
            flash = flashmatch.Flash_t()
            flash.idx = f_idx
            pe_reco_v = f['pe_v']
            pe_true_v = f['pe_true_v']
            flash.pe_v.resize(self.det.NOpDets(),0.)
            flash.pe_err_v.resize(self.det.NOpDets(),0.)
            flash.pe_true_v.resize(self.det.NOpDets(),0.)
            for pmt in range(self.det.NOpDets()):
                flash.pe_v[pmt] = pe_reco_v[pmt]
                flash.pe_err_v[pmt] = 0.
                flash.time = f['time']
                flash.pe_true_v[pmt] = pe_true_v[pmt]
                flash.time_true = f['time_true']
            if np.sum(flash.pe_v) > 0:
                flash_v.append(flash)
        return flash_v


    def make_flashmatch_input(self, event):
        """
        Make sample from ROOT files
        """
        result = FlashMatchInput()
        # Find the list of sim::MCTrack entries for this event
        particles = self.get_particle(event)
        result.raw_qcluster_v = self.make_qcluster(particles,select_pdg=[13])
        result.qcluster_v = [flashmatch.QCluster_t(tpc) for tpc in result.raw_qcluster_v]
        # If configured, shift X (for MCTrack to imitate reco)
        if self._shift_tpc:
            for i, qcluster in enumerate(result.qcluster_v):
                qcluster.xshift = qcluster.time_true * self.det.DriftVelocity()
                if qcluster.min_x() - self.det.ActiveVolume().Min()[0] < self.det.ActiveVolume().Max()[0] - qcluster.max_x():
                    result.qcluster_v[i] = qcluster + qcluster.xshift #+ (qcluster.min_x() - self.det.ActiveVolume().Min()[0])
                else:
                    result.qcluster_v[i] = qcluster - qcluster.xshift #- (self.det.ActiveVolume().Max()[0] - qcluster.min_x())

        if self._truncate_tpc_readout:
            # Define allowed X recording regions
            min_tpcx, max_tpcx = [t * self.det.DriftVelocity() for t in self._periodTPC]
            for tpc in result.qcluster_v: tpc.drop(min_tpcx,max_tpcx)

        # Find the list of recob::OpFlash entries for this event
        result.flash_v = self.make_flash(self.get_opflash(event))

        # compute dt to previous and next flash
        for pmt in result.flash_v:
            pmt.dt_prev,pmt.dt_next = flashmatch.kINVALID_DOUBLE,flashmatch.kINVALID_DOUBLE
        for pmt in result.flash_v:
            for pmt2 in result.flash_v:
                if pmt2.idx == pmt.idx: continue
                if pmt2.time < pmt.time and abs(pmt2.time - pmt.time) < pmt.dt_prev:
                    pmt.dt_prev = abs(pmt2.time - pmt.time)
                if pmt2.time > pmt.time and abs(pmt.time - pmt2.time) < pmt.dt_next:
                    pmt.dt_next = abs(pmt.time - pmt2.time)
            #print(pmt.time, pmt.dt_prev, pmt.dt_next)

        # Exclude flashes too close apart
        if self._exclude_reflashing:
            selected = []
            for pmt in result.flash_v:
                if pmt.dt_prev > dt_threshold and pmt.dt_next > dt_threshold:
                    selected.append(pmt)
                else:
                    print('dropping', pmt.dt_prev, pmt.dt_next)
            result.flash_v = selected

        # Assign idx now based on true timings
        tpc_matched = []
        pmt_matched = []
        for pmt in result.flash_v:
            for tpc in result.qcluster_v:
                if tpc.idx in tpc_matched: continue
                dt = abs(pmt.time_true - tpc.time_true)
                if dt < self._matching_window:
                    result.true_match.append((pmt.idx,tpc.idx))
                    tpc_matched.append(tpc.idx)
                    pmt_matched.append(pmt.idx)
                    break
        # Assign idx based on opflash timings for tpc that have not been matched
        pmt_matched_second = []
        for tpc in result.qcluster_v:
            if tpc.idx in tpc_matched: continue
            for pmt in result.flash_v:
                if pmt.idx in pmt_matched or pmt.idx in pmt_matched_second: continue
                if pmt.time_true > 1e5: # this opflash has no mcflash
                    dt = (pmt.time - tpc.time_true)
                    if dt > self._matching_window_opflash[0] and dt < self._matching_window_opflash[1]:
                        result.true_match.append((pmt.idx, tpc.idx))
                        pmt_matched_second.append(pmt.idx)
        return result
