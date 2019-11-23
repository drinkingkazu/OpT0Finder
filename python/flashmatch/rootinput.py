import numpy as np
from flashmatch import flashmatch, geoalgo
import sys, ast
from .utils import FlashMatchInput

class ROOTInput:
    def __init__(self, particleana, opflashana, cfg_file=None):
        self.det = flashmatch.DetectorSpecs.GetME()
        self.geoalgo = geoalgo.GeoAlgo()
        self._qcluster_algo  = None
        self._periodTPC = [-1000,1000]
        if not cfg_file is None:
            self.configure(cfg_file)

        from root_numpy import root2array
        self._particles = root2array(particleana, self._tpc_tree_name)
        self._opflash = root2array(opflashana, self._pmt_tree_name)
        # Check data consistency
        self._entries_to_event = np.unique(self._particles['event']).astype(np.int32)

    def event_id(self,entry_id):
        return self._entries_to_event[entry_id]

    def entry_id(self,event_id):
        return np.where(self._entries_to_event == event_id)[0][0]

    def __len__(self):
        return len(self._entries_to_event)

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
        # Set seed if there is any specified
        if pset.contains_value('NumpySeed'):
            seed = int(pset.get('NumpySeed'))
            if seed < 0:
                import time
                seed = int(time.time())
            np.random.seed(seed)

    def num_entries(self):
        return len(self._entries_to_event)

    def get_entry(self,entry):
        event = self._entries_to_event[entry]
        return self.get_event(event)

    def get_event(self,event):
        return self._particles[self._particles['event'] == event], self._opflash[self._opflash['event'] == event]

    def make_qcluster(self,particles,select_pdg=[]):
        qcluster_v = []
        for p_idx, p in enumerate(particles):
            # Only use specified PDG code if pdg_code is available
            if len(select_pdg) and int(p['pdg_code']) not in select_pdg:
                continue
            xyzs = np.column_stack([p['x_v'],p['y_v'],p['z_v']]).astype(np.float64)
            traj = flashmatch.as_geoalgo_trajectory(xyzs)
            # Need at least 2 points to make QCluster
            if traj.size() < 2: continue;
            # Make QCluster
            qcluster = self._qcluster_algo.MakeQCluster(traj)
            # If configured, truncate the physical boundary here
            if self._truncate_tpc_active:
                bbox = self.det.ActiveVolume()
                pt_min, pt_max = bbox.Min(), bbox.Max()
                qcluster.drop(pt_min[0],pt_min[1],pt_min[2],
                              pt_max[0],pt_max[1],pt_max[2])
            # need to be non-empty
            if qcluster.empty(): continue
            ts=p['time_v']
            qcluster.time_true = np.min(ts) * 1.e-3

            # if qcluster.min_x() < -365:
            #     print("touching ", qcluster.min_x(), qcluster.time_true)
            if qcluster.min_x() < -370:
                print('** wrong *** ', qcluster.min_x(), qcluster.max_x(), qcluster.time_true)

            # Assign the index number of a particle
            qcluster.idx = p_idx

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


    def make_flashmatch_input(self, entry):
        """
        Make sample from ROOT files
        """
        if entry > len(self):
            raise IndexError
        result = FlashMatchInput()
        # Find the list of sim::MCTrack entries for this event
        particles, opflash = self.get_entry(entry)
        result.raw_qcluster_v = self.make_qcluster(particles,select_pdg=[13])
        result.qcluster_v = [flashmatch.QCluster_t(tpc) for tpc in result.raw_qcluster_v]
        # If configured, shift X (for MCTrack to immitate reco)
        if self._shift_tpc:
            for i, qcluster in enumerate(result.qcluster_v):
                qcluster.xshift = qcluster.time_true * self.det.DriftVelocity()
                result.qcluster_v[i] = qcluster + qcluster.xshift
        if self._truncate_tpc_readout:
            # Define allowed X recording regions
            min_tpcx, max_tpcx = [t * self.det.DriftVelocity() for t in self._periodTPC]
            for tpc in result.qcluster_v: tpc.drop(min_tpcx,max_tpcx)

        # Find the list of recob::OpFlash entries for this event
        result.flash_v = self.make_flash(opflash)

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
        for pmt in result.flash_v:
            for tpc in result.qcluster_v:
                if tpc.idx in tpc_matched: continue
                dt = abs(pmt.time_true - tpc.time_true)
                if dt < self._matching_window:
                    result.true_match.append((pmt.idx,tpc.idx))
                    tpc_matched.append(tpc.idx)
                    break

        return result
