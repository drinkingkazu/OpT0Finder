import numpy as np
from flashmatch import flashmatch, geoalgo
import sys, ast

class ToyMC:
    def __init__(self, cfg_file=None, opflashana=None,particleana=None):
        self.mgr = flashmatch.FlashMatchManager()
        self.det = flashmatch.DetectorSpecs.GetME()
        self._qcluster_algo  = None
        self._flash_algo     = None
        self._time_algo      = None
        self._track_algo     = None
        self._ly_variation = 0.0
        self._pe_variation = 0.0
        self._periodTPC = [-1000,1000]
        self._periodPMT = [-1000,1000]
        self._min_flash_pe = 0.0
        self._min_track_length = 0.0
        if not cfg_file is None:
            self.configure(cfg_file)
        if particleana is not None and opflashana is not None:
            from root_numpy import root2array
            # import glob
            # particleana_list = glob.glob(particleana)
            # opflashana_list = glob.glob(opflashana)
            # assert len(particleana_list) == list(opflashana_list)
            # #particles = root2array(particleana, 'largeant_particletree')
            # for i in range(len(particleana_list)):
            self._particles = root2array(particleana, 'mctracktree')
            self._opflash = root2array(opflashana, self._tree_name)
            #self._opflash = root2array(opflashana, 'cheatflash_flashtree')

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
        # Period in micro-seconds
        self._periodTPC = ast.literal_eval(pset.get('PeriodTPC'))
        self._periodPMT = ast.literal_eval(pset.get('PeriodPMT'))
        # Minimum track length to be registered for matching
        self._min_track_length = float(pset.get("MinTrackLength"))
        # Minimum PE to be registered for matching
        self._min_flash_pe = float(pset.get("MinFlashPE"))
        # Truncate TPC tracks (readout effect)
        self._truncate_tpc = int(pset.get("TruncateTPC"))
        # Time window to match MCFlash and MCTrack
        self._matching_window = float(pset.get("MatchingWindow"))
        # Whether to exclude flashes too close to each other
        self._exclude_reflashing = int(pset.get("ExcludeReflashing"))
        self._tree_name = str(pset.get("TreeName"))
        # Set seed if there is any specified
        if pset.contains_value('NumpySeed'):
            seed = int(pset.get('NumpySeed'))
            if seed < 0:
                import time
                seed = int(time.time())
            np.random.seed(seed)


    def read_input(self, event):
        """
        Make sample from ROOT files
        """
        # from ROOT import TChain
        # ch_particle = TChain('largeant_particletree')
        # ch_particle.AddFile(particleana)
        # ch_opflash = TChain('opflash_flashtree')
        # ch_opflash.AddFile(opflashana)
        # # ch_cheatflash = TChain('cheatflash_flashtree')
        # num_tracks = ch_particle.GetEntries()
        # for e in range(num_tracks):
        #     ch_particle.GetEntry(e)
        #     ch_opflash.GetEntry(e)
        #     br_particle = getattr(ch_particle, 'largeant_particlebranch')
        #     br_opflash = getattr(ch_opflash, 'opflash_flashbranch')


        particles = self._particles[(self._particles['event'] == event) & (self._particles['pdg_code'] == 13)]
        opflash = self._opflash[self._opflash['event'] == event]

        num_tracks = max(len(particles), len(opflash))
        print("Particle ", len(particles), " Flashes ", len(opflash))
        # for p in particles:
        #     print(p['start_x'], p['start_y'], p['start_z'])
        #     print(p['end_x'], p['end_y'], p['end_z'])
        for f in opflash:
            print("Opflash sum = ", f['pe_v'][:180].sum(), "; MCFlash sum = ", f['pe_true_v'][:180].sum())
        #assert len(particles) == len(opflash)
        #xt_v  = self.gen_xt_shift(num_tracks)

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
        print('Loaded %d tracks' % num_tracks)
        print('drift velocity = ', self.det.DriftVelocity())
        track_root_idx, flash_root_idx = -1, -1
        for i in range(num_tracks):
            flash = None
            if i < len(opflash):
                flash_root_idx += 1
                flash = flashmatch.Flash_t()
                #print(len(opflash[i]['pe_v']))
                for x in opflash[i]['pe_v'][:180]:  # FIXME hardcoded NOpDets
                    flash.pe_v.push_back(x)
                    flash.pe_err_v.push_back(0.)
                for x in opflash[i]['pe_true_v'][:180]:
                    flash.pe_true_v.push_back(x)
                flash.x = opflash[i]['x']
                flash.y = opflash[i]['y']
                flash.z = opflash[i]['z']
                flash.xerr, flash.yerr, flash.zerr = 0., 0., 0.
                flash.time = opflash[i]['time']
                flash.time_true =  opflash[i]['time_true']
                flash.idx = i
                flash.ROOT_idx = flash_root_idx

                #print('flash time', opflash[i]['time'], opflash[i]['time_true'])

            #raw_qcluster = flashmatch.QCluster_t()
            #raw_qcluster.time = particles[i]['time_v'][0]
            #raw_qcluster.idx = i
            track, qcluster, raw_qcluster = None, None, None
            if i < len(particles):
                track_root_idx += 1
                track = geoalgo.Trajectory(len(particles[i]['x_v']), 3)
                track_time = 10000000000000000  # will be in ns
                for j in range(len(particles[i]['x_v'])):
                    # print(particles[i]['time_v'][j])
                    # qpoint = flashmatch.QPoint_t(
                    #     particles[i]['x_v'][j],
                    #     particles[i]['y_v'][j],
                    #     particles[i]['z_v'][j],
                    #     particles[i]['energy_v'][j]*24000
                    # )
                    #print(i, j, qpoint.x)
                    #raw_qcluster.push_back(qpoint)
                    track[j] = geoalgo.Vector(
                        particles[i]['x_v'][j],
                        particles[i]['y_v'][j],
                        particles[i]['z_v'][j]
                    )
                    if particles[i]['time_v'][j] < track_time:
                        track_time = particles[i]['time_v'][j]
                if track.size() == 0:
                    track = None
                if track is not None:
                    raw_qcluster = self.make_qcluster(track)
                    raw_qcluster.time_true = track_time*1e-3
                    raw_qcluster.ROOT_idx = track_root_idx
                    qcluster = raw_qcluster
                    # Drop QCluster points that are outside the recording range
                    if self._truncate_tpc:
                        qcluster.drop(min_tpcx,max_tpcx)

            # check if this is an orphan
            orphan = qcluster is not None and flash is not None and qcluster.size() and track.Length() > self._min_track_length and np.sum(flash.pe_v) > self._min_flash_pe
            #print(i, orphan, qcluster, flash)
            if orphan:
                # set IDs
                qcluster.idx = match_id
                flash.idx    = match_id
                match_id += 1
                pmt_v.push_back(flash)
                tpc_v.push_back(qcluster)
                raw_tpc_v.push_back(raw_qcluster)
                true_track_v.append(track)
            else:
                if flash is not None and np.sum(flash.pe_v) > self._min_flash_pe:
                    orphan_pmt.append(flash)
                if qcluster is not None and qcluster.size() and track is not None and track.Length() > self._min_track_length:
                    orphan_tpc.append(qcluster)
                    orphan_raw.append(raw_qcluster)

        # append orphans
        for orphan in orphan_pmt:
            orphan.idx = pmt_v.size()
            pmt_v.push_back(orphan)
        for orphan in orphan_tpc:
            orphan.idx = tpc_v.size()
            tpc_v.push_back(orphan)
        for orphan in orphan_raw:
            orphan.idx = raw_tpc_v.size()
            raw_tpc_v.push_back(orphan)

        # compute dt to previous and next flash
        for pmt in pmt_v:
            for pmt2 in pmt_v:
                if pmt2.idx != pmt.idx:
                    if pmt2.time < pmt.time and abs(pmt2.time - pmt.time) < pmt.dt_prev:
                        pmt.dt_prev = abs(pmt2.time - pmt.time)
                    if pmt2.time > pmt.time and abs(pmt.time - pmt2.time) < pmt.dt_next:
                        pmt.dt_next = abs(pmt.time - pmt2.time)
            #print(pmt.time, pmt.dt_prev, pmt.dt_next)

        # Exclude flashes too close apart
        if self._exclude_reflashing:
            new_pmt_v = flashmatch.FlashArray_t()
            dt_threshold = 15 # us
            for pmt in pmt_v:
                if pmt.dt_prev > dt_threshold and pmt.dt_next > dt_threshold:
                    new_pmt_v.push_back(pmt)
                else:
                    print('dropping', pmt.dt_prev, pmt.dt_next)
            pmt_v = new_pmt_v

        # Assign idx now based on true timings
        for pmt in pmt_v:
            pmt.idx = 0
        for tpc in tpc_v:
            tpc.idx = 0
        global_idx = 1
        time_threshold = self._matching_window # 100ns
        for pmt in pmt_v:
            for tpc in tpc_v:
                if tpc.idx == 0 and abs(pmt.time_true - tpc.time_true) < time_threshold:
                    pmt.idx = global_idx
                    tpc.idx = global_idx
                    print(global_idx, pmt.time_true, tpc.time_true, tpc.min_x())
                    global_idx += 1

        for tpc in tpc_v:
            if tpc.idx == 0:
                tpc.idx = global_idx
                global_idx += 1
        for pmt in pmt_v:
            if pmt.idx == 0:
                pmt.idx = global_idx
                global_idx += 1
        for i, tpc in enumerate(raw_tpc_v):
            tpc.idx = tpc_v[i].idx

        return true_track_v, pmt_v, tpc_v, raw_tpc_v

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
            # PMT pe spectrum
            flash = self.make_flash(raw_qcluster)
            # Apply x shift and set flash time
            ftime,dx = xt_v[idx]
            flash.time = ftime
            qcluster = raw_qcluster + dx
            # Drop QCluster points that are outside the recording range
            if self._truncate_tpc:
                qcluster.drop(min_tpcx,max_tpcx)
            # check if this is an orphan
            orphan = qcluster.size() and track.Length() > self._min_track_length and np.sum(flash.pe_v) > self._min_flash_pe
            if not orphan:
                # set IDs
                qcluster.idx = match_id
                flash.idx    = match_id
                match_id += 1
                pmt_v.push_back(flash)
                tpc_v.push_back(qcluster)
                raw_tpc_v.push_back(raw_qcluster)
                true_track_v.append(track)
            else:
                if np.sum(flash.pe_v) > self._min_flash_pe:
                    orphan_pmt.append(flash)
                if qcluster.size() and track.Length() > self._min_track_length:
                    orphan_tpc.append(qcluster)
                    orphan_raw.append(raw_qcluster)

        # append orphans
        for orphan in orphan_pmt:
            orphan.idx = pmt_v.size()
            pmt_v.push_back(orphan)
        for orphan in orphan_tpc:
            orphan.idx = tpc_v.size()
            tpc_v.push_back(orphan)
        for orphan in orphan_raw:
            orphan.idx = raw_tpc_v.size()
            raw_tpc_v.push_back(orphan)

        return true_track_v, pmt_v, tpc_v, raw_tpc_v

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

def demo(cfg_file,repeat=17,num_tracks=None,out_file='',particleana=None,opflashana=None):
    """
    Run function for ToyMC
    ---------
    Arguments
      cfg_file:   string for the location of a config file
      repeat:     int, number of times to run the toy MC simulation
      out_file:   string for an output analysis csv file path (optional)
      num_tracks: int for number of tracks to be generated (optional)
    """
    mgr = ToyMC(cfg_file, particleana=particleana, opflashana=opflashana)
    # dump config
    sys.stdout.write(mgr.cfg.dump())
    sys.stdout.flush()
    # override number of tracks to simulate
    event_list = range(repeat)
    if num_tracks is not None:
        num_tracks = int(num_tracks)
    if particleana is not None and opflashana is not None:
        event_list = np.unique(mgr._particles['event'])
        print('Found %d events' % len(event_list))

    np_result = None
    for event in event_list:
        sys.stdout.write('Event %d/%d\n' %(event,len(event_list)))
        sys.stdout.flush()
        # Generate samples
        if particleana is None or opflashana is None:
            track_v, pmt_v, tpc_v, raw_tpc_v = mgr.gen_input(num_tracks)
        else:
            track_v, pmt_v, tpc_v, raw_tpc_v = mgr.read_input(event)
        # Run matching
        match_v = mgr.match(tpc_v, pmt_v)

        #tpc_v = [flashmatch.as_ndarray(tpc) for tpc in tpc_v]
        #pmt_v = [flashmatch.as_ndarray(pmt) for pmt in pmt_v]

        # Either report or store output
        # true_xmin_v = []
        # for idx in range(tpc_v.size()):
        #     tpc = tpc_v[idx]
        #     pmt = pmt_v[idx]
        #     true_xmin_v.append(tpc.min_x() - pmt.time * mgr.det.DriftVelocity())

        # If no analysis output saving option given, return
        if not out_file:
            print('Number of match result',len(match_v))

        all_matches = []
        for idx, match in enumerate(match_v):
            # FIXME is id same as order in tpc_v? YES
            tpc = tpc_v[match.tpc_id]
            pmt = pmt_v[match.flash_id]
            raw = raw_tpc_v[match.tpc_id]
            truncation   = raw.front().dist(raw.back()) - tpc.front().dist(tpc.back())
            if raw.front().dist(raw.back()) > 0:
                truncation_frac = truncation / raw.front().dist(raw.back())
            else:
                truncation_frac = -1
            #if particleana is None or opflashana is None:
            true_minx = tpc.min_x() #- pmt.time * mgr.det.DriftVelocity() #true_xmin_v[match.flash_id]
            true_maxx = tpc.max_x()
            reco_maxx = match.tpc_point.x + (tpc.max_x() - tpc.min_x())
            correct_match = tpc.idx == pmt.idx and tpc.idx < len(track_v)

            if not out_file:
                print('Match ID',idx)
                msg = '  TPC/PMT IDs %d/%d Score %f Min-X %f PE sum %f ... true Min-X %f true PE sum %f truncation %f (%f%%)'
                msg = msg % (tpc.idx,pmt.idx,match.score,match.tpc_point.x,np.sum(match.hypothesis),
                             true_minx,np.sum(pmt.pe_v),truncation,truncation_frac*100.)
                print(msg)
                continue
            store = np.array([[
                event,
                match.score,
                pmt.ROOT_idx,
                tpc.ROOT_idx,
                true_minx,
                true_maxx,
                reco_maxx,
                match.tpc_point.x,
                match.tpc_point.y,
                match.tpc_point.z,
                tpc.front().x,
                tpc.front().y,
                tpc.front().z,
                tpc.back().x,
                tpc.back().y,
                tpc.back().z,
                raw.front().x,
                raw.front().y,
                raw.front().z,
                raw.back().x,
                raw.back().y,
                raw.back().z,
                truncation,
                truncation_frac,
                len(tpc),
                int(tpc.idx == pmt.idx),
                tpc.sum(),
                tpc.length(),
                tpc.time_true,
                np.sum(match.hypothesis),
                pmt.TotalPE(),
                pmt.TotalTruePE(),
                pmt.time,
                pmt.time_true,
                pmt.dt_prev,
                pmt.dt_next,
                match.duration,
                match.num_steps,
                match.minimizer_min_x,
                match.minimizer_max_x
            ]])
            #)]], dtype=[
            #    ('event', 'i4'),
            #    ('score', 'f4'),
            #    ('true_min_x', 'f4'),
            #    ('tpc_point_x', 'f4'),
            #    ('tpc_point_y', 'f4'),
            #    ('tpc_point_z', 'f4'),
            #    ('start_point_x', 'f4'),
            #    ('start_point_y', 'f4'),
            #    ('start_point_z', 'f4'),
            #    ('end_point_x', 'f4'),
            #    ('end_point_y', 'f4'),
            #    ('end_point_z', 'f4'),
            #    ('raw_start_point_x', 'f4'),
            #    ('raw_start_point_y', 'f4'),
            #    ('raw_start_point_z', 'f4'),
            #    ('raw_end_point_x', 'f4'),
            #    ('raw_end_point_y', 'f4'),
            #    ('raw_end_point_z', 'f4'),
            #    ('truncation', 'f4'),
            #    ('truncation_fraction', 'f4'),
            #    ('qcluster_num_points', 'f4'),
            #    ('matched', 'B'),
            #    ('qcluster_sum', 'f4'),
            #    ('flash_sum', 'f4'),
            #    ('flash_time', 'f4'),
            #    ('duration', 'f4')
            #])
            #print(store)
            #print(store.shape)
            all_matches.append(store)
        if out_file and len(all_matches):
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
        'flash_ROOT_idx',
        'track_ROOT_idx',
        'true_min_x',
        'true_max_x',
        'reco_max_x',
        'tpc_point_x',
        'tpc_point_y',
        'tpc_point_z',
        'start_point_x',
        'start_point_y',
        'start_point_z',
        'end_point_x',
        'end_point_y',
        'end_point_z',
        'raw_start_point_x',
        'raw_start_point_y',
        'raw_start_point_z',
        'raw_end_point_x',
        'raw_end_point_y',
        'raw_end_point_z',
        'truncation',
        'truncation_fraction',
        'qcluster_num_points',
        'matched',
        'qcluster_sum',
        'qcluster_length',
        'qcluster_time_true',
        'hypothesis_sum',  # Hypothesis flash sum
        'flash_sum', # OpFlash Sum
        'flash_true_sum', # MCFlash sum
        'flash_time',
        'flash_time_true',
        'flash_dt_prev',
        'flash_dt_next',
        'duration',
        'num_steps',
        'minimizer_min_x',
        'minimizer_max_x'
    ]
    np.savetxt(out_file, np_result, delimiter=',', header=','.join(names))

if __name__ == '__main__':
    import os

    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                            'dat', 'flashmatch.cfg')
    num_tracks = None
    particleana, opflashana = None, None
    outfile = ''
    if len(sys.argv) > 1:
        if len(sys.argv) == 2:
            num_tracks = int(sys.argv[1])
        else:
            particleana = sys.argv[1]
            opflashana = sys.argv[2]
            print(particleana, opflashana)
            if len(sys.argv) > 3:
                outfile = sys.argv[3]
        # for argv in sys.argv[1:]:
        #     if argv.isdigit():
        #         num_tracks = int(argv)
        #     else:
        #         cfg_file = sys.argv[1]

    # run demo
    if num_tracks is not None:
        demo(cfg_file,num_tracks=num_tracks)
    else:
        print('particleana', particleana)
        print('opflashana', opflashana)
        print('outfile', outfile)
        demo(cfg_file, particleana=particleana, opflashana=opflashana, out_file=outfile)
