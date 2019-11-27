import sys
import numpy as np
import plotly.graph_objs as go
from flashmatch import flashmatch, AnalysisManager
from .vis_icarus import vis_icarus, icarus_layout3d

class DataManager:
    """
    DataManager class "holds" data retrieved using AnalysisManager. In particular, it holds:
      - "current" entry and event id
      - "cpp" QCluster,Flash_t, etc. in flashmatch.FlashMatchInput format (this is self.cpp attribute)
      - numpy array representation of what's stored in "cpp"
      - Read & update local data attributes is done via self.update function
    TODO: useful properties of DataManager is incremental to AnalysisManager. One can
          either merge two classes or make DataManager inherit from AnalysisManager, if wished.
    """
    def __init__(self):
        self.entry = -1
        self.event = -1
        self.hypothesis_made = False

    def update(self,manager,entry,is_entry=True,make_hypothesis=False):
        """
        Given flashmatch.AnalysisManager, which can interface input data stream and OpT0Finder tools,
        this function "read and update" local data attributes if needed.
        """
        if not is_entry:
            try:
                entry = manager.entry_id(entry)
            except IndexError:
                return False
        if entry < 0 or entry >= len(manager.entries()):
            return False
        if not self.entry == entry:
            self.cpp = manager.make_flashmatch_input(entry)
            self.np_qcluster_v     = [flashmatch.as_ndarray(qcluster) for qcluster in self.cpp.qcluster_v    ]
            self.np_raw_qcluster_v = [flashmatch.as_ndarray(qcluster) for qcluster in self.cpp.raw_qcluster_v]
            self.np_all_pts_v      = [flashmatch.as_ndarray(qcluster) for qcluster in self.cpp.all_pts_v     ]
            self.np_flash_v        = [flashmatch.as_ndarray(flash)    for flash    in self.cpp.flash_v       ]
            pmt_maxid,tpc_maxid=0,0
            for tpc in self.cpp.qcluster_v: tpc_maxid = max(tpc_maxid,tpc.idx)
            for pmt in self.cpp.flash_v: pmt_maxid = max(pmt_maxid,pmt.idx)
            self.np_match_tpc2pmt  = np.ones(shape=(tpc_maxid+1),dtype=np.int32) * -1
            self.np_match_pmt2tpc  = np.ones(shape=(pmt_maxid+1),dtype=np.int32) * -1
            for match in self.cpp.true_match:
                self.np_match_tpc2pmt[match[1]]=match[0]
                self.np_match_pmt2tpc[match[0]]=match[1]
            self.entry = entry
            self.event = manager.event_id(entry)
            self.hypothesis_made = False
        if make_hypothesis and not self.hypothesis_made:
            self.cpp.hypo_v = [manager.flash_hypothesis(qcluster) for qcluster in self.cpp.qcluster_v]
            self.np_hypo_v  = [flashmatch.as_ndarray(flash)       for flash    in self.cpp.hypo_v    ]
            for idx in range(len(self.np_hypo_v)):
                self.cpp.hypo_v[idx].idx = self.cpp.qcluster_v[idx].idx
        return True


class AppManager:
    """
    AppManager serves 2 purposes: 0) it holds data objects which state must be tracked
    through the lifetime of an (stateless) dash html app, and 1) it implements useful
    data parsing methods (e.g. creating a figure, update dropdown options, etc.).
    """
    def __init__(self,cfg,geo,data_particle,data_opflash,dark_mode):
        """
        INPUT:
          - cfg ... OpT0Finder configuration file, passed onto AnalysisManager
          - geo ... vis_icarus (geometry) data file in either PSet or yaml format
          - data_particle ... particle data file stored by ICARUSParticleAna_module
          - data_opflash ... OpFlash data file stored by ICARUSOpFlashAna_module
          - dark_mode ... plotly_dark to be consistent with figure/layout
        """
        self.dark_mode = bool(dark_mode)
        self.ana_manager = AnalysisManager(cfg=cfg, particleana=data_particle, opflashana=data_opflash)
        assert(len(self.ana_manager.entries()))
        self.dat_manager = DataManager()
        self.vis_icarus = vis_icarus(geo)
        self.detector_trace = self.vis_icarus.get_trace_detector()
        self.detector_trace_no_pmt = self.vis_icarus.get_trace_detector(draw_pmts=False)
        self.layout = icarus_layout3d(self.vis_icarus.data(),set_camera=False,dark=self.dark_mode)
        self.empty_view = go.Figure(self.detector_trace,layout=self.layout)
        sys.stdout.flush()


    def qll_score(self,hypothesis,flash):
        res = self.ana_manager.matcher.GetAlgo(flashmatch.kFlashMatch).QLL(hypothesis,flash)
        sys.stdout.flush()
        return np.power(10,-1.*res)

    def run_flashmatch(self,qcluster,flash):
        self.ana_manager.matcher.Reset()
        self.ana_manager.matcher.Add(qcluster)
        self.ana_manager.matcher.Add(flash)
        res = self.ana_manager.matcher.Match()
        sys.stdout.flush()
        if not len(res) == 1: return None
        else: return res[0]

    def valid_data_entry(self,entry,is_entry):
        """
        Returns data entry given a valid data index (exist in file) or None
        INPUTS:
          - entry ... either data index ... entry or event id
          - is_entry ... False if the input data_index is event id
        """
        if not is_entry:
            try:
                entry = self.ana_manager.entry_id(entry)
                return entry
            except IndexError:
                return None
        if entry < 0 or entry >= len(self.ana_manager.entries()):
            return None
        return entry

    def current_data_index(self,is_entry):
        """
        Returns current "data index"
        INPUTS:
          - is_entry ... a boolean True=return entry, False=return event id
        OUTPUT:
          - data index integer, either entry or event id
        """
        return self.dat_manager.entry if is_entry else self.dat_manager.event

    def dropdown_qcluster(self,data_index,is_entry):
        """
        Generates a dropdown menu for a dash app to select a QCluster
        INPUTS:
          - data_index ... an integer to specify which entry/event to read
          - is_entry ... True=entry, False=event id
        OUTPUT:
          A list of dicts to be consumed by dash app dropdown options
        """
        if not self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry):
            return None
        label_v = []
        for idx, qc in enumerate(self.dat_manager.cpp.qcluster_v):
            label = 'Track %02d (%dpts Match=%d T=%.3fus)'
            label = label % (qc.idx, qc.size(), self.dat_manager.np_match_tpc2pmt[qc.idx], qc.time_true)
            label_v.append(label)
        dropdown_qcluster = [dict(label=label_v[idx], value=idx)
                             for idx,qcluster in enumerate(self.dat_manager.np_qcluster_v)]
        dropdown_qcluster += [dict(label='All tracks',value=len(self.dat_manager.np_qcluster_v))]
        return dropdown_qcluster

    def dropdown_flash(self,data_index,is_entry,mode_flash):
        """
        Generates a dropdown menu for a dash app to select a Flash
        INPUTS:
          - data_index ... an integer to specify which entry/event to read
          - is_entry ... True=entry, False=event id
          - mode_flash ... True=Flash, False=Hypothesis
        OUTPUT:
          A list of dicts to be consumed by dash app dropdown options
        """
        if not self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry, make_hypothesis=(not mode_flash)):
            return None
        target_v = self.dat_manager.np_flash_v if mode_flash else self.dat_manager.np_hypo_v
        label_v = []
        if mode_flash:
            for idx, pmt in enumerate(self.dat_manager.cpp.flash_v):
                pesum = self.dat_manager.np_flash_v[idx].sum()
                label = 'Flash %02d (PE=%de2 Match=%d T=%.3fus MC-T=%.2fus)'
                label = label % (pmt.idx, pesum/100.,
                                 self.dat_manager.np_match_pmt2tpc[pmt.idx], pmt.time, pmt.time_true)
                label_v.append(label)
        else:
            for idx, tpc in enumerate(self.dat_manager.cpp.qcluster_v):
                pesum = self.dat_manager.np_hypo_v[idx].sum()
                label = 'Hypothesis %02d (PE=%de2 MC-T=%.3fus)'
                label = label % (tpc.idx, pesum/100., tpc.time_true)
                label_v.append(label)
        if mode_flash: time_v = [flash.time for flash in self.dat_manager.cpp.flash_v]
        dropdown_flash  = [dict(label=label_v[idx],value=idx)
                            for idx,flash in enumerate(target_v)]
        dropdown_flash += [dict(label='All flashes',value=len(target_v))]
        sys.stdout.flush()
        return dropdown_flash

    def event_display(self, data_index, qcluster_idx_v, flash_idx_v,
                      is_entry, mode_flash, mode_qcluster, use_all_pts, pmt_range):
        """
        Generate 3D display for change in event/entry and/or selection of flash/qcluster
        """
        if not self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry, make_hypothesis=(not mode_flash)):
            return None
        data = []
        # QCluster
        # if mode_qcluster == 0, then only show qcluster with color-per-cluster
        # if mode_qcluster == 1, then show both qcluster and raw version with 2 colors (red + green)
        # if mode_qcluster == 2, then only show raw qcluster with color-per-cluster

        if qcluster_idx_v is not None and len(qcluster_idx_v):
            if len(self.dat_manager.np_qcluster_v) in qcluster_idx_v:
                qcluster_idx_v = range(len(self.dat_manager.np_qcluster_v))
            idx_v = [self.dat_manager.cpp.qcluster_v[idx].idx for idx in qcluster_idx_v]
            for idx in qcluster_idx_v:
                if mode_qcluster in [0,1]:
                    xyz = self.dat_manager.np_qcluster_v[idx]
                    name = 'Track %02d (%d pts)' % (self.dat_manager.cpp.qcluster_v[idx].idx,len(xyz))
                    trace = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
                                         name=name,
                                         #template='plotly' if not self.dark_mode else 'plotly_dark',
                                         marker = dict(size=2, opacity=0.5, color=None if mode_qcluster==0 else 'orange')
                                        )
                    data.append(trace)
                if mode_qcluster in [1,2]:
                    xyz = self.dat_manager.np_raw_qcluster_v[idx] if not use_all_pts else self.dat_manager.np_all_pts_v[idx]
                    name = 'Track (Raw) %02d (%d pts)' % (self.dat_manager.cpp.qcluster_v[idx].idx,len(xyz))
                    trace = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
                                         name=name,
                                         #template='plotly' if not self.dark_mode else 'plotly_dark',
                                         marker = dict(size=2, opacity=0.5, color=None if mode_qcluster==2 else 'cyan')
                                        )
                    data.append(trace)
        target_v = self.dat_manager.np_flash_v if mode_flash else self.dat_manager.np_hypo_v
        cpp_target_v = self.dat_manager.cpp.flash_v if mode_flash else self.dat_manager.cpp.hypo_v
        if flash_idx_v and len(flash_idx_v):
            pmt_pos = self.vis_icarus.data()['pmts']
            pmt_val = None
            if len(target_v) in flash_idx_v: pmt_val = np.sum(target_v,axis=0)
            else: pmt_val = np.sum(np.column_stack([target_v[idx] for idx in flash_idx_v]),axis=1)
            pmt_range[0] = np.min(pmt_val) if pmt_range[0].lower() == 'min' else float(pmt_range[0])
            pmt_range[1] = np.max(pmt_val) if pmt_range[1].lower() == 'max' else float(pmt_range[1])
            name = 'Flash (%dE2 PEs)' % int(np.sum(pmt_val)/100.)
            if(len(flash_idx_v)==1 and mode_flash):
                name = 'Flash (%dE2 PEs @ %.2f us)' % (int(np.sum(pmt_val)/100.), int(cpp_target_v[flash_idx_v[0]].time))
            trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                 name=name,
                                 #template='plotly' if not self.dark_mode else 'plotly_dark',
                                 marker = dict(size=6, color=pmt_val, opacity=0.5, cmin=pmt_range[0], cmax=pmt_range[1])
                                )
            data.append(trace)
        elif len(data):
            # no pmt data => show default pmt view
            return go.Figure(self.detector_trace + data, layout=self.layout)

        sys.stdout.flush()
        if len(data)<1:
            return self.empty_view
        else:
            return go.Figure(self.detector_trace_no_pmt + data, layout=self.layout)

    def hypothesis_display(self, data_index, qcluster_idx, flash_idx,
                           pe_range, normalize_pe, is_entry, xpos):
        """
        Generate 3D display of a flash hypothesis for a given QCluster and x-offset
        """
        if not self.dat_manager.update(self.ana_manager, entry=data_index,
                                       is_entry=is_entry, make_hypothesis=False):
            return None
        data = []
        hypo_val,flash_val,diff_val=None,None,None
        hypo_name,flash_name,diff_name=None,None,None
        if qcluster_idx is not None:
            # Find QCluster
            qcluster = self.dat_manager.cpp.qcluster_v[qcluster_idx]
            raw_qcluster = self.dat_manager.cpp.raw_qcluster_v[qcluster_idx]
            # find x span allowed
            xspan = qcluster.max_x() - qcluster.min_x()
            # if xoffset is larger than allowed range, truncate
            active_bb  = flashmatch.DetectorSpecs.GetME().ActiveVolume()
            max_xpos = min(active_bb.Min()[0] + xspan, active_bb.Max()[0])
            xpos     = max(min(xpos,max_xpos),active_bb.Min()[0])
            # shift to xoffset
            offset_qcluster = flashmatch.QCluster_t(qcluster) - qcluster.min_x() + xpos
            # Make hypothesis
            hypothesis = self.ana_manager.flash_hypothesis(offset_qcluster)
            hypo_val = flashmatch.as_ndarray(hypothesis)
            hypo_name ='Hypothesis (%dE2 PEs)' % int(np.sum(hypo_val)/100.)
            # make a trace for tpc
            xyzv = flashmatch.as_ndarray(offset_qcluster)
            tpc_trace = go.Scatter3d(x=xyzv[:,0],y=xyzv[:,1],z=xyzv[:,2],mode='markers',
                                     name='QCluster (X %.2f)' % xpos,
                                     #template='plotly' if not self.dark_mode else 'plotly_dark',
                                     marker = dict(size=2,opacity=0.5,color='orange')
                                    )
            data.append(tpc_trace)
            # also plot raw all points
            xyzv = self.dat_manager.np_all_pts_v[qcluster_idx]
            tpc_trace = go.Scatter3d(x=xyzv[:,0],y=xyzv[:,1],z=xyzv[:,2],mode='markers',
                                     name='True all points',
                                     #template='plotly' if not self.dark_mode else 'plotly_dark',
                                     marker = dict(size=2,opacity=0.5,color='green')
                                    )
            data.append(tpc_trace)

        if flash_idx is not None:
            pmt_pos = self.vis_icarus.data()['pmts']
            flash_val = self.dat_manager.np_flash_v[flash_idx]
            flash_name = 'Flash (%dE2 PEs)' % int(np.sum(flash_val)/100.)

        if hypo_val is not None and flash_val is not None:
            #numerator = (flash_val - hypo_val)
            #denominator = (flash_val + hypo_val)
            #diff_val = np.zeros(shape=numerator.shape,dtype=np.float32)
            #where = np.where(denominator>0)
            #diff_val[where] = numerator[where] / denominator[where] * 2.0
            diff_val = flash_val - hypo_val
            diff_name = 'Diff (Flash-Hypothesis)'

        # Now decide the color range and normalization
        flash_range,hypo_range = (1.e20,-1.e20), (1.e20,-1.e20)
        flash_range = (1.e20,-1.e20) if flash_val is None else (flash_val.min(),flash_val.max())
        hypo_range  = (1.e20,-1.e20) if hypo_val  is None else (hypo_val.min(), hypo_val.max() )
        if normalize_pe == 0:
            if pe_range[0].lower() == 'min':
                pe_range[0] = min(flash_range[0],hypo_range[0])
            else:
                pe_range[0] = float(pe_range[0])
            if pe_range[1].lower() == 'max':
                pe_range[1] = max(flash_range[1],hypo_range[1])
            else:
                pe_range[1] = float(pe_range[1])
        if normalize_pe == 1 and not flash_val is None:
            pe_range = [0.,1.]
            norm = flash_val.sum()
            if not diff_val  is None: diff_val   = diff_val / flash_val
            if not flash_val is None: flash_val /= norm
            if not hypo_val  is None: hypo_val  /= norm
        if normalize_pe == 2 and not hypo_val is None:
            pe_range = [0.,1.]
            norm = hypo_val.sum()
            if not diff_val  is None: diff_val   = diff_val / hypo_val
            if not flash_val is None: flash_val /= norm
            if not hypo_val  is None: hypo_val  /= norm

        # make a trace for hypothesis
        pmt_pos = self.vis_icarus.data()['pmts']
        if not hypo_val is None:
            hypo_trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                      name=hypo_name,
                                      #template='plotly' if not self.dark_mode else 'plotly_dark',
                                      marker = dict(size=6,
                                                    color=hypo_val,
                                                    colorscale="OrRd",
                                                    cmin=pe_range[0],cmax=pe_range[1],
                                                    opacity=0.5),
                                      hoverinfo = ['x','y','z','text'],
                                      hovertext = ['%.2f' % val for val in hypo_val]
                                     )
            data.append(hypo_trace)

        if not flash_val is None:
            flash_trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                       name=flash_name,
                                       #template='plotly' if not self.dark_mode else 'plotly_dark',
                                       marker = dict(size=6,
                                                     color=flash_val,
                                                     colorscale = 'Greens',
                                                     cmin=pe_range[0], cmax=pe_range[1],
                                                     opacity=0.5
                                                    ),
                                      hoverinfo = ['x','y','z','text'],
                                      hovertext = ['%.2f' % val for val in flash_val]
                                      )
            data.append(flash_trace)

        if not diff_val is None:
            diff_trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                      name=diff_name,
                                      #template='plotly' if not self.dark_mode else 'plotly_dark',
                                      marker = dict(size=6,
                                                    color=diff_val,
                                                    colorscale=None,
                                                    cmin=pe_range[1]*(-1), cmax=pe_range[1],
                                                    opacity=0.5,),
                                      hoverinfo = ['x','y','z','text'],
                                      hovertext = ['%.2f' % val for val in diff_val]
                                     )
            data.append(diff_trace)

        sys.stdout.flush()
        if len(data)<1:
            return self.empty_view
        return go.Figure(self.detector_trace_no_pmt + data, layout=self.layout)

    def run_qll(self, data_index, qcluster_idx, flash_idx, is_entry, xpos):
        """
        Generate 3D display of a flash hypothesis for a given QCluster and x-offset
        """
        if not self.dat_manager.update(self.ana_manager, entry=data_index,
                                       is_entry=is_entry, make_hypothesis=False):
            return None
        data = []
        if qcluster_idx is None or flash_idx is None:
            return None
        # Find QCluster
        qcluster = self.dat_manager.cpp.qcluster_v[qcluster_idx]
        raw_qcluster = self.dat_manager.cpp.raw_qcluster_v[qcluster_idx]
        # find x span allowed
        xspan = qcluster.max_x() - qcluster.min_x()
        # if xoffset is larger than allowed range, truncate
        active_bb  = flashmatch.DetectorSpecs.GetME().ActiveVolume()
        max_xpos = min(active_bb.Min()[0] + xspan, active_bb.Max()[0])
        xpos     = max(min(xpos,max_xpos),active_bb.Min()[0])
        # shift to xoffset
        offset_qcluster = flashmatch.QCluster_t(qcluster) - qcluster.min_x() + xpos
        # Make hypothesis
        hypothesis = self.ana_manager.flash_hypothesis(offset_qcluster)
        # Find flash
        flash = self.dat_manager.cpp.flash_v[flash_idx]

        # Run QLL
        score_xoffset = self.qll_score(hypothesis,flash)
        match = self.run_flashmatch(qcluster,flash)
        sys.stdout.flush()
        return raw_qcluster.min_x(), score_xoffset, match
