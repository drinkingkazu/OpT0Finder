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
            entry = manager.entry_id(entry)
        if entry < 0 or entry >= len(manager.entries()): return
        if not self.entry == entry:
            self.cpp = manager.make_flashmatch_input(entry)
            self.np_qcluster_v     = [flashmatch.as_ndarray(qcluster) for qcluster in self.cpp.qcluster_v    ]
            self.np_raw_qcluster_v = [flashmatch.as_ndarray(qcluster) for qcluster in self.cpp.raw_qcluster_v]
            self.np_flash_v        = [flashmatch.as_ndarray(flash)    for flash    in self.cpp.flash_v       ]
            self.entry = entry
            self.event = manager.event_id(entry)
            self.hypothesis_made = False
        if make_hypothesis and not self.hypothesis_made:
            self.cpp.hypo_v = [manager.flash_hypothesis(qcluster) for qcluster in self.cpp.qcluster_v]
            self.np_hypo_v  = [flashmatch.as_ndarray(flash)       for flash    in self.cpp.hypo_v    ]
            for idx in range(len(self.np_hypo_v)):
                self.cpp.hypo_v[idx].idx = self.cpp.qcluster_v[idx].idx

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
        self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry)
        idx_v = [qcluster.idx for qcluster in self.dat_manager.cpp.qcluster_v]
        dropdown_qcluster = [dict(label='Track %02d (%d pts)' % (idx_v[idx],len(qcluster)),
                                  value=idx)
                             for idx,qcluster in enumerate(self.dat_manager.np_qcluster_v)]
        dropdown_qcluster += [dict(label='All tracks',value=len(self.dat_manager.np_qcluster_v))]
        return dropdown_qcluster

    def dropdown_flash(self,data_index,is_entry,mode_flash):
        """
        Generates a dropdown menu for a dash app to select a Flash
        INPUTS:
          - data_index ... an integer to specify which entry/event to read
          - is_entry ... True=entry, False=event id
        OUTPUT:
          A list of dicts to be consumed by dash app dropdown options
        """
        self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry, make_hypothesis=(not mode_flash))
        target_v = self.dat_manager.np_flash_v if mode_flash else self.dat_manager.np_hypo_v
        idx_v  = [flash.idx for flash in self.dat_manager.cpp.flash_v]
        time_v = np.zeros(shape=[len(target_v)])
        if mode_flash: time_v = [flash.time for flash in self.dat_manager.cpp.flash_v]
        dropdown_flash  = [dict(label='Flash %02d (t=%f us  PE=%d)' % (idx_v[idx],
                                                                       time_v[idx],
                                                                       int(flash.sum()/100.)
                                                                      ),
                                value=idx)
                            for idx,flash in enumerate(target_v)]
        dropdown_flash += [dict(label='All flashes',value=len(target_v))]
        return dropdown_flash

    def event_display(self, data_index, qcluster_idx_v, flash_idx_v, is_entry, mode_flash, mode_qcluster):
        """
        Generate 3D display for change in event/entry and/or selection of flash/qcluster
        """
        self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry, make_hypothesis=(not mode_flash))
        data = []
        print('qcluster',mode_qcluster)
        # QCluster
        # if mode_qcluster == 0, then only show qcluster with color-per-cluster
        # if mode_qcluster == 1, then show both qcluster and raw version with 2 colors (red + green)
        # if mode_qcluster == 2, then only show raw qcluster with color-per-cluster

        if qcluster_idx_v is not None and len(qcluster_idx_v):
            if len(self.dat_manager.np_qcluster_v) in qcluster_idx_v:
                qcluster_idx_v = range(len(_blob.np_qcluster_v))
            idx_v = [self.dat_manager.cpp.qcluster_v[idx].idx for idx in qcluster_idx_v]
            for idx in qcluster_idx_v:
                if mode_qcluster in [0,1]:
                    xyz = self.dat_manager.np_qcluster_v[idx]
                    name = 'Track %02d (%d pts)' % (self.dat_manager.cpp.qcluster_v[idx].idx,len(xyz))
                    trace = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
                                         name=name,
                                         #template='plotly' if not self.dark_mode else 'plotly_dark',
                                         marker = dict(size=2, opacity=0.5, color=None if mode_qcluster==0 else 'Red')
                                        )
                    data.append(trace)
                if mode_qcluster in [1,2]:
                    xyz = self.dat_manager.np_raw_qcluster_v[idx]
                    name = 'Track (Raw) %02d (%d pts)' % (self.dat_manager.cpp.qcluster_v[idx].idx,len(xyz))
                    trace = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
                                         name=name,
                                         #template='plotly' if not self.dark_mode else 'plotly_dark',
                                         marker = dict(size=2, opacity=0.5, color=None if mode_qcluster==2 else 'Blue')
                                        )
                    data.append(trace)
        target_v = self.dat_manager.np_flash_v if mode_flash else self.dat_manager.np_hypo_v
        cpp_target_v = self.dat_manager.cpp.flash_v if mode_flash else self.dat_manager.cpp.hypo_v
        if flash_idx_v and len(flash_idx_v):
            pmt_pos = self.vis_icarus.data()['pmts']
            pmt_val = None
            if len(target_v) in flash_idx_v: pmt_val = np.sum(target_v,axis=0)
            else: pmt_val = np.sum(np.column_stack([target_v[idx] for idx in flash_idx_v]),axis=1)
            name = 'Flash (%f PEs)' % np.sum(pmt_val)
            if(len(flash_idx_v)==1):
                name = 'Flash (%f PEs @ %f us)' % (np.sum(pmt_val), cpp_target_v[flash_idx_v[0]].time)
            trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                 name=name,
                                 #template='plotly' if not self.dark_mode else 'plotly_dark',
                                 marker = dict(size=6, color=pmt_val, opacity=0.5)
                                )
            data.append(trace)
        elif len(data):
            # no pmt data => show default pmt view
            return go.Figure(self.detector_trace + data, layout=self.layout)

        if len(data)<1:
            return self.empty_view
        else:
            return go.Figure(self.detector_trace_no_pmt + data, layout=self.layout)

    def hypothesis_display(self, data_index, qcluster_idx_v, is_entry, xoffset):
        """
        Generate 3D display of a flash hypothesis for a given QCluster and x-offset
        """
        self.dat_manager.update(self.ana_manager, entry=data_index, is_entry=is_entry, make_hypothesis=False)
        data = []
        if qcluster_idx_v is None or len(qcluster_idx_v)<1: return self.empty_view
        if len(self.dat_manager.np_qcluster_v) in qcluster_idx_v:
            qcluster_idx_v = range(len(self.dat_manager.np_qcluster_v))
        qcluster = flashmatch.QCluster_t()
        for idx in qcluster_idx_v:
            qcluster += self.dat_manager.cpp.qcluster_v[idx]
        # find x span allowed
        xmin,xmax=qcluster.min_x(),qcluster.max_x()
        # if xoffset is larger than xmax-xmin, truncate
        xoffset = min(xoffset,xmax-xmin)
        # shift
        qcluster += (xoffset - xmin)
        # make hypothesis
        hypothesis = self.ana_manager.flash_hypothesis(qcluster)
        # make a trace for tpc
        xyzv = flashmatch.as_ndarray(qcluster)
        tpc_trace = go.Scatter3d(x=xyzv[:,0],y=xyzv[:,1],z=xyzv[:,2],mode='markers',
                                 name='Track (offset %f)' % xoffset,
                                 #template='plotly' if not self.dark_mode else 'plotly_dark',
                                 marker = dict(size=2,opacity=0.5)
                                )
        # make a trace for pmt
        pmt_val = flashmatch.as_ndarray(hypothesis)
        pmt_pos = self.vis_icarus.data()['pmts']
        pmt_trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                 name='Hypothesis (%f PEs)' % np.sum(pmt_val),
                                 #template='plotly' if not self.dark_mode else 'plotly_dark',
                                 marker = dict(size=6,
                                               color=pmt_val,
                                               colorscale=None,
                                               opacity=0.5)
                                )
        return go.Figure(self.detector_trace_no_pmt + [pmt_trace] + [tpc_trace], layout=self.layout)
