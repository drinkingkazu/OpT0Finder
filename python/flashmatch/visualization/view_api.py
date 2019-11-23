import numpy as np
import plotly.graph_objs as go
from flashmatch import flashmatch, AnalysisManager
from .vis_icarus import vis_icarus, icarus_layout3d

class DataManager:
    """
    DataManager class is meant to hold a blob data for OpT0Finder. In particular, it holds:
     - "current" entry and event id
     - "cpp" QCluster,Flash_t, etc. in flashmatch.FlashMatchInput format (this is self.cpp attribute)
     - numpy array representation of what's stored in "cpp"
     - Read & update local data attributes is done via self.update function
    """
    def __init__(self):
        self.entry = -1
        self.event = -1
    def update(self,manager,entry,is_entry=True):
        """
        Given flashmatch.AnalysisManager, which can interface input data stream and OpT0Finder tools, this
        function "read data and update local data attributes".
        """
        if not is_entry:
            entry = manager.entry_id(entry)
        if entry < 0 or entry >= len(manager.entries()): return
        if self.entry == entry: return
        self.cpp = manager.make_flashmatch_input(entry)
        self.cpp.hypo_v    = [manager.flash_hypothesis(qcluster) for qcluster in self.cpp.qcluster_v]
        self.np_qcluster_v = [flashmatch.as_ndarray(qcluster)    for qcluster in self.cpp.qcluster_v]
        self.np_flash_v    = [flashmatch.as_ndarray(flash)       for flash    in self.cpp.flash_v   ]
        self.np_hypo_v     = [flashmatch.as_ndarray(flash)       for flash    in self.cpp.hypo_v    ]
        self.entry = entry
        self.event = manager.event_id(entry)

class VisManager:
    """
    VisManager is a mere object instance holder that are useful for data visualization 
    """
    def __init__(self,geo):
        self.vis_icarus = vis_icarus(geo)
        self.detector_trace = self.vis_icarus.get_trace_detector()#draw_pmts=False)
        self.layout = icarus_layout3d(self.vis_icarus.data(),set_camera=False,dark=True)
        self.empty_view = go.Figure(self.detector_trace,layout=self.layout)
