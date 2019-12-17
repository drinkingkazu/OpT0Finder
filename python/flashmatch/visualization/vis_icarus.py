import numpy as np
import plotly.graph_objs as go

def load_geometry_from_pset(cfg):
    import ast
    from flashmatch import flashmatch
    pset = flashmatch.CreatePSetFromFile(cfg).get('flashmatch::PSet')('DetectorSpecs')
    det = {}
    det['MinPosition'] = ast.literal_eval(pset.get('PhotonLibraryVolumeMin'))
    det['MaxPosition'] = ast.literal_eval(pset.get('PhotonLibraryVolumeMax'))
    det['MinPosition_TPC0'] = ast.literal_eval(pset.get('MinPosition_TPC0'))
    det['MaxPosition_TPC0'] = ast.literal_eval(pset.get('MaxPosition_TPC0'))
    det['MinPosition_TPC1'] = ast.literal_eval(pset.get('MinPosition_TPC1'))
    det['MaxPosition_TPC1'] = ast.literal_eval(pset.get('MaxPosition_TPC1'))
    pmt = 0
    while 1:
        if not pset.contains_value('PMT%d' % pmt):
            break
        det['PMT%d' % pmt] = ast.literal_eval(pset.get('PMT%d' % pmt))
        pmt+=1
    return det

def load_geometry_from_yaml(cfg):
    import yaml
    return yaml.load(open(cfg).read(),Loader=yaml.Loader)['DetectorSpecs']

def load_geometry(data=None):
    """
    Produces geometry information dictionary with a given input data file
    INPUTS
        - data is a path to a geometry data file that can be either yaml or larcv/opt0finder PSet format.
          Necessary elements include MinPosition_TPCX, MaxPosition_TPCX for X=0 and 1 (for cryostat 0).
          Also PMTX for X 0 to 179. Finally the MinPosition_Cryo0, MaxPosition_Cryo0 denoting
          cryostat 0 min and max points.
    OUTPUT
        - geo is a dictionary with keys 'anode0', 'anode1', 'cathode0', 'cryo0' with an associated
          numpy arrays of shape (2,3) for min and max xyz point coordinates. Finally 'pmt' is associated
          to a numpy array of shape (180,3) for 180 PMT xyz coordinates.
    """
    import os
    if data is None: data = os.path.join(os.environ['VIS_ICARUS'],'dat','icarus.yaml')
    det=None
    if data.endswith('.yaml') or data.endswith('.yml'):
        det = load_geometry_from_yaml(data)
    elif data.endswith('.cfg') or data.endswith('.fcl'):
        det = load_geometry_from_pset(data)
    else:
        print('Config file not in supported format:',data)
        raise TypeError
    geo={}
    # create a table for PMT
    geo['pmts']=np.array([det['PMT%d' % i] for i in range(180)]).astype(np.float32)
    # create a table for anode for tpc0
    a0_min=det['MinPosition_TPC0']
    a0_max=det['MaxPosition_TPC0']
    geo['tpc0']=np.array([a0_min,a0_max]).astype(np.float32)
    # create a table for anode for anode0
    geo['anode0']=np.array([[a0_min[0],a0_min[1],a0_min[2]],[a0_min[0],a0_max[1],a0_max[2]]]).astype(np.float32)
    # create a table for anode for tpc1
    a1_min=det['MinPosition_TPC1']
    a1_max=det['MaxPosition_TPC1']
    geo['tpc1']=np.array([a1_min,a1_max]).astype(np.float32)
    # create a table for anode for anode1
    geo['anode1']=np.array([[a1_max[0],a1_min[1],a1_min[2]],[a1_max[0],a1_max[1],a1_max[2]]]).astype(np.float32)
    # assert a0 and a1 yz bounds
    assert a0_min[1] == a1_min[1] and a0_min[2] == a1_min[2]
    assert a0_max[1] == a1_max[1] and a0_max[2] == a1_max[2]
    # create a table for cathode for cryo0
    cathode_x = a0_max[0]+(a1_min[0]-a0_max[0])/2.
    geo['cathode0']=np.array([[cathode_x,a0_min[1],a0_min[2]],[cathode_x,a0_max[1],a0_max[2]]]).astype(np.float32)
    # create a table for cryo0
    c0_min=det['MinPosition']
    c0_max=det['MaxPosition']
    geo['cryo0']=np.array([c0_min,c0_max]).astype(np.float32)
    return geo


def icarus_layout3d(ranges=None,set_camera=True,dark=False):
    """
    Produces go.Layout object for a ICARUS event display.
    INPUTS
        - ranges can be used to specify the plot region in (x,y,z) directions.
          The default (None) will determine the range to include all points.
          Alternatively can be an array of shape (3,2) specifying (x,y,z) axis (min,max) range for a display,
          or simply a list of points with shape (N,3+) where [:,0],[:,1],[:,2] correspond to (x,y,z) values and
          the plotting region is decided by measuring the min,max range in each coordinates. This last option
          is useful if one wants to define the region based on a set of points that is not same as what's plotted.
    OUTPUTS
        - The return is go.Layout object that can be given to go.Figure for visualization (together with traces)

    """
    xrange,yrange,zrange=None,None,None
    if ranges is None:
        ranges=[None,None,None]
    elif type(ranges) == type(dict()):
        xrange,yrange,zrange=ranges['cryo0'][:,0],ranges['cryo0'][:,1],ranges['cryo0'][:,2]
    elif np.shape(ranges) == (3,2):
        xrange,yrange,zrange=ranges
    else:
        xrange = (np.min(ranges[:,0]),np.max(ranges[:,0]))
        yrange = (np.min(ranges[:,1]),np.max(ranges[:,1]))
        zrange = (np.min(ranges[:,2]),np.max(ranges[:,2]))

    scene = dict(xaxis = dict(nticks=10, range = xrange, showticklabels=True,
                              title='x',
                              #backgroundcolor="lightgray", gridcolor="rgb(255, 255, 255)",
                              #showbackground=True,
                             ),
                 yaxis = dict(nticks=10, range = yrange, showticklabels=True,
                              title='y',
                              #backgroundcolor="lightgray", gridcolor="rgb(255, 255, 255)",
                              #showbackground=True
                             ),
                 zaxis = dict(nticks=10, range = zrange, showticklabels=True,
                              title='z',
                              #backgroundcolor="lightgray", gridcolor="rgb(255, 255, 255)",
                              #showbackground=True,
                             ),
                 aspectmode='data',
                )

    if set_camera: scene['camera'] = dict(up=dict(x=0, y=1, z=0),
                                          center=dict(x=0, y=0, z=-1.),
                                          eye=dict(x=0.0, y=0.0, z=1.0),
                                          )
    layout = go.Layout(
        showlegend=True,
        legend=dict(x=1.01,y=0.95),
        width=1024,
        height=768,
        hovermode='closest',
        margin=dict(l=0,r=0,b=0,t=0),
        #plot_bgcolor  = '#111111',
        #paper_bgcolor = '#111111',
        #font = dict(color = '#7FDBFF'),
        template='plotly_dark' if dark else None,
        #uirevision = False,
        #uirevision = 'dataset',
        uirevision = 'same',
        scene = scene,
    )

    return layout

class vis_icarus(object):

    def __init__(self,data=None):
        self._geo = load_geometry(data)

    def data(self):
        return self._geo

    def get_trace_anode0(self):
        """
        Returns plotly 3D trace for ICARUS anode plane for TPC0
        INPUTS:
        OUTPUTS:
          a length 1 list of trace
        """
        data = self._geo['anode0']
        return go.Surface(x=data[:,0], y=data[:,1], z=[data[:,2],data[:,2]],name='anode0',
                          showscale=False,opacity=0.5,surfacecolor=[0,0])

    def get_trace_anode1(self):
        """
        Returns plotly 3D trace for ICARUS anode plane for TPC1
        INPUTS:
        OUTPUTS:
          a length 1 list of trace
        """
        data = self._geo['anode1']
        return go.Surface(x=data[:,0], y=data[:,1], z=[data[:,2],data[:,2]],name='anode1',
                          showscale=False,opacity=0.5,surfacecolor=[0,0])

    def get_trace_cathode0(self):
        """
        Returns plotly 3D trace for ICARUS cathode plane between TPC0 and TPC1
        INPUTS:
        OUTPUTS:
          a length 1 list of trace
        """
        data = self._geo['cathode0']
        return go.Surface(x=data[:,0], y=data[:,1], z=[data[:,2],data[:,2]],name='cathode0',
                          showscale=False,opacity=0.5,surfacecolor=[0,0])

    def get_trace_pmt(self,pmt_color='gray',pmt_opacity=0.5,pmt_colorscale=None,name=None):
        """
        Returns plotly 3D trace for ICARUS PMTs (for both TPC0 and TPC1)
        INPUTS:
        - pmt_color is either string or a 1D array of length 180 filled with color values
        - pmt_colorscale is a string, useful when providing color values
        - pmt_opacity is a floating point value specifying the opacity of PMTs
        OUTPUTS:
          a length 1 list of trace
        """
        data = self._geo['pmts']
        return go.Scatter3d(x=data[:,0],y=data[:,1],z=data[:,2],mode='markers',
                            name=name,
                            marker = dict(
                                size = 6,
                                color=pmt_color,
                                colorscale=pmt_colorscale,
                                opacity=pmt_opacity,
                            ))

    def get_trace_detector(self,draw_anode0=True,draw_anode1=True,draw_cathode0=True,draw_pmts=True,
                           pmt_value=None, pmt_opacity=0.5, pmt_colorscale=None):
        """
        Returns plotly 3D trace for ICARUS TPC0+TPC1+PMT in cryostat 0 structure.
        INPUTS:
        - draw_anode0 ... draws anode0 if True
        - draw_anode1 ... draws anode1 if True
        - draw_cathode0 ... draws cathode if True
        - draw_pmts ... draws PMTs if True
        - pmt_value ... an array of length 180 for PMT temperature to be drawn
        - pmt_colorscale ... a cmap for PMT heat color
        - pmt_opacity ... the opacity of PMTs
        OUTPUT:
          a list of traces
        """
        trace  = []
        if draw_anode0: trace += [self.get_trace_anode0()]
        if draw_anode1: trace += [self.get_trace_anode1()]
        if draw_cathode0: trace += [self.get_trace_cathode0()]
        if draw_pmts: trace += [self.get_trace_pmt(pmt_value,pmt_opacity,pmt_colorscale,name='PMTs')]
        return trace

def get_trace_icarus(pmt_value=None,geo_yaml=None,
                     draw_anode0=True,draw_anode1=True,draw_cathode0=True,draw_pmts=True,
                     return_geo=False):

    vis = vis_icarus(geo_yaml)
    trace = vis.get_trace_detector(draw_anode0=draw_anode0, draw_anode1=draw_anode1, draw_cathode0=draw_cathode0,
                                   draw_pmts=draw_pmts, pmt_value=pmt_value)

    return trace if not return_geo else (trace,vis.data())
