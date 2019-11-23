import numpy as np
import plotly.graph_objs as go
from flashmatch import flashmatch, AnalysisManager
from .vis_icarus import vis_icarus, icarus_layout3d
import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from .view_api import VisManager, DataManager

def view_data(cfg,geo,data_particle,data_opflash):
    _manager = AnalysisManager(cfg=cfg, particleana=data_particle, opflashana=data_opflash)
    _vis = VisManager(geo)
    assert(len(_manager.entries()))
    _entry = 0
    _xmin = _vis.vis_icarus.data()['tpc0'][0][0]
    _xmax = _vis.vis_icarus.data()['tpc1'][1][0]
    _blob = DataManager()
    _blob.update(_manager,_entry)

    #
    # Create app
    #
    app = dash.Dash('Flash Match Event Viewer')#,external_stylesheets=[dbc.themes.DARKLY])
    header_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=26)
    label_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=20)

    app.layout = html.Div([
        html.H2('Flash Matching Data Viewer', style=header_text_style),

        html.Div([html.Label('View data entry: (0-%d)' % (len(_manager.entries())-1),
                             style=label_text_style)
                 ],style={'padding': '5px'}),

        html.Div([dcc.RadioItems(id='entry_or_event',
                                 options=[dict(label='Entry Number',value='entry'),
                                          dict(label='Event ID',value='event')],
                                 value='entry'),
                  dcc.Input(id='data_index', value=str(0), type='int')],
                 style={'width': '100%', 'display': 'inline-block', 'padding': '5px'}),

        html.Div([html.Label('View PMTs with OpFlash or Hypothesis?',
                             style=label_text_style),
                  dcc.RadioItems(id='flash_or_hypo',
                                 options=[dict(label='Flash',value='flash'),
                                          dict(label='Hypothesis',value='hypothesis')],
                                 value='Flash'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),

        html.Div([html.Label('Select QCluster_t to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_qcluster',multi=True),
                 ],style={'padding': '5px'}),

        html.Div([html.Label('Select Flash_t to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_flash',multi=True),
                 ],style={'padding': '5px'}),

        html.Div([dcc.Graph(id='visdata',figure=_vis.empty_view)],style={'padding': '5px'}),
        html.H2('Hypothesis playground: move QCluster along x!',style=header_text_style),
        html.Div([html.Label('X Offset', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Slider(id='xoffset',min=_xmin,max=_xmax,step=1.,value=_xmin),
                 ],style={'padding': '5px'}),
        html.Div([html.Label('Select QCluster to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='target_qcluster',multi=True),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Graph(id='playdata',figure=_vis.empty_view)],style={'padding': '5px'}),
    ],style={'width': '80%', 'display': 'inline-block', 'vertical-align': 'middle'})

    # call backs

    @app.callback(dash.dependencies.Output("select_qcluster","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def update_dropdown_select_qcluster(entry_or_event,data_index):
        data_index=int(data_index)
        _blob.update(_manager, data_index, is_entry=(entry_or_event=="entry"))
        dropdown_qcluster = [dict(label='Track %02d (%d pts)' % (_blob.cpp.qcluster_v[idx].idx,
                                                                 len(qcluster)
                                                                ),
                                  value=idx)
                             for idx,qcluster in enumerate(_blob.np_qcluster_v)]
        dropdown_qcluster += [dict(label='All tracks',value=len(_blob.np_qcluster_v))]
        return dropdown_qcluster

    @app.callback(dash.dependencies.Output("target_qcluster","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def update_dropdown_target_qcluster(entry_or_event,data_index):
        data_index=int(data_index)
        _blob.update(_manager, data_index, is_entry=(entry_or_event=="entry"))
        dropdown_qcluster = [dict(label='Track %02d (%d pts)' % (_blob.cpp.qcluster_v[idx].idx,
                                                                 len(qcluster)
                                                                ),
                                  value=idx)
                             for idx,qcluster in enumerate(_blob.np_qcluster_v)]
        dropdown_qcluster += [dict(label='All tracks',value=len(_blob.np_qcluster_v))]
        return dropdown_qcluster

    @app.callback(dash.dependencies.Output("select_flash","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value"),
                   dash.dependencies.Input("flash_or_hypo","value")])
    def update_dropdown_qcluster(entry_or_event,data_index,flash_or_hypo):
        data_index=int(data_index)
        _blob.update(_manager, data_index, is_entry=(entry_or_event=="entry"))
        target_v = _blob.np_hypo_v if flash_or_hypo == 'hypothesis' else _blob.np_flash_v
        dropdown_flash  = [dict(label='Flash %02d (%dE2 p.e.)' % (_blob.cpp.flash_v[idx].idx,
                                                                  int(flash.sum()/100.)
                                                                 ),
                                value=idx)
                            for idx,flash in enumerate(target_v)]
        dropdown_flash += [dict(label='All flashes',value=len(target_v))]
        return dropdown_flash

    @app.callback(dash.dependencies.Output("visdata","figure"),
                 [dash.dependencies.Input("entry_or_event","value"),
                  dash.dependencies.Input("select_qcluster","value"),
                  dash.dependencies.Input("select_flash","value"),
                  dash.dependencies.Input("flash_or_hypo","value")]
                 )
    def update_static(entry_or_event, select_qcluster, select_flash, flash_or_hypo):
        _blob.update(_manager, _entry, is_entry=(entry_or_event=="entry"))
        data = []
        if select_qcluster is not None and len(select_qcluster):
            if len(_blob.np_qcluster_v) in select_qcluster:
                select_qcluster = range(len(_blob.np_qcluster_v))
            for idx in select_qcluster:
                xyz = _blob.np_qcluster_v[idx]
                trace = go.Scatter3d(x=xyz[:,0],y=xyz[:,1],z=xyz[:,2],mode='markers',
                                     name='Track %02d (%d pts)' % (_blob.cpp.qcluster_v[idx].idx,
                                                                   len(xyz)
                                                                  ),
                                     marker = dict(size=2, opacity=0.5)
                                    )
                data.append(trace)
        target_v = _blob.np_hypo_v if flash_or_hypo == 'hypothesis' else _blob.np_flash_v
        if select_flash and len(select_flash):
            pmt_pos = _vis.vis_icarus.data()['pmts']
            pmt_val = None
            if len(target_v) in select_flash: pmt_val = np.sum(target_v,axis=0)
            else: pmt_val = np.sum(np.column_stack([target_v[idx] for idx in select_flash]),axis=1)
            trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                 name='Flash (%f PEs)' % np.sum(pmt_val),
                                 marker = dict(size=6, color=pmt_val, opacity=0.5)
                                )
            data.append(trace)

        if len(data)<1:
            return _vis.empty_view
        else:
            return go.Figure(_vis.detector_trace + data, layout=_vis.layout)


    @app.callback(dash.dependencies.Output("playdata","figure"),
                 [dash.dependencies.Input("entry_or_event","value"),
                  dash.dependencies.Input("target_qcluster","value"),
                  dash.dependencies.Input("xoffset","value")]
                 )
    def update_dynamic(entry_or_event,target_qcluster,xoffset):
        _blob.update(_manager, _entry, is_entry=(entry_or_event=="entry"))
        data = []
        if not target_qcluster or len(target_qcluster)<1: return _vis.empty_view
        if len(_blob.np_qcluster_v) in target_qcluster:
            target_qcluster = range(len(_blob.np_qcluster_v))
        qcluster = flashmatch.QCluster_t()
        for idx in target_qcluster:
            qcluster += _blob.cpp.qcluster_v[idx]
        # find x span allowed
        xmin,xmax=qcluster.min_x(),qcluster.max_x()
        # if xoffset is larger than xmax-xmin, truncate
        xoffset = min(xoffset,xmax-xmin)
        # shift
        qcluster += (xoffset - xmin)
        # make hypothesis
        hypothesis = _manager.flash_hypothesis(qcluster)
        # make a trace for tpc
        xyzv = flashmatch.as_ndarray(qcluster)
        tpc_trace = go.Scatter3d(x=xyzv[:,0],y=xyzv[:,1],z=xyzv[:,2],mode='markers',
                                 name='Track (offset %f)' % xoffset,
                                 marker = dict(size=2,
                                               #color='red',
                                               opacity=0.5)
                                )
        # make a trace for pmt
        pmt_val = flashmatch.as_ndarray(hypothesis)
        pmt_pos = _vis.vis_icarus.data()['pmts']
        pmt_trace = go.Scatter3d(x=pmt_pos[:,0],y=pmt_pos[:,1],z=pmt_pos[:,2],mode='markers',
                                 name='Hypothesis (%f PEs)' % np.sum(pmt_val),
                                 marker = dict(size=6,
                                               color=pmt_val,
                                               colorscale=None,
                                               opacity=0.5)
                                )
        return go.Figure(_vis.detector_trace + [pmt_trace] + [tpc_trace], layout=_vis.layout)

    import sys
    sys.stdout.write('Wait for PhotonLibrary loading...\n')
    sys.stdout.write('If you are developing view_data, consider doing in Jupyter notebook w/o re-loading library each time...\n')
    sys.stdout.flush()
    app.server.run()
