import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from .view_api import AppManager

def view_data(cfg,geo,data_particle,data_opflash,dark_mode=True,port=5000):
    """
    Visualized OpT0Finder input data as well as FlashHypothesis for varying x-position
    INPUT:
      - cfg ... OpT0Finder configuration file, passed onto AnalysisManager
      - geo ... vis_icarus (geometry) data file in either PSet or yaml format
      - data_particle ... particle data file stored by ICARUSParticleAna_module
      - data_opflash ... OpFlash data file stored by ICARUSOpFlashAna_module
      - dark_mode ... use dark theme dash app
      - port ... set the custom port id for dash app
    """
    _manager = AppManager(cfg=cfg, geo=geo, data_particle=data_particle, data_opflash=data_opflash, dark_mode=dark_mode)
    entry = 0
    xmin = _manager.vis_icarus.data()['tpc0'][0][0]
    xmax = _manager.vis_icarus.data()['tpc1'][1][0]

    #
    # Create app
    #
    #app = dash.Dash('Flash Match Event Viewer')#,external_stylesheets=[dbc.themes.DARKLY])
    app = dash.Dash('Flash Match Event Viewer',
                    external_stylesheets= ['https://cdnjs.cloudflare.com/ajax/libs/animate.css/3.5.2/animate.min.css)']
                   )
    header_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=26)
    label_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=20)

    app.layout = html.Div([
        html.H2('Flash Matching Data Viewer', style=header_text_style),

        html.Div([html.Label('View data entry: (0-%d)' % (len(_manager.ana_manager.entries())-1),
                             style=label_text_style)
                 ],style={'padding': '5px'}),

        html.Div([dcc.RadioItems(id='entry_or_event',
                                 options=[dict(label='Entry Number',value='entry'),
                                          dict(label='Event ID',value='event')],
                                 value='entry'),
                  dcc.Input(id='data_index', value=str(entry), type='int')],
                 style={'width': '100%', 'display': 'inline-block', 'padding': '5px'}),

        dcc.Loading(id="loading", children=[html.Div(id="wait_out")], type="default"),

        html.Div([html.Label('Show also raw QClusters?',
                             style=label_text_style),
                  dcc.RadioItems(id='mode_qcluster',
                                 options=[dict(label='No',value='no_raw_qcluster'),
                                          dict(label='With Raw QCluster',value='with_raw_qcluster'),
                                          dict(label='Only Raw QCluster',value='only_raw_qcluster')],
                                 value='no_raw_qcluster'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div([html.Label('Raw QCluster to be shown should be:',
                             style=label_text_style),
                  dcc.RadioItems(id='use_all_pts',
                                 options=[dict(label='Before X-shift',value='raw_qcluster'),
                                          dict(label='... and before active BB cut',value='all_pts')],
                                 value='raw_qcluster'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div([html.Label('Select QCluster_t to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_qcluster',
                               options=_manager.dropdown_qcluster(entry,is_entry=True),
                               multi=True),
                 ],style={'padding': '5px'}),

        html.Div([html.Label('View PMTs with OpFlash or Hypothesis?',
                             style=label_text_style),
                  dcc.RadioItems(id='mode_flash',
                                 options=[dict(label='Flash',value='flash'),
                                          dict(label='Hypothesis',value='hypothesis')],
                                 value='flash'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div([html.Label('Select Flash_t to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_flash',
                               options=_manager.dropdown_flash(entry,is_entry=True,mode_flash=True),
                               multi=True),
                 ],style={'padding': '5px'}),
        #html.Div([html.Label('Event display (all selected flash combined)')])
        html.Div([dcc.Graph(id='visdata',figure=_manager.empty_view)],style={'padding': '5px'}),
        #
        # Hypothesis event display
        #
        html.H2('Hypothesis playground: move QCluster along x!',style=header_text_style),

        html.Div([html.Label('X Offset', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Slider(id='xoffset',min=xmin,max=xmax,step=1.,value=xmin),
                 ],style={'padding': '5px'}),

        html.Div([html.Label('Select QCluster to display', style=label_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='target_qcluster',
                               options=_manager.dropdown_qcluster(entry,is_entry=True),
                               multi=True),
                 ],style={'padding': '5px'}),
        #dcc.Loading(id="loading-playdata",
        #            children=[html.Div([dcc.Graph(id='playdata',figure=_manager.empty_view)],
        #                               style={'padding': '5px'})],
        #            type="circle"),
        html.Div([dcc.Graph(id='playdata',figure=_manager.empty_view)],style={'padding': '5px'}),
    ],style={'width': '80%', 'display': 'inline-block', 'vertical-align': 'middle'})

    #
    # call backs
    #
    @app.callback(dash.dependencies.Output("wait_out", "children"),
                  [dash.dependencies.Input("mode_flash", "value")])
    def load_plib(mode):
        if mode == 'hypothesis':
            from flashmatch import phot
            print('hypothesis mode chosen, loading photon library...')
            phot.PhotonVisibilityService.GetME().LoadLibrary()
        return ""

    @app.callback(dash.dependencies.Output("select_qcluster","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def update_dropdown_select_qcluster(entry_or_event,data_index):
        """
        Update drop-down "select_qcluster" options for QCluster when event/entry is changed
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None or int(data_index) == _manager.current_data_index(is_entry=is_entry):
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_qcluster(data_index=int(data_index),is_entry=is_entry)
        if options is None:
            print('Cannot find data for',entry_or_event,data_index)
            raise dash.exceptions.PreventUpdate
        else: return options


    @app.callback(dash.dependencies.Output("target_qcluster","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def update_dropdown_target_qcluster(entry_or_event,data_index):
        """
        Update drop-down "target_qcluster" options for QCluster when event/entry is changed
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None or int(data_index) == _manager.current_data_index(is_entry=is_entry):
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_qcluster(data_index=int(data_index),is_entry=is_entry)
        if options is None: raise dash.exceptions.PreventUpdate
        else: return options


    @app.callback(dash.dependencies.Output("select_flash","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value"),
                   dash.dependencies.Input("mode_flash","value")])
    def update_dropdown_flash(entry_or_event,data_index,mode_flash):
        """
        Update drop-down "select_flash" options for Flash when event/entry is changed
        """
        is_entry = entry_or_event == 'entry'
        mode_flash = mode_flash == 'hypothesis'
        if data_index is None or (int(data_index) == _manager.current_data_index(is_entry=is_entry)
                                  and not mode_flash):
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_flash(data_index=int(data_index),
                                          is_entry=is_entry,
                                          mode_flash=mode_flash)
        if options is None: raise dash.exceptions.PreventUpdate
        else: return options


    @app.callback(dash.dependencies.Output("visdata","figure"),
                  [dash.dependencies.Input("select_qcluster","value"),
                   dash.dependencies.Input("select_flash","value"),
                   dash.dependencies.Input("mode_flash","value"),
                   dash.dependencies.Input("mode_qcluster","value"),
                   dash.dependencies.Input("use_all_pts","value"),
                   dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")]
                 )
    def update_static(select_qcluster, select_flash, mode_flash, mode_qcluster, use_all_pts, entry_or_event, data_index):
        """
        Update "visdata" 3D display for change in event/entry and/or selection of flash/qcluster
        """
        is_entry = entry_or_event == 'entry'
        mode_flash = mode_flash == 'flash'
        if mode_qcluster == 'no_raw_qcluster': mode_qcluster = 0
        elif mode_qcluster == 'with_raw_qcluster': mode_qcluster=1
        elif mode_qcluster == 'only_raw_qcluster': mode_qcluster=2
        else: raise ValueError
        use_all_pts = True if use_all_pts == 'all_pts' else False
        if data_index is None: raise dash.exceptions.PreventUpdate
        fig = _manager.event_display(int(data_index), select_qcluster, select_flash, is_entry, mode_flash, mode_qcluster, use_all_pts)
        if fig is None: raise dash.exceptions.PreventUpdate
        else: return fig


    @app.callback(dash.dependencies.Output("playdata","figure"),
                  [dash.dependencies.Input("target_qcluster","value"),
                   dash.dependencies.Input("xoffset","value"),
                   dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")]
                 )
    def update_dynamic(target_qcluster, xoffset, entry_or_event, data_index):
        """
        Update "playdata" 3D display for change in event/entry and/or selection of flash/qcluster
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None or target_qcluster is None or len(target_qcluster)<1:
            raise dash.exceptions.PreventUpdate
        fig = _manager.hypothesis_display(int(data_index), target_qcluster, is_entry, float(xoffset))
        if fig is None: raise dash.exceptions.PreventUpdate
        else: return fig


    app.server.run(port=port)
