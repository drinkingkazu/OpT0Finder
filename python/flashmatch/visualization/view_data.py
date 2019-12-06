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
    span_text_style = dict(fontFamily='Georgia', fontWeight=600, fontSize=24, color='navy')
    label_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=20)
    para_text_style = dict(fontFamily='Georgia', fontWeight=300, fontSize=16)
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
        html.Div(id='error_entry'),
        #
        # TPC options
        #
        html.Div([html.Span('Select QCluster_t to display', style=span_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_qcluster',
                               options=_manager.dropdown_qcluster(entry,is_entry=True),
                               multi=True),
                 ],style={'padding': '5px'}),
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

        #
        # PMT options
        #
        html.Div([html.Span('Select Flash_t to display', style=span_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='select_flash',
                               options=_manager.dropdown_flash(entry,is_entry=True,mode_flash=True),
                               multi=True),
                 ],style={'padding': '5px'}),
        html.Div([html.Label('View PMTs with OpFlash or Hypothesis?',
                             style=label_text_style),
                  dcc.RadioItems(id='mode_flash',
                                 options=[dict(label='Flash',value='flash'),
                                          dict(label='Hypothesis',value='hypothesis')],
                                 value='flash'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div([html.Label('PMT color scale: ', style=label_text_style),
                  html.Br(),
                  dcc.Input(id='pmt_cmin', value='min', type='str'),# style={'padding':'5px'}),
                  html.Label(' to ', style=label_text_style),
                  dcc.Input(id='pmt_cmax', value='max', type='str'),# style={'padding':'5px'}),
                  html.Label(' (type "min" and "max" for automatic range setting)', style=para_text_style),
                 ], style={'width': '100%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div(id='error_pmt_range'),
        dcc.Loading(id="loading_static_display", children=[html.Div(id="wait_static_display")], type="default"),
        html.Div([dcc.Graph(id='visdata',figure=_manager.empty_view)],
                  style={'padding':'5px'}),
        #
        # Hypothesis event display
        #
        html.H2('Hypothesis playground: move QCluster along x!',style=header_text_style),

        html.Div([html.Span('Select QCluster to display', style=span_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='target_qcluster',
                               options=_manager.dropdown_qcluster(entry,is_entry=True),
                               multi=False),
                 ],style={'padding': '5px'}),
        html.Div([html.Span('Select Flash_t to match', style=span_text_style),
                 ],style={'padding': '5px'}),
        html.Div([dcc.Dropdown(id='target_flash',
                               options=_manager.dropdown_flash(entry,is_entry=True,mode_flash=True),
                               multi=False),
                 ],style={'padding': '5px'}),
        html.Div([html.Label('Normalize PE scale for score calculation',
                             style=label_text_style),
                  dcc.RadioItems(id='normalize_pe',
                                 options=[dict(label='No',value='0'),
                                          dict(label='Normalize to hypothesis',value='1'),
                                          dict(label='Normalize to flash',value='2')],
                                 value='0'),],
                 style={'width': '50%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div([html.Label('PMT color scale: ', style=label_text_style),
                  html.Br(),
                  dcc.Input(id='flash_cmin', value='min', type='str'),# style={'padding':'5px'}),
                  html.Label(' to ', style=label_text_style),
                  dcc.Input(id='flash_cmax', value='max', type='str'),# style={'padding':'5px'}),
                  html.Label(' (type "min" and "max" for automatic range setting)', style=para_text_style),
                 ], style={'width': '100%', 'display': 'inline-block', 'padding': '5px'}),
        html.Div(id='error_flash_range'),
        html.Div([html.Label('X position ', style=label_text_style),
                  #dcc.Input(id='xoffset_text', value='min', type='str')
                 ],style={'padding': '5px'}),
        html.Div([dcc.Slider(id='xoffset',min=xmin,max=xmax,step=.1,value=xmin),
                 ],style={'padding': '5px'}),
        #dcc.Loading(id="loading-playdata",
        #            children=[html.Div([dcc.Graph(id='playdata',figure=_manager.empty_view)],
        #                               style={'padding': '5px'})],
        #            type="circle"),
        html.Div(id="match_result"),
        dcc.Loading(id="loading_dynamic_display", children=[html.Div(id="wait_dynamic_display")], type="default"),
        html.Div([dcc.Graph(id='playdata',figure=_manager.empty_view)],style={'padding': '5px'}),
    ],style={'width': '80%', 'display': 'inline-block', 'vertical-align': 'middle'})

    #
    # call backs
    #
    @app.callback(dash.dependencies.Output("error_entry", "children"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def error_entry(entry_or_event,data_index):
        is_entry = entry_or_event == 'entry'
        if data_index is None:
            raise dash.exceptions.PreventUpdate
        if not data_index.isdigit():
            return [html.Label("Invalid %s ID %s (not an integer)" % (entry_or_event,data_index),
                               style={'color':'red'})]
        data_index = int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            return [html.Label("Invalid %s ID %d (data does not exist)" % (entry_or_event,data_index),
                               style={'color':'red'})]
        else:
            return ""

    @app.callback(dash.dependencies.Output("error_pmt_range", "children"),
                  [dash.dependencies.Input("pmt_cmin","value"),
                   dash.dependencies.Input("pmt_cmax","value")])
    def error_entry(pmt_cmin,pmt_cmax):
        pmt_cmin = str(pmt_cmin).lower()
        pmt_cmax = str(pmt_cmax).lower()
        cmin_good,cmax_good=True,True
        try:
            if not pmt_cmin == "min":
                pmt_cmin = float(pmt_cmin)
        except ValueError:
            return [html.Label('Invalid PMT min value: %s (either numerical value or "min")' % pmt_cmin,
                                style={'color':'red'})]
        try:
            if not pmt_cmax == "max":
                pmt_cmax = float(pmt_cmax)
        except ValueError:
            return [html.Label('Invalid PMT max value: %s (either numerical value or "min")' % pmt_cmax,
                                style={'color':'red'})]
        return ""

    @app.callback(dash.dependencies.Output("error_flash_range", "children"),
              [dash.dependencies.Input("flash_cmin","value"),
               dash.dependencies.Input("flash_cmax","value")])
    def error_entry(flash_cmin,flash_cmax):
        flash_cmin = str(flash_cmin).lower()
        flash_cmax = str(flash_cmax).lower()
        cmin_good,cmax_good=True,True
        try:
            if not flash_cmin == "min":
                flash_cmin = float(flash_cmin)
        except ValueError:
            return [html.Label('Invalid PMT min value: %s (either numerical value or "min")' % flash_cmin,
                                style={'color':'red'})]
        try:
            if not flash_cmax == "max":
                flash_cmax = float(flash_cmax)
        except ValueError:
            return [html.Label('Invalid PMT max value: %s (either numerical value or "min")' % flash_cmax,
                                style={'color':'red'})]
        return ""


    @app.callback(dash.dependencies.Output("wait_static_display", "children"),
                  [dash.dependencies.Input("mode_flash", "value")])
    def wait_static_display(mode):
        if mode == 'hypothesis':
            from flashmatch import phot
            print('hypothesis mode chosen, loading photon library...')
            phot.PhotonVisibilityService.GetME().LoadLibrary()
        return ""

    @app.callback(dash.dependencies.Output("wait_dynamic_display", "children"),
                  [dash.dependencies.Input("target_qcluster", "value")])
    def wait_dynamic_display(target_qcluster):
        if target_qcluster is not None:
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
        if data_index is None:
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
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
        if data_index is None:
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_qcluster(data_index=int(data_index),is_entry=is_entry)
        if options is None: raise dash.exceptions.PreventUpdate
        else: return options[:-1]


    @app.callback(dash.dependencies.Output("select_flash","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value"),
                   dash.dependencies.Input("mode_flash","value")])
    def update_dropdown_select_flash(entry_or_event,data_index,mode_flash):
        """
        Update drop-down "select_flash" options for Flash when event/entry is changed
        """
        is_entry = entry_or_event == 'entry'
        mode_flash = mode_flash == 'flash'
        if data_index is None:
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_flash(data_index=int(data_index),
                                          is_entry=is_entry,
                                          mode_flash=mode_flash)
        if options is None: raise dash.exceptions.PreventUpdate
        else: return options

    @app.callback(dash.dependencies.Output("target_flash","options"),
                  [dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")])
    def update_dropdown_target_flash(entry_or_event,data_index):
        """
        Update drop-down "target_flash" options for Flash when event/entry is changed
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None:
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        options = _manager.dropdown_flash(data_index=int(data_index),
                                          is_entry=is_entry,
                                          mode_flash=True)
        if options is None: raise dash.exceptions.PreventUpdate
        else: return options[:-1]

    @app.callback(dash.dependencies.Output("visdata","figure"),
                  [dash.dependencies.Input("select_qcluster","value"),
                   dash.dependencies.Input("select_flash","value"),
                   dash.dependencies.Input("mode_flash","value"),
                   dash.dependencies.Input("mode_qcluster","value"),
                   dash.dependencies.Input("use_all_pts","value"),
                   dash.dependencies.Input("pmt_cmin","value"),
                   dash.dependencies.Input("pmt_cmax","value"),
                   dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")]
                 )
    def update_static(select_qcluster, select_flash, mode_flash, mode_qcluster, use_all_pts,
                      pmt_cmin, pmt_cmax,
                      entry_or_event, data_index):
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
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        pmt_range=[0,0]
        if pmt_cmin.lower()=='min':
            pmt_range[0]='min'
        else:
            try:
                pmt_range[0]=str(float(pmt_cmin))
            except ValueError:
                raise dash.exceptions.PreventUpdate
        if pmt_cmax.lower()=='max': pmt_range[1]='max'
        else:
            try:
                pmt_range[1]=str(float(pmt_cmax))
            except ValueError:
                raise dash.exceptions.PreventUpdate
        fig = _manager.event_display(int(data_index), select_qcluster, select_flash,
                                     is_entry, mode_flash, mode_qcluster, use_all_pts, pmt_range)
        if fig is None: raise dash.exceptions.PreventUpdate
        else: return fig


    @app.callback(dash.dependencies.Output("playdata","figure"),
                  [dash.dependencies.Input("target_qcluster","value"),
                   dash.dependencies.Input("target_flash","value"),
                   dash.dependencies.Input("xoffset","value"),
                   dash.dependencies.Input("flash_cmin","value"),
                   dash.dependencies.Input("flash_cmax","value"),
                   dash.dependencies.Input("normalize_pe","value"),
                   dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value")]
                 )
    def update_dynamic(target_qcluster, target_flash, xoffset,
                       flash_cmin, flash_cmax, normalize_pe,
                       entry_or_event, data_index):
        """
        Update "playdata" 3D display for change in event/entry and/or selection of flash/qcluster
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None or (target_qcluster is None and target_flash is None):
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        flash_range=[0,0]
        if flash_cmin.lower()=='min': flash_range[0]='min'
        else:
            try:
                flash_range[0]=str(float(flash_cmin))
            except ValueError:
                raise dash.exceptions.PreventUpdate
        if flash_cmax.lower()=='max': flash_range[1]='max'
        else:
            try:
                flash_range[1]=str(float(flash_cmax))
            except ValueError:
                raise dash.exceptions.PreventUpdate
        fig = _manager.hypothesis_display(int(data_index), target_qcluster, target_flash,
                                          flash_range, int(normalize_pe),
                                          is_entry, float(xoffset))
        if fig is None: raise dash.exceptions.PreventUpdate
        else: return fig

    @app.callback(dash.dependencies.Output("match_result","children"),
                  [dash.dependencies.Input("target_qcluster","value"),
                   dash.dependencies.Input("target_flash","value"),
                   dash.dependencies.Input("xoffset","value"),
                   dash.dependencies.Input("entry_or_event","value"),
                   dash.dependencies.Input("data_index","value"),]
                 )
    def update_match(target_qcluster, target_flash, xoffset, entry_or_event, data_index):
        """
        Update "playdata" 3D display for change in event/entry and/or selection of flash/qcluster
        """
        is_entry = entry_or_event == 'entry'
        if data_index is None or target_qcluster is None or target_flash is None:
            raise dash.exceptions.PreventUpdate
        data_index=int(data_index)
        if _manager.valid_data_entry(data_index,is_entry) is None:
            raise dash.exceptions.PreventUpdate
        qll = _manager.run_qll(int(data_index), target_qcluster, target_flash,
                               is_entry, float(xoffset))
        show_style_text = dict(fontFamily='Georgia', fontWeight=300, fontSize=20)
        show_style_num  = dict(fontFamily='Georgia', fontWeight=600, fontSize=20, color='blue')

        xoffset = float(xoffset)
        res = []
        if qll is None:
            raise dash.exceptions.PreventUpdate
        else:
            # xmin and score_xoffset should not be None
            xmin,score_xoffset,match = qll
            res += [html.Label('QLL joint probability ',style=show_style_text),
                    html.Label('%.4f' % score_xoffset,style=show_style_num),
                    html.Label(' at X ',style=show_style_text),
                    html.Label('%.2f' % xoffset,style=show_style_num),
                    html.Label(' ... true X ',style=show_style_text),
                    html.Label('%.2f' % xmin,style=show_style_num),
                    html.Br()]
            if match is not None:
                res += [html.Label('QLL match result: score ',style=show_style_text),
                        html.Label('%.4f' % match.score, style=show_style_num),
                        html.Label(' at X ',style=show_style_text),
                        html.Label('%.2f' % match.tpc_point.x,style=show_style_num),
                        html.Label(' ... %.1fms %dsteps (scanned %.2f to %.2f)' % (match.duration/1.e6,
                                                                                   match.num_steps,
                                                                                   match.minimizer_min_x,
                                                                                   match.minimizer_max_x))]

        return res

    app.server.run(port=port)
