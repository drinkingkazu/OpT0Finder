from flashmatch import flashmatch,AnalysisManager
import numpy as np

def attribute_names():

    return [
        'event',
        'entry',
        'score',
        'touch_match',
        'flash_idx',
        'track_idx',
        'true_min_x',
        'true_max_x',
        'qcluster_min_x',
        'qcluster_max_x',
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
        'flash_max_pe',
        'flash_min_pe',
        'flash_max_pe_true',
        'flash_min_pe_true',
        'hypothesis_max_pe',
        'hypothesis_min_pe',
        'duration',
        'num_steps',
        'minimizer_min_x',
        'minimizer_max_x'
    ]

def demo(cfg_file,repeat=1,num_tracks=None,out_file='',particleana=None,opflashana=None,start_entry=0,num_entries=-1):
    """
    Run function for ToyMC
    ---------
    Arguments
      cfg_file:   string for the location of a config file
      repeat:     int, number of times to run the toy MC simulation
      out_file:   string for an output analysis csv file path (optional)
      num_tracks: int for number of tracks to be generated (optional)
    """
    # create & configure manager API
    mgr = AnalysisManager(cfg_file,particleana,opflashana)
    cfg = flashmatch.CreatePSetFromFile(cfg_file)
    # dump config
    sys.stdout.write(mgr.dump_config())
    sys.stdout.flush()
    # override number of tracks to simulate
    if num_tracks is not None:
        num_tracks = int(num_tracks)

    entries=mgr.entries()

    toymc=False
    if len(entries) < 1:
        toymc=True
        entries = np.arange(repeat)
    else: 

        if start_entry <0:
            return
        else:
            entries = np.arange(start_entry,max(entries)+1)
            if len(entries) < 1:
                return 

        if num_entries >0:
            entries = np.arange(start_entry,min(start_entry+num_entries,max(entries)+1))
            if len(entries)<1:
                return

    if out_file:
        import os
        if os.path.isfile(out_file):
            print('Output file',out_file,'already exists. Exiting...')
            return 

    np_result = None
    counter = 0
    #entries = [14]
    for entry in entries:
        sys.stdout.write('Entry %d/%d\n' %(entry,len(entries)))
        sys.stdout.write('Event %d\n' % mgr.event_id(entry))
        sys.stdout.flush()
        # Generate samples
        generator_arg = entry if not toymc else num_tracks
        print(generator_arg)
        match_input = mgr.make_flashmatch_input(generator_arg)
        match_v = mgr.run_flash_match(match_input)

        # If no analysis output saving option given, return
        if not out_file:
            print('Number of match result',len(match_v))

        all_matches = []
        for idx, match in enumerate(match_v):
            tpc = match_input.qcluster_v[match.tpc_id]
            pmt = match_input.flash_v[match.flash_id]
            raw = match_input.raw_qcluster_v[match.tpc_id]
            truncation   = raw.front().dist(raw.back()) - tpc.front().dist(tpc.back())
            if raw.front().dist(raw.back()) > 0:
                truncation_frac = truncation / raw.front().dist(raw.back())
            else:
                truncation_frac = -1
            #if particleana is None or opflashana is None:
            true_minx = raw.min_x()  #true_xmin_v[match.flash_id]
            true_maxx = raw.max_x()
            # print(true_minx, raw.min_x())
            # if (particleana is None and opflashana is None) or (hasattr(ihandler, '_shift_tpc') and ihandler._shift_tpc):
            #     true_minx -= tpc.time_true * ihandler.det.DriftVelocity()
            #     true_maxx -= tpc.time_true * ihandler.det.DriftVelocity()
            #print(true_minx)
            reco_maxx = match.tpc_point.x + (tpc.max_x() - tpc.min_x())
            correct_match = (pmt.idx,tpc.idx) in match_input.true_match

            if not out_file:
                print('Match ID',idx)

                msg = '  TPC/PMT IDs %d/%d Correct? %s Touch match? %d F-Score %.3f T-Score %.3f Trunc. %.3f ... reco vs. true: X %.3f vs. (%.3f, %.3f), PE %f vs. %f, Time %f vs. %f'
                msg = msg % (tpc.idx, pmt.idx, correct_match, match.touch_match, match.score, match.touch_score, truncation, match.tpc_point.x, true_minx, true_maxx,
                             np.sum(match.hypothesis), np.sum(pmt.pe_v), tpc.time_true, pmt.time_true)
                print(msg)
                continue
            #print(match.hypothesis)
            store = np.array([[
                mgr.event_id(entry),
                entry,
                match.score,
                match.touch_match,
                pmt.idx,
                tpc.idx,
                true_minx,
                true_maxx,
                tpc.min_x(),
                tpc.max_x(),
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
                int(correct_match),
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
                pmt.max_pe(),
                pmt.min_pe(),
                pmt.max_pe_true(),
                pmt.min_pe_true(),
                np.max(match.hypothesis),
                np.min(match.hypothesis),
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

        if entry%100 == 0 and out_file:

            with open(out_file,'a') as f:

                if entry < 1:
                    np.savetxt(f,np_result,delimiter=',',header=','.join(attribute_names()))
                else:
                    np.savetxt(f,np_result,delimiter=',')

            np_result=None

        counter += len(match_input.true_match)
        #break

    print("counter ", counter)
    if not out_file:
        return
    #print(x.shape)

    np.savetxt(out_file, np_result, delimiter=',', header=','.join(attribute_names()))

if __name__ == '__main__':
    import os,sys

    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                            'dat', 'flashmatch.cfg')
    num_tracks = 5
    particleana, opflashana = None, None
    start_entry = 0
    num_entries = -1
    outfile = ''
    if len(sys.argv) > 1:
        if len(sys.argv) == 2:
            num_tracks = int(sys.argv[1])
        else:
            num_tracks = None
            particleana = sys.argv[1]
            opflashana = sys.argv[2]
            if len(sys.argv) > 3:
                for v in sys.argv[3:]:
                    if v.startswith('out='):
                        outfile = v.replace('out=','')
                    if v.startswith('start='):
                        start_entry = int(v.replace('start=',''))
                    if v.startswith('num='):
                        num_entries = int(v.replace('num=',''))
                    if v.startswith('end='):
                        num_entries = int(v.replace('end=','')) - start_entry + 1
                        if num_entries <= 0:
                            print('start=%d and end =%d cannot run' % (start_entry,int(v.replace('end=',''))))
                            sys.exit(1)

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
        demo(cfg_file, particleana=particleana, opflashana=opflashana, out_file=outfile, start_entry=start_entry, num_entries=num_entries)
