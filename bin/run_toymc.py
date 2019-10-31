#!/usr/bin/env python
import os
import sys
import yaml
from os import environ

current_directory = os.path.dirname(os.path.abspath(__file__))
current_directory = os.path.dirname(current_directory)
sys.path.insert(0, current_directory)
from flashmatch import flashmatch, toymc



cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],
                        'dat', 'flashmatch.cfg')
out_file = ''
num_tracks = None
if len(sys.argv) > 1:
    for argv in sys.argv[1:]:
        if argv.isdigit():
            num_tracks = int(argv)
        elif argv.endswith('.cfg'):
            cfg_file = argv
        elif argv.endswith('.csv'):
            out_file = argv

# run demo
toymc.demo(cfg_file,num_tracks,out_file)
