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
if len(sys.argv) > 1:
    for argv in sys.argv[1:]:
        if argv.startswith('repeat='):
            repeat = int(argv.replace('repeat=',''))
        elif argv.startswith('cfg='):
            cfg_file = argv.replace('cfg=','')
        elif argv.startswith('out='):
            out_file = argv.replace('out=','')

# run demo
toymc.demo(cfg_file,repeat=repeat,out_file=out_file)
