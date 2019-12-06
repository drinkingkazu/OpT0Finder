import sys,os
import argparse

parser = argparse.ArgumentParser(description="view_data.py is a simple event display executer using plotly/dash on your web-browser.")
#subparsers = self.parser.add_subparsers(title="Modules", description="Valid subcommands", dest='script', help="aho")

# train parser
if not 'FMATCH_DIR' in os.environ:
    sys.stderr.write('FMATCH_DIR not in your shell environment variable (did you run configure.sh?)\n')
    sys.stderr.flush()
    sys.exit(1)

cfg = os.path.join(os.environ['FMATCH_BASEDIR'],'dat/flashmatch.cfg')
geo = os.path.join(os.environ['FMATCH_BASEDIR'],'dat/detector_specs.cfg' )
dark_mode = True
port= 5000

parser.add_argument('-od','--opflash_data', type=str,
                    help='OpFlash data file (ICARUSOpFlashAna_module output)')
parser.add_argument('-pd','--particle_data', type=str,
                    help='Particle data file (ICARUSParticleAna_module output)')
parser.add_argument('-p','--port',type=int, default=port,
                    help='Port id number [default %s]' % port)
parser.add_argument('-c','--cfg', type=str, default=cfg,
                    help='Flash matching (OpT0Finder) PSet config file [default: %s]' % cfg)
parser.add_argument('-g','--geo_cfg', type=str, default=geo,
                    help='Detector geometry data yaml/PSet data file [default: %s]' % geo)
parser.add_argument('-d','--dark_mode', type=str, default=str(dark_mode),
                    help='Dark mode in plotting [default: %s]' % dark_mode)

args=parser.parse_args()
import ast
cfg = args.cfg
geo = args.geo_cfg
data_particle = args.particle_data
data_opflash = args.opflash_data
dark_mode = bool(ast.literal_eval(args.dark_mode))
port = int(args.port)
if data_particle is None:
    sys.stderr.write('-p --particle_data (particle data file) argument is required!\n')
    sys.stderr.flush()
    sys.exit(1)
if data_opflash is None:
    sys.stderr.write('-o --opflash_data (opflash data file) argument is required!\n')
    sys.stderr.flush()
    sys.exit(1)
print('Using OpT0Finder config  :',cfg)
print('Using geometry config    :',cfg)
print('Using Particle data file :',data_particle)
print('Using OpFlash data file  :',data_opflash)
print('Using a port ID:',port)
print('Dark mode:',dark_mode)

#
# Import flashmatch stuff AFTER parsing, otherwise argparse is overriden by ROOT argument parsing
# (problematic for --help)
#
try:
    import flashmatch.visualization as vis
except ImportError:
    sys.stderr.write('failed to import visualization library (plotly or dash missing?)\n')
    sys.stderr.flush()
import flashmatch.visualization as vis

vis.view_data(cfg,geo,data_particle,data_opflash,dark_mode,port=port)
sys.exit(0)
