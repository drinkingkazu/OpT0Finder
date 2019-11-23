import sys,os
import argparse
try:
    import flashmatch.visualization as vis
except ImportError:
    sys.stderr.write('failed to import visualization library (plotly or dash missing?)\n')
    sys.stderr.flush()
import flashmatch.visualization as vis

parser = argparse.ArgumentParser(description="view_data.py argument parser")
#subparsers = self.parser.add_subparsers(title="Modules", description="Valid subcommands", dest='script', help="aho")

# train parser
if not 'FMATCH_DIR' in os.environ:
    sys.stderr.write('FMATCH_DIR not in your shell environment variable (did you run configure.sh?)\n')
    sys.stderr.flush()
    sys.exit(1)

cfg = os.path.join(os.environ['FMATCH_BASEDIR'],'dat/flashmatch.cfg')
geo = os.path.join(os.environ['FMATCH_BASEDIR'],'dat/detector_specs.cfg' )
print(cfg)
parser.add_argument('-o','--opflash_data', type=str,
                          help='OpFlash data file (ICARUSOpFlashAna_module output)')
parser.add_argument('-p','--particle_data', type=str,
                          help='Particle data file (ICARUSParticleAna_module output)')
parser.add_argument('-c','--cfg', type=str, default=cfg,
                          help='Flash matching (OpT0Finder) PSet config file [default: %s]' % cfg)
parser.add_argument('-g','--geo_cfg', type=str, default=geo,
                          help='Detector geometry data yaml/PSet data file [default: %s]' % geo)

args=parser.parse_args()

data_particle = args.particle_data
data_opflash = args.opflash_data
if data_particle is None:
    sys.stderr.write('-p --particle_data (particle data file) argument is required!\n')
    sys.stderr.flush()
    sys.exit(1)
if data_opflash is None:
    sys.stderr.write('-o --opflash_data (opflash data file) argument is required!\n')
    sys.stderr.flush()
    sys.exit(1)
vis.view_data(cfg,geo,data_particle,data_opflash)
sys.exit(0)
