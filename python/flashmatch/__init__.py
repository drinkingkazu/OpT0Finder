from ROOT import flashmatch, phot, sim, geoalgo
#force loading C functions in dict by instantiating a class
c=flashmatch.PSet
c=flashmatch.load_pyutil
#from .demo import demo
from .toymc import ToyMC
from .rootinput import ROOTInput
from .utils import AnalysisManager
