import numpy as np
from ROOT import flashmatch

def dump_geo():
    '''
    A function to print out geometry used in OpT0Finder package using flashana::DetectorSpecs
    '''
    s = flashmatch.DetectorSpecs.GetME()
    xyz_range = [float()] * 6
    bbox=s.ActiveVolume()
    print('Active Volume')
    print('X:',bbox.Max()[0],bbox.Min()[0])
    print('Y:',bbox.Max()[1],bbox.Min()[1])
    print('Z:',bbox.Max()[2],bbox.Min()[2])

    print('Visibility @ Center per PMT')
    cx=(bbox.Max()[0]-bbox.Min()[0])/2. + bbox.Min()[0]
    cy=(bbox.Max()[1]-bbox.Min()[1])/2. + bbox.Min()[1]
    cz=(bbox.Max()[2]-bbox.Min()[2])/2. + bbox.Min()[2]
    for opch in range(s.NOpDets()):
        print('PMT',opch,'Visibility',s.GetVisibility(cx,cy,cz,opch))
    
dump_geo()
