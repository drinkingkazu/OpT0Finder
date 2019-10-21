import numpy as np
from ROOT import flashmatch

def dump_geo():
    '''
    A function to print out geometry used in OpT0Finder package using flashana::DetectorSpecs
    '''
    s = flashmatch.DetectorSpecs.GetME()
    xyz_range = [float()] * 6
    s.MaxPosition(xyz_range[0],xyz_range[2],xyz_range[4])
    s.MinPosition(xyz_range[1],xyz_range[3],xyz_range[5])
    print('Active Volume')
    print('X:',xyz_range[0:2])
    print('Y:',xyz_range[2:4])
    print('Z:',xyz_range[4:6])
    
dump_geo()
