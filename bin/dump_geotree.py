from ROOT import TChain
import sys
ch=TChain("geotree")
ch.AddFile(sys.argv[1])
ch.GetEntry(0)

pmtX,pmtY,pmtZ = ch.pmtX, ch.pmtY, ch.pmtZ
minX,minY,minZ = ch.minX, ch.minY, ch.minZ
maxX,maxY,maxZ = ch.maxX, ch.maxY, ch.maxZ

tpcid=0

print('DetectorSpecs: {')
print('  MaxPosition: [%f,%f,%f]' % (maxX[tpcid],maxY[tpcid],maxZ[tpcid]))
print('  MinPosition: [%f,%f,%f]' % (minX[tpcid],minY[tpcid],minZ[tpcid]))
for pmt in range(pmtX.size()):
    print('  PMT%d: [%f,%f,%f]' % (pmt, pmtX[pmt], pmtY[pmt], pmtZ[pmt]))
print('}')
