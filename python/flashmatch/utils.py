from flashmatch import flashmatch
import numpy as np

class ManagerAPI:
    """
    FlashMatchManager simple wrapper
    """
    def __init__(self,cfg=None):
        self.mgr=flashmatch.FlashMatchManager()
        self.cfg=None
        if cfg is not None:
            self.configure(cfg)
            
    def configure(self,cfg):
        """
        Configure FlashMatchManager
        """
        self.cfg = flashmatch.CreatePSetFromFile(cfg)
        self.mgr.Configure(self.cfg)

    def dump_config(self):
        """
        Dump configuration
        """
        if self.cfg is None: return
        return self.cfg.dump()

    def FlashHypothesis(self,qcluster):
        """
        Fill flash hypothesis for a provided QCluster
        """
        flash_algo = self.mgr.GetAlgo(flashmatch.kFlashHypothesis)
        flash = flashmatch.Flash_t()
        flash_algo.FillEstimate(qcluster, flash)
        return flash

    def QCluster(self,xs,ys,zs):
        """
        If LightPath algorithm is available, use it to generate QCluster_t
        """
        xyzs = np.column_stack([xs,ys,zs]).astype(np.float32)
        traj = as_geoalgo_trajectory(xyzs)
        qcluster_algo = self.mgr.GetCustomAlgo("LightPath")
        return qcluster_algo.MakeQCluster(traj)

    def OffsetQCluster(self, qcluster, xoffset):
        xoffset = qcluster.min_x() - xoffset - flashmatch.DetectorSpecs.GetME().ActiveVolume().Min()[0]
        qcluster = qcluster - xoffset
        return qcluster
    
class FlashMatchInput:
    def __init__(self):
        # input array of flashmatch::Flash_t
        self.flash_v = []
        # input array of flashmatch::QCluster_t
        self.qcluster_v = []
        # "RAW" flashmatch::QCluster_t (optional, may not be present)
        self.raw_qcluster_v = []
        # True matches, an array of integer-pairs.
        self.true_match = []
        
def print_history(self,mgr):
    """
    Print history of the last minimization if stored
    """
    qll=mgr.GetAlgo(flashmatch.kFlashMatch)
    hist_llhd = qll.HistoryLLHD()
    hist_chi2 = qll.HistoryChi2()
    hist_xpos = qll.HistoryX()
    if hist_llhd.size():
        print('QLLMatch history')
        #np.savetxt('hist_llhd.csv', np.array(hist_llhd), delimiter=',')
        #np.savetxt('hist_chi2.csv', np.array(hist_chi2), delimiter=',')
        #np.savetxt('hist_xpos.csv', np.array(hist_xpos), delimiter=',')
        for i in range(hist_xpos.size()):
            print('Step %d X %f Chi2 %f LLHD %f' % (i,hist_xpos[i],hist_chi2[i],hist_llhd[i]))



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

