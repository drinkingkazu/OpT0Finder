from flashmatch import flashmatch
import numpy as np

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


class AnalysisManager(object):
    """
    FlashMatchManager simple wrapper
    """
    def __init__(self, cfg, particleana='', opflashana=''):
        self.configure(cfg, particleana, opflashana)


    def configure(self,cfg,particleana='',opflashana=''):
        """
        Configure FlashMatchManager
        """
        import os
        self.cfg = flashmatch.CreatePSetFromFile(cfg)
        self.matcher=flashmatch.FlashMatchManager()
        self.matcher.Configure(self.cfg)
        if not particleana and not opflashana:
            from .toymc import ToyMC
            self.reader=ToyMC(cfg)
        else:
            from .rootinput import ROOTInput
            if not os.path.isfile(particleana):
                print('File not found:',particleana)
                raise FileNotFoundError
            if not os.path.isfile(opflashana):
                print('File not found:',opflashana)
                raise FileNotFoundError
            self.reader=ROOTInput(particleana,opflashana,cfg)


    def dump_config(self):
        """
        Dump configuration
        """
        if self.cfg is None: return
        return self.cfg.dump()


    def make_flashmatch_input(self, entry):
        """
        Create flash matching input
        INPUT:
          - entry ... integer key to generate data (a sequential data entry index
            for ROOTInput, number of trajectory to generate for ToyMC)
        RETURN:
          - FlashMatchInput class instance
        """
        return self.reader.make_flashmatch_input(entry)


    def event_id(self,entry_id):
        """
        Convert "entry id" to "event id"
        """
        from .rootinput import ROOTInput
        if not isinstance(self.reader, ROOTInput):
            return -1
        return self.reader.event_id(entry_id)


    def entry_id(self,event_id):
        """
        Convert "event id" to "entry id"
        """
        from .rootinput import ROOTInput
        if not isinstance(self.reader, ROOTInput):
            return -1
        return self.reader.entry_id(event_id)


    def entries(self):
        """
        Returns a list of "entry ids" (only for ROOTInput)
        """
        from .rootinput import ROOTInput
        if not isinstance(self.reader, ROOTInput):
            return []
        return np.arange(len(self.reader))


    def flash_hypothesis(self,qcluster):
        """
        Fill flash hypothesis for a provided QCluster
        """
        flash_algo = self.matcher.GetAlgo(flashmatch.kFlashHypothesis)
        flash = flashmatch.Flash_t()
        flash_algo.FillEstimate(qcluster, flash)
        return flash


    def offset_qcluster(self, qcluster, xoffset):
        """
        Shifts a qcluster by the specified offset, returns a new qcluster
        """
        xoffset = qcluster.min_x() - xoffset - flashmatch.DetectorSpecs.GetME().ActiveVolume().Min()[0]
        qcluster = qcluster - xoffset
        return qcluster


    def run_flash_match(self,match_input):
        """
        Runs flash matching, the input type should be FlashMatchInput
        """
        # Register for matching
        self.matcher.Reset()
        for pmt in match_input.flash_v:
            self.matcher.Add(pmt)
        for tpc in match_input.qcluster_v:
            self.matcher.Add(tpc)
        # Run matching
        return self.matcher.Match()


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
