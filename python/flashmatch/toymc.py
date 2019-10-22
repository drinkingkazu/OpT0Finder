import numpy as np
from flashmatch import flashmatch,geoalgo
import numpy as np
import sys

# dumb function to generate a stupid trajectory
def gen_trajectories(num_tracks,det):
    res = []
    for i in range(num_tracks):
        trj = geoalgo.Trajectory(2,3) # 2 points, 3d
        trj[0] = det.ActiveVolume().Min()
        trj[1] = det.ActiveVolume().Max()
        res.append(trj)
    return res


# Create a set ot TPC (flashmatch::QCluster) and PMT (flashmatch::Flash_t) interactions
def gen_input(num_tracks,blob):
    track_v = gen_trajectories(num_tracks,blob.det)
    tpc_v = flashmatch.QClusterArray_t()
    pmt_v = flashmatch.FlashArray_t()
    
    tpc_algo = blob.mgr.GetCustomAlgo("LightPath")
    pmt_algo = blob.mgr.GetAlgo(flashmatch.kFlashHypothesis)
    for idx,track in enumerate(track_v):
        # Decide time range to generate, completely made up for now
        tmin,tmax=(-1000,1000) # in micro-seconds
        time = np.random.random() * (tmax - tmin) + tmin
        # TPC xyz & light info
        qcluster = tpc_algo.FlashHypothesis(track)
        qcluster.idx = idx
        # PMT pe spectrum
        flash = flashmatch.Flash_t()
        flash.idx  = idx
        flash.time = time
        pmt_algo.FillEstimate(qcluster,flash)
        # save
        tpc_v.push_back(qcluster)
        pmt_v.push_back(flash)
        
    return tpc_v,pmt_v


def demo(cfg_file,num_tracks):

    class BLOB:
        pass

    blob=BLOB()
    blob.mgr=flashmatch.FlashMatchManager()
    blob.det=flashmatch.DetectorSpecs.GetME()
    blob.cfg=flashmatch.CreatePSetFromFile(cfg_file)

    # configure
    blob.mgr.Configure(blob.cfg)
    return gen_input(num_tracks, blob)

if __name__ == '__main__':
    import os
    cfg_file = os.path.join(os.environ['FMATCH_BASEDIR'],'dat','flashmatch.cfg')
    if len(sys.argv)>1:
        cfg_file=sys.argv[1]
    
    demo(cfg_file,10)

