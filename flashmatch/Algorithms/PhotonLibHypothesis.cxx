#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"
#include <assert.h>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include <omp.h>
#define NUM_THREADS 12
#endif

//using namespace std::chrono;

namespace flashmatch {

    static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

    PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
    {
        #if USING_LARSOFT == 0
        omp_set_num_threads(NUM_THREADS); 
        #endif
    }

    void PhotonLibHypothesis::_Configure_(const Config_t &pset)
    {

        _global_qe = pset.get<double>("GlobalQE");
        _global_qe_refl = pset.get<double>("GlobalQERefl", -1);
        _extend_tracks = pset.get<bool>("ExtendTracks", false);
        _threshold_extend_track = pset.get<double>("ThresholdExtendTrack", 5.0);
        _segment_size = pset.get<double>("SegmentSize", 0.5);

        _qe_v.clear();
        _qe_v = pset.get<std::vector<double> >("CCVCorrection",_qe_v);
        if(_qe_v.empty()) _qe_v.resize(DetectorSpecs::GetME().NOpDets(),1.0);
        if(_qe_v.size() != DetectorSpecs::GetME().NOpDets()) {
          FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
          << " != number of opdet (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
          throw OpT0FinderException();
        }
    }

    bool PhotonLibHypothesis::InspectTouchingEdges(const QCluster_t& trk, size_t& min_idx, size_t& max_idx) const
    {
        double min_x = kINVALID_DOUBLE; double max_x = -kINVALID_DOUBLE;
        min_idx = max_idx = 0;
        for (size_t pt_index = 0; pt_index < trk.size(); ++pt_index) {
            if (trk[pt_index].x < min_x) {
                min_x = trk[pt_index].x;
                min_idx = pt_index;
            }
            if (trk[pt_index].x > max_x) {
                max_x = trk[pt_index].x;
                max_idx = pt_index;
            }
        } 
        if (((trk[min_idx].x - DetectorSpecs::GetME().ActiveVolume().Min()[0]) < _threshold_extend_track) ||
            ((trk[min_idx].y - DetectorSpecs::GetME().ActiveVolume().Min()[1]) < _threshold_extend_track) ||
            ((trk[min_idx].z - DetectorSpecs::GetME().ActiveVolume().Min()[2]) < _threshold_extend_track) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[0] - trk[max_idx].x) < _threshold_extend_track) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[1] - trk[max_idx].y) < _threshold_extend_track) ||
            ((DetectorSpecs::GetME().ActiveVolume().Max()[2] - trk[max_idx].z) < _threshold_extend_track)) {
                //std::cout << "*** Extending track " << trk.size() << " " << min_x << " " << max_x << std::endl;
                //std::cout << trk.front().x << " " << trk.back().x << std::endl;
            return true;
        }
        return false;
    }

    QCluster_t PhotonLibHypothesis::TrackExtension(const QCluster_t& in_trk, const size_t min_idx, const size_t max_idx) const
    {
        QCluster_t trk = in_trk;
        // Extend the track
        // Compute coordinates of final point first
        geoalgo::Vector A(trk[max_idx].x, trk[max_idx].y, trk[max_idx].z);
        geoalgo::Vector B(trk[min_idx].x, trk[min_idx].y, trk[min_idx].z);
        geoalgo::Vector AB = B - A;

        double x_C = DetectorSpecs::GetME().PhotonLibraryVolume().Min()[0];
        double lengthBC = (x_C - B[0])/(B[0] - A[0]) * AB.Length();
        geoalgo::Vector C = B + AB / AB.Length() * lengthBC;

                // Add to _var_trk the part betwen boundary and C
                //_custom_algo->MakeQCluster(C, B, _var_trk, -1);
        geoalgo::Vector unit = (C - B).Dir();

        QPoint_t q_pt;
        geoalgo::Vector current = B;
        int num_pts = int(lengthBC / _segment_size);
        trk.reserve(trk.size() + num_pts);
        for (int i = 0; i < num_pts+1; i++) {
            double current_segment_size = (i < num_pts ? _segment_size : (lengthBC - _segment_size*num_pts));
            current = current + unit * current_segment_size/2.0;
            q_pt.x = current[0];
            q_pt.y = current[1];
            q_pt.z = current[2];
            q_pt.q = current_segment_size * DetectorSpecs::GetME().LightYield() * DetectorSpecs::GetME().MIPdEdx();
            if (trk.front().x < trk.back().x) {
                trk.insert(trk.begin(), q_pt);
            }
            else {
                trk.emplace_back(q_pt);
            }
            current = current + unit * current_segment_size/2.0;
            //std::cout << "Adding point " << current  << " " << i << " " << num_pts << std::endl;
        }
        //std::cout << " done " << trk.size() << std::endl;
        //std::cout << trk.front().x << " " << trk.back().x << std::endl;

        return trk;
    }


    void PhotonLibHypothesis::FillEstimate(const QCluster_t& tpc_trk, Flash_t &flash) const
    {
        size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
        if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
        if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);
        if(flash.pe_true_v.empty()) flash.pe_true_v.resize(n_pmt);

        assert(flash.pe_v.size()     == n_pmt);
        assert(flash.pe_true_v.size() == n_pmt);
        assert(flash.pe_err_v.size() == n_pmt);

        for (auto& v : flash.pe_v      ) {v = 0;}
        for (auto& v : flash.pe_err_v  ) {v = 0;}
        for (auto& v : flash.pe_true_v ) {v = 0;}

        size_t min_idx, max_idx;
        bool extend_tracks = (_extend_tracks && this->InspectTouchingEdges(tpc_trk,min_idx,max_idx));

        if(!extend_tracks) 
            this->BuildHypothesis(tpc_trk,flash);
        else {
            auto trk = this->TrackExtension(tpc_trk,min_idx,max_idx);
            this->BuildHypothesis(trk,flash);
        }

    }

    void PhotonLibHypothesis::BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const
    {

        size_t n_pmt = DetectorSpecs::GetME().NOpDets();
        //auto const& lib_data = DetectorSpecs::GetME().GetPhotonLibraryData();
        
        #if USING_LARSOFT == 0
        #pragma omp parallel
        #endif
        {
            size_t num_pts  = trk.size();
            size_t start_pt = 0;

            #if USING_LARSOFT == 0
            size_t thread_id = omp_get_thread_num();
            size_t num_threads = omp_get_num_threads();
            if(num_threads == 1 || num_threads>num_pts) {
                if(thread_id > 0) {
                    start_pt = 0;
                }
            }else{
                num_pts = trk.size() / (num_threads-1);
                start_pt = num_pts * thread_id;
                if(thread_id+1 == num_threads) num_pts = (trk.size() % (num_threads-1));
                //sleep(thread_id);
                //std::cout<<"Start " <<start_pt << " num pts " << num_pts << " / " << trk.size() << std::endl;
            }
            #endif

            auto const& vox_def = DetectorSpecs::GetME().GetVoxelDef();
            std::vector<double> local_pe_v(n_pmt,0.);
            std::vector<double> local_pe_refl_v(n_pmt,0.);
            double pos[3];
            size_t nproc0,nproc1; // for debug cout
            nproc0 = nproc1 = 0;  // for debug cout

            for(size_t ipt = start_pt; ipt < start_pt+num_pts; ++ipt) {

                auto const& pt = trk[ipt];
                pos[0] = pt.x;
                pos[1] = pt.y;
                pos[2] = pt.z;
                nproc0 += 1; // for debug cout
                //auto const& vis_pmt = lib_data[vox_id];
                //int vox_id = vox_def.GetVoxelID(pt.x,pt.y,pt.z);
                int vox_id = vox_def.GetVoxelID(pos);
                if (vox_id < 0) continue;

                nproc1 += 1; // for debug cout
                for(size_t ipmt=0; ipmt < n_pmt; ++ipmt) {

                    if(_channel_mask[ipmt]) continue;

                    if(_global_qe_refl > 0.)
                        local_pe_refl_v[ipmt] += pt.q * DetectorSpecs::GetME().GetVisibilityReflected(vox_id,ipmt);

                    if(!_uncoated_pmt_list[ipmt])
                        local_pe_v[ipmt] += pt.q * DetectorSpecs::GetME().GetVisibility(vox_id, ipmt);

                    //local_pe_v[ipmt] += pt.q * vis_pmt[ipmt];
                    //std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q0 << "," << q1 << "]" << std::endl;
                }
            }

            #if USING_LARSOFT == 0
            #pragma omp critical
            #endif
            {
                double qsum = 0.; // for debug cout
                for(size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
                    double q0 = (local_pe_v[ipmt] * _global_qe / _qe_v[ipmt]);
                    double q1 = (local_pe_refl_v[ipmt] * _global_qe_refl * _qe_v[ipmt]);
                    double q = q0 + q1;
                    flash.pe_v[ipmt] +=  q;
                    qsum += q; // for debug cout
                }

                // for debug cout
                //double qsum2 = 0.;
                //for(auto const& v: flash.pe_v) qsum2 += v;
                //std::cout<<"Thread ID " << thread_id << " ... " << start_pt << " => " << start_pt+num_pts << " ... " << nproc0 << "/" << nproc1 << " qsum " << qsum << " total sum " << qsum2 << std::endl<<std::flush;
            }

        }
        // for debug cout
        //double charge = 0.;
        //for(auto const& v : flash.pe_v) { charge+=v; }
        //if(charge>0)
        //    std::cout << std::endl << "Track size " << trk.size() << " points, total pe " << charge << std::endl << std::flush;
        //sleep(3);
        //throw OpT0FinderException();
        return;
    }

}
#endif
