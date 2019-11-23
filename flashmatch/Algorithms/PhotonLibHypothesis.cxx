#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX
#include <cassert>
#include "PhotonLibHypothesis.h"
#include "flashmatch/Base/OpT0FinderException.h"
#include "flashmatch/Base/FMWKInterface.h"
#include <chrono>

#include <omp.h>
#define NUM_THREADS 4

using namespace std::chrono;
namespace flashmatch {

  static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

  PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
  {}

  void PhotonLibHypothesis::_Configure_(const Config_t &pset)
  {
    omp_set_num_threads(NUM_THREADS);
    _global_qe = pset.get<double>("GlobalQE");
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

  void PhotonLibHypothesis::FillEstimate(const QCluster_t& old_trk, Flash_t &flash) const
  {
      QCluster_t trk = old_trk;
    if (_extend_tracks) {
        double min_x = kINVALID_DOUBLE; double max_x = -kINVALID_DOUBLE;
        size_t min_idx_x = 0; size_t max_idx_x = 0;
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
            for (size_t i = 0; i < num_pts+1; i++) {
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

        }
    }

    size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
    if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
    if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);
    if(flash.pe_true_v.empty()) flash.pe_true_v.resize(n_pmt);

    assert(flash.pe_v.size()      == n_pmt);
    assert(flash.pe_true_v.size() == n_pmt);
    assert(flash.pe_err_v.size()  == n_pmt);

    for (auto& v : flash.pe_v     ) v = 0;
    for (auto& v : flash.pe_err_v ) v = 0;
    for (auto& v : flash.pe_true_v ) v = 0;

    auto det = DetectorSpecs::GetME();

    auto const& lib_data = DetectorSpecs::GetME().GetPhotonLibraryData();

    //start = high_resolution_clock::now();
    #pragma omp parallel
    {
      size_t thread_id = omp_get_thread_num();
      size_t num_threads = omp_get_num_threads();
      size_t num_pts = trk.size() / num_threads;
      size_t start_pt = num_pts * thread_id;
      if(thread_id+1 == num_threads) num_pts += (trk.size() % num_threads);

      auto const& vox_def = DetectorSpecs::GetME().GetVoxelDef();
      // auto s = vox_def.GetVoxelSize();
      // auto s1 = vox_def.GetRegionLowerCorner();
      // auto s2 = vox_def.GetRegionUpperCorner();
      // std::cout << s[0] << " " << s[1] << " " << s[2] << std::endl;
      // std::cout << s1[0] << " " << s2[1] << " " << s1[2] << std::endl;
      // std::cout << s2[0] << " " << s2[1] << " " << s2[2] << std::endl;
      std::vector<double> local_pe_v(n_pmt,0);
      int vox_id;
      for( size_t ipt = start_pt; ipt < start_pt + num_pts; ++ipt) {
	auto const& pt = trk[ipt];
	vox_id = vox_def.GetVoxelID(pt.x,pt.y,pt.z);
	if (vox_id < 0) continue;
	auto const& vis_pmt = lib_data[vox_id];
	for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
	  local_pe_v[ipmt] += pt.q * vis_pmt[ipmt];
	}
      }
      #pragma omp critical
      for(size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {
	flash.pe_v[ipmt] += local_pe_v[ipmt] * _global_qe / _qe_v[ipmt];
      }

    }
    return;
  }
}
#endif
