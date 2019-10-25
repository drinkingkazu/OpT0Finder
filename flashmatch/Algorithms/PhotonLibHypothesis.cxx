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
    _qe_v.clear();
    _qe_v = pset.get<std::vector<double> >("CCVCorrection",_qe_v);
    if(_qe_v.empty()) _qe_v.resize(DetectorSpecs::GetME().NOpDets(),1.0);
    if(_qe_v.size() != DetectorSpecs::GetME().NOpDets()) {
      FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
		       << " != number of opdet (" << DetectorSpecs::GetME().NOpDets() << ")!" << std::endl;
      throw OpT0FinderException();
    }
  }

  void PhotonLibHypothesis::FillEstimate(const QCluster_t& trk, Flash_t &flash) const
  {
    size_t n_pmt = DetectorSpecs::GetME().NOpDets();//n_pmt returns 0 now, needs to be fixed
    if(flash.pe_v.empty()) flash.pe_v.resize(n_pmt);
    if(flash.pe_err_v.empty()) flash.pe_err_v.resize(n_pmt);

    assert(flash.pe_v.size()     == n_pmt);
    assert(flash.pe_err_v.size() == n_pmt);

    for (auto& v : flash.pe_v     ) v = 0;
    for (auto& v : flash.pe_err_v ) v = 0;

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
      
      std::vector<double> local_pe_v(n_pmt,0);
      int vox_id;
      for( size_t ipt = start_pt; ipt < start_pt + num_pts; ++ipt) {
	auto const& pt = trk[ipt];
	vox_id = vox_def.GetVoxelID(pt.x,pt.y,pt.z);
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
