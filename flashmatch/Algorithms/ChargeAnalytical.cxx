#ifndef CHARGEANALYTICAL_CXX
#define CHARGEANALYTICAL_CXX

#include "ChargeAnalytical.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#include <math.h>

namespace flashmatch {

  static ChargeAnalyticalFactory __global_ChargeAnalyticalFactory__;
  
  ChargeAnalytical::ChargeAnalytical(const std::string name)
    : BaseFlashHypothesis(name)
  {}

  void ChargeAnalytical::_Configure_(const Config_t &pset)
  {
    _global_qe = pset.get<double>("GlobalQE");
    _qe_v      = pset.get<std::vector<double> >("CCVCorrection");
    if(_qe_v.size() != flashmatch::NOpDets()) {
      FLASH_CRITICAL() << "CCVCorrection factor array has size " << _qe_v.size()
		       << " != number of opdet (" << flashmatch::NOpDets() << ")!" << std::endl;
      throw OpT0FinderException();
    }
  }

  void ChargeAnalytical::FillEstimate(const QCluster_t &track,
				      Flash_t &flash) const
  {
    
    size_t n_pmt = flashmatch::NOpDets();
    
    for (size_t i = 0; i < n_pmt; ++i) {
      flash.pe_v[i] = 0;
    }

    double x,y,z;
    
    for (size_t pmt_index = 0; pmt_index < n_pmt; ++pmt_index) {
      
      for (size_t pt_index = 0; pt_index < track.size(); ++pt_index) {
	
	auto const &pt = track[pt_index];

	flashmatch::PMTPosition(pmt_index,x,y,z);
	double dx = x - pt.x;
	double dy = y - pt.y;
	double dz = z - pt.z;
	
	double r2 = (pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
	
	double angle = dx / sqrt(r2);
	
	if (angle < 0) angle *= -1;
	
	flash.pe_v[pmt_index] += pt.q * angle / r2 * _global_qe / _qe_v[pmt_index];
	
      }
    }
  }
}
#endif

