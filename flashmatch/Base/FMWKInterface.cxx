#ifndef __OPT0FINDERFMWKINTERFACE_CXX__
#define __OPT0FINDERFMWKINTERFACE_CXX__

#include "FMWKInterface.h"
#include <assert.h>

namespace flashmatch{
  DetectorSpecs* DetectorSpecs::_me = nullptr;
}

#if USING_LARSOFT == 0
#include "flashmatch/Base/FMWKTools/PSetUtils.h"
#include "flashmatch/Base/FMWKTools/PhotonVisibilityService.h"
namespace flashmatch{
  
  DetectorSpecs::DetectorSpecs(std::string filename) {

    assert(!filename.empty());
    if(filename.find("/") != 0)
      filename = std::string(getenv("FMATCH_DATADIR")) + "/" + filename;

    auto cfg = CreatePSetFromFile(filename,"cfg");
    auto const& p = cfg.get<::flashmatch::Config_t>("DetectorSpecs");

    auto max_pt = p.get<std::vector<double> >("MaxPosition");
    auto min_pt = p.get<std::vector<double> >("MinPosition");
    assert(max_pt.size() == 3);
    assert(min_pt.size() == 3);
    assert(max_pt[0] >= min_pt[0] &&
	   max_pt[1] >= min_pt[1] &&
	   max_pt[2] >= min_pt[2]);
    _bbox = geoalgo::AABox(min_pt[0],min_pt[1],min_pt[2],max_pt[0],max_pt[1],max_pt[2]);
    std::cout<<_bbox.Min()[0]<<" "<<_bbox.Min()[1]<<" "<<_bbox.Min()[2]<<std::endl;
    std::cout<<_bbox.Max()[0]<<" "<<_bbox.Max()[1]<<" "<<_bbox.Max()[2]<<std::endl;
    size_t ch=0;
    _pmt_v.clear();
    while(1) {
      std::string key = "PMT" + std::to_string(ch);
      if(!p.contains_value(key)) break;
      geoalgo::Point_t pmt(p.get<std::vector<double> >(key));
      assert(pmt.size()==3);
      _pmt_v.push_back(pmt);
      ch++;
    }

    _drift_velocity = p.get<double>("DriftVelocity");
  }

  float DetectorSpecs::GetVisibility(double x, double y, double z, unsigned int opch) const
  { return phot::PhotonVisibilityService::GetME().GetVisibility(x,y,z,opch); }

#else
  DetectorSpecs::DetectorSpecs(std::string filename)
  {}
#endif
}

#endif

