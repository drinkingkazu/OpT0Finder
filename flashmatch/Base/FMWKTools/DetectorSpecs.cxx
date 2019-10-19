#ifndef __DETECTOR_SPECS_CXX__
#define __DETECTOR_SPECS_CXX__
#include "DetectorSpecs.h"
#include "PSetUtils.h"
#include <assert.h>
namespace flashmatch{

  DetectorSpecs* DetectorSpecs::_me = nullptr;

  DetectorSpecs::DetectorSpecs(std::string filename) {

    assert(!filename.empty());
    if(filename.find("/") != 0)
      filename = std::string(getenv("FMATCH_DATADIR")) + "/" + filename;

    auto p = CreatePSetFromFile(filename,"cfg");

    auto max_pt = p.get<std::vector<double> >("ActiveMax");
    auto min_pt = p.get<std::vector<double> >("ActiveMin");

    assert(max_pt.size() == 3);
    _xmax = max_pt[0];
    _ymax = max_pt[1];
    _zmax = max_pt[2];
    
    assert(min_pt.size() == 3);
    _xmin = min_pt[0];
    _ymin = min_pt[1];
    _zmin = min_pt[2];

    size_t ch=0;
    while(1) {
      std::string key = "PMT" + std::to_string(ch);
      if(!p.contains_value(key)) break;
      auto pt = p.get<std::vector<double> >(key);
      assert(pt.size()==3);
      _xpos_v.push_back(pt[0]);
      _ypos_v.push_back(pt[1]);
      _zpos_v.push_back(pt[2]);
    }

    _drift_velocity = p.get<double>("DriftVelocity");
  }

  void DetectorSpecs::PMTPosition(size_t ch, double& x, double& y, double& z) const
  {
    assert(ch < _xpos_v.size());
    x = _xpos_v[ch];
    y = _ypos_v[ch];
    z = _zpos_v[ch];
  }
  

}

#endif
