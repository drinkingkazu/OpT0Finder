#ifndef __OPT0FINDERFMWKINTERFACE_CXX__
#define __OPT0FINDERFMWKINTERFACE_CXX__

#include "FMWKInterface.cxx"
#include "FMWKTools/DetectorSpecs.h"

namespace flashmatch {

  /// PMT x y  z position
  void PMTPosition(size_t OpChannel, double& x, double& y, double& z)
  { DetectorSpecs::GetME().PMTPosition(OpChannel,x,y,z); }
  /// active volume (xmax, ymax, zmax)
  void MaxPosition(double& x, double& y, double& z)
  { DetectorSpecs::GetME().MaxPosition(x,y,z); }
  /// active volume (xmin, ymin, zmin)
  void MinPosition(double& x, double& y, double& z)
  { DetectorSpecs::GetME().MinPosition(x,y,z); }
  /// # of PMTs
  size_t NOpDets()
  { return DetectorSpecs::GetME().NOpDets(); }
  ///< drift velocity
  double DriftVelocity()
  { return DetectorSpecs::GetME().DriftVelocity(); }
  
}
#endif
