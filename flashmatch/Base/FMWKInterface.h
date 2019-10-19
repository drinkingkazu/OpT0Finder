#ifndef __OPT0FINDERFMWKINTERFACE_H__
#define __OPT0FINDERFMWKINTERFACE_H__

#include "FMWKTools/ConfigManager.h"
#include "FMWKTools/PhotonVisibilityService.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = flashmatch::PSet;
  /// PMT x y  z position
  void PMTPosition(size_t OpChannel, double& x, double& y, double& z);
  /// active volume (xmax, ymax, zmax)
  void MaxPosition(double& x, double& y, double& z);
  /// active volume (xmin, ymin, zmin)
  void MinPosition(double& x, double& y, double& z);
  /// # of PMTs
  size_t NOpDets();
  ///< drift velocity
  double DriftVelocity();

}
#endif
