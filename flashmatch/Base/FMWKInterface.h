#ifndef __OPT0FINDERFMWKINTERFACE_H__
#define __OPT0FINDERFMWKINTERFACE_H__

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "FMWKTools/ConfigManager.h"
#include "FMWKTools/PhotonVoxels.h"
namespace flashmatch {
  /// Configuration object
  using Config_t = flashmatch::PSet;
}
#else
namespace flashmatch{
  /// Configuration object
  using Config_t = fhicl::FhiclParameterSet;
}
#endif

#include "flashmatch/GeoAlgo/GeoAABox.h"
namespace flashmatch {

  class DetectorSpecs {
    
  public:
    DetectorSpecs(std::string filename="specs.cfg");
    ~DetectorSpecs(){}
    
    inline static DetectorSpecs& GetME(std::string filename="detector_specs.cfg")
    {
      if(!_me) _me = new DetectorSpecs(filename);
      return *_me;
    }

    /// PMT XYZ position filler
    inline const geoalgo::Point_t& PMTPosition(size_t opch) { return _pmt_v.at(opch); }

    /// Detector active volume
    inline const geoalgo::AABox& ActiveVolume() const { return _bbox; }
    
    /// # of PMTs
    inline size_t NOpDets() const { return _pmt_v.size(); }

    /// Drift velocity
    inline double DriftVelocity() const { return _drift_velocity; }

    /// Visibility
    float GetVisibility(double x, double y, double z, unsigned int opch) const;

    /// Photon Library data access FIXME
    const std::vector<std::vector<float > >& GetPhotonLibraryData() const;

    /// Voxel definition
    inline const sim::PhotonVoxelDef& GetVoxelDef() const { return _voxel_def; }

  private:
    static DetectorSpecs* _me;
    std::vector<geoalgo::Point_t> _pmt_v;
    geoalgo::AABox _bbox;
    double _drift_velocity;
    sim::PhotonVoxelDef _voxel_def;
  };

}
#endif
