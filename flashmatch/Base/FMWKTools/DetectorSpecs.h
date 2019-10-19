#ifndef __DETECTORSPECS_H__
#define __DETECTORSPECS_H__

#include "ConfigManager.h"

namespace flashmatch {

  class DetectorSpecs {
    
  public:
    DetectorSpecs(std::string filename="specs.cfg");
    ~DetectorSpecs(){}
    
    static DetectorSpecs& GetME(std::string filename="detector_specs.cfg")
    {
      if(!_me) _me = new DetectorSpecs(filename);
      return *_me;
    }

    /// PMT X positions
    inline const std::vector<double>& PMTPositionX() const { return _xpos_v; }

    /// PMT Y positions
    inline const std::vector<double>& PMTPositionY() const { return _ypos_v; }

    /// PMT Z positions
    inline const std::vector<double>& PMTPositionZ() const { return _zpos_v; }

    /// PMT XYZ position filler
    void PMTPosition(size_t OpChannel, double& x, double& y, double& z) const;

    /// Detector max X
    inline double MaxX() const { return _xmax; }
    /// Detector max Y
    inline double MaxY() const { return _ymax; }
    /// Detector max Z
    inline double MaxZ() const { return _zmax;}

    /// Detector min X
    inline double MinX() const { return _xmin; }
    /// Detector min Y
    inline double MinY() const { return _ymin; }
    /// Detector min Z
    inline double MinZ() const { return _zmin; }

    /// Detector max position filler
    inline void MaxPosition(double& x, double& y, double& z) const { x=_xmax; y=_ymax; z=_zmax; }
    /// Detector min position filler
    inline void MinPosition(double& x, double& y, double& z) const { x=_xmin; y=_ymin; z=_zmin; }
    
    /// # of PMTs
    inline size_t NOpDets() const { return _xpos_v.size(); }

    /// Drift velocity
    inline double DriftVelocity() const { return _drift_velocity; }

  private:
    static DetectorSpecs* _me;
    size_t _num_pmts;
    std::vector<double> _xpos_v;
    std::vector<double> _ypos_v;
    std::vector<double> _zpos_v;
    double _xmax, _ymax, _zmax;
    double _xmin, _ymin, _zmin;
    double _drift_velocity;
  };

}
#endif
