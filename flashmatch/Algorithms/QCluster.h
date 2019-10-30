/**
 * \file QCluster.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class QCluster
 *
 * @author Rui
 */

/** \addtogroup Algorithms

    @{*/
#ifndef QCluster_H
#define QCluster_H

#include <iostream>
#include <numeric>
#include <functional>
#include <algorithm>

#include "flashmatch/GeoAlgo/GeoTrajectory.h"
#include "flashmatch/Base/BaseAlgorithm.h"
#include "flashmatch/Base/CustomAlgoFactory.h"
#include <TRandom.h>

namespace flashmatch{
/**
   \class QCluster
   User defined class QCluster ... these comments are used to generate
   doxygen documentation!
 */

  class QCluster : public flashmatch::BaseAlgorithm {

  public:

    /// Default constructor
    QCluster(const std::string name="QCluster");

    /// Default destructor
    ~QCluster(){}

    // Setter function
    double Set_Gap      ( double x) { _gap   =x;      return _gap;}

    // Flash Hypothesis for Trajectory (Track)
    flashmatch::QCluster_t FlashHypothesis(const ::geoalgo::Trajectory& trj) const;

    void MakeQCluster(const ::geoalgo::Vector& pt_1,
                  const ::geoalgo::Vector& pt_2,
                  flashmatch::QCluster_t& Q_cluster,
		  double dedx=-1) const;

    // Getter for light yield configured paramater
    double GetLightYield() const { return _light_yield; }
    double GenerateLightYield(double dedx) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _gap;
    double _light_yield;
    double _dEdxMIP;
    double _sigma_dedx;
    TRandom* _trandom;
  };

  /**
     \class flashmatch::QClusterFactory
  */
  class QClusterFactory : public CustomAlgoFactoryBase {
  public:
    /// ctor
    QClusterFactory() { CustomAlgoFactory::get().add_factory("QCluster",this); }
    /// dtor
    ~QClusterFactory() {}
    /// creation method
    BaseAlgorithm* create(const std::string instance_name) { return new QCluster(instance_name); }
  };
}

#endif
/** @} */ // end of doxygen group
