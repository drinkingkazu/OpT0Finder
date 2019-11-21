/**
 * \file PhotonLibHypothesis.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class PhotonLibHypothesis
 *
 * @author yuntse
 */

/** \addtogroup Algorithms

    @{*/

#ifndef PHOTONLIBHYPOTHESIS_H
#define PHOTONLIBHYPOTHESIS_H

#include <iostream>
#include "flashmatch/Base/BaseFlashHypothesis.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"

namespace flashmatch {
  /**
     \class PhotonLibHypothesis
     User custom analysis class made by SHELL_USER_NAME
   */
  class PhotonLibHypothesis : public BaseFlashHypothesis {

  public:

    /// Default constructor
    PhotonLibHypothesis(const std::string name="PhotonLibHypothesis");

    /// Default destructor
    virtual ~PhotonLibHypothesis(){}

    void FillEstimate(const QCluster_t&, Flash_t&) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _global_qe;             ///< Global QE
    double _sigma_qe;              ///< Sigma for Gaussian centered on Global QE
    std::vector<double> _qe_v;     ///< PMT-wise relative QE
    double _segment_size;
    bool _extend_tracks;
    double _threshold_extend_track;
  };

  /**
     \class flashmatch::PhotonLibHypothesisFactory
  */
  class PhotonLibHypothesisFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    PhotonLibHypothesisFactory() { FlashHypothesisFactory::get().add_factory("PhotonLibHypothesis",this); }
    /// dtor
    ~PhotonLibHypothesisFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new PhotonLibHypothesis(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group
