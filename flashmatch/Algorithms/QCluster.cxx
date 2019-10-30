#ifndef QCluster_CXX
#define QCluster_CXX

#include "QCluster.h"
#include "flashmatch/GeoAlgo/GeoTrajectory.h"

namespace flashmatch {

  static QClusterFactory __global_QClusterFactory__;

  QCluster::QCluster(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
    , _gap         ( 0.5    )
    , _light_yield ( 40000. )
    , _dEdxMIP     ( 2.07   ) //1.42[Mev*cm^2*g]*1.4[g/cm^3]=2.004MeV/cm
  {}

  void QCluster::_Configure_(const Config_t &pset)
  {
    _gap          = pset.get< double > ( "SegmentSize" );
    _light_yield  = pset.get< double > ( "LightYield"  );
    _dEdxMIP      = pset.get< double > ( "MIPdEdx"     );
    _sigma_dedx   = pset.get< double > ( "SigmadEdx", 0);
    _trandom      = new TRandom();
  }

  double QCluster::GenerateLightYield(double dedx) const {
    // Generate smeared light yield using Poisson distribution
    double smeared_light_yield = _trandom->PoissonD(_light_yield);
    double smeared_dedx = dedx;
    if (_sigma_dedx > 0) {
      smeared_dedx = std::max(_trandom->Gaus(dedx, _sigma_dedx), 0.0);
    }
    // std::cout << smeared_light_yield << " " << _light_yield << std::endl;
    // std::cout << smeared_dedx << " " << dedx << std::endl;
    return smeared_dedx * smeared_light_yield;
  }

  void QCluster::MakeQCluster(const ::geoalgo::Vector& pt_1,
                           const ::geoalgo::Vector& pt_2,
                           QCluster_t& Q_cluster,
                           double dedx) const {

    if(dedx < 0) dedx = _dEdxMIP;

    double dist = pt_1.Dist(pt_2);
    QPoint_t q_pt;
    FLASH_INFO() << "Filling points between (" << pt_1[0] << "," << pt_1[1] << "," << pt_1[2] << ")"
		 << " => (" << pt_2[0] << "," << pt_2[1] << "," << pt_2[2] << ") ... dist="<<dist<<std::endl;
    if (dist <= _gap) {
      ::geoalgo::Vector mid_pt((pt_1 + pt_2) / 2.);
      q_pt.x = mid_pt[0];
      q_pt.y = mid_pt[1];
      q_pt.z = mid_pt[2];

      double smeared_factor = QCluster::GenerateLightYield(dedx);

      q_pt.q = smeared_factor * dist;
      FLASH_DEBUG() << "Smaller than gap threshold (" << _gap << ")" << std::endl
		   << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
      Q_cluster.emplace_back(q_pt);
      return;
    }

    int num_div = int(dist / _gap);

    ::geoalgo::Vector direct = (pt_1 - pt_2).Dir();

    Q_cluster.reserve(Q_cluster.size() + num_div);

    for (int div_index = 0; div_index < num_div + 1; div_index++) {
      if (div_index < num_div) {
        auto const mid_pt = pt_2 + direct * (_gap * div_index + _gap / 2.);
        q_pt.x = mid_pt[0] ;
        q_pt.y = mid_pt[1];
        q_pt.z = mid_pt[2];
        double smeared_factor = QCluster::GenerateLightYield(dedx);
        q_pt.q = _gap * smeared_factor;
	FLASH_DEBUG() << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
        Q_cluster.emplace_back(q_pt);
      }
      else {
        double weight = (dist - int(dist / _gap) * _gap);
        auto const mid_pt = pt_2 + direct * (_gap * div_index + weight / 2.);
        q_pt.x = mid_pt[0] ;
        q_pt.y = mid_pt[1];
        q_pt.z = mid_pt[2];
        double smeared_factor = QCluster::GenerateLightYield(dedx);
        q_pt.q = weight * smeared_factor;
	FLASH_DEBUG() << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
        Q_cluster.emplace_back(q_pt);
      }//Last segment less than gap
    }
  }

  QCluster_t QCluster::FlashHypothesis(const ::geoalgo::Trajectory& trj) const {

    QCluster_t result;
    result.clear();

    for (size_t i = 0; i < trj.size() - 1; i++) {
      auto const& this_loc(trj[i]);
      auto const& last_loc(trj[i + 1]);
      QCluster::MakeQCluster(this_loc, last_loc, result);
    }

    FLASH_INFO() << result << std::endl;

    return result;
  }

}


#endif
