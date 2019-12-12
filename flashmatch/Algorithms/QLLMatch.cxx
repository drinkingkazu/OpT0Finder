#ifndef QLLMATCH_CXX
#define QLLMATCH_CXX

#include "QLLMatch.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#include <TMinuit.h>
#include <cmath>
#include <numeric>
#include <TMath.h>
#include <cassert>
#include <chrono>
#include <climits>
using namespace std::chrono;
namespace flashmatch {

  static QLLMatchFactory __global_QLLMatchFactory__;

  QLLMatch *QLLMatch::_me = nullptr;

  void MIN_vtx_qll(Int_t &, Double_t *, Double_t &, Double_t *, Int_t);

  QLLMatch::QLLMatch(const std::string name)
    : BaseFlashMatch(name), _mode(kChi2), _record(false), _normalize(false), _poisson("QLLMatch_poisson","TMath::Poisson(x,y)"), _minuit_ptr(nullptr)
  { _current_llhd = _current_chi2 = -1.0; }

  QLLMatch::QLLMatch()
  { throw OpT0FinderException("Use QLLMatch::GetME() to obtain singleton pointer!"); }

  void QLLMatch::_Configure_(const Config_t &pset) {
    _record = pset.get<bool>("RecordHistory");
    _check_touching_track = pset.get<bool>("CheckTouchingTrack");
    _normalize = pset.get<bool>("NormalizeHypothesis");
    _mode   = (QLLMode_t)(pset.get<unsigned short>("QLLMode"));
    _pe_observation_threshold = pset.get<double>("PEObservationThreshold", 1.e-6);
    _pe_hypothesis_threshold  = pset.get<double>("PEHypothesisThreshold", 1.e-6);
    _migrad_tolerance         = pset.get<double>("MIGRADTolerance", 0.1);
    _offset                   = pset.get<double>("Offset", 0.0);
    _touching_track_window    = pset.get<double>("TouchingTrackWindow", 5.0);
    // _custom_algo              = pset.get<std::string>("CustomAlgo", "");
    // if (!_custom_algo.empty()) _alg_custom_algo = CustomAlgoFactory::get().create(_custom_algo, _custom_algo);
    // if (_alg_custom_algo) {
    //
    // }
    _penalty_threshold_v = pset.get<std::vector<double> >("PEPenaltyThreshold");
    _penalty_value_v = pset.get<std::vector<double> >("PEPenaltyValue");

    _recox_penalty_threshold = pset.get<double>("XPenaltyThreshold");
    _recoz_penalty_threshold = pset.get<double>("ZPenaltyThreshold");

    _onepmt_score_threshold = pset.get<double>("OnePMTScoreThreshold");
    _onepmt_xdiff_threshold = pset.get<double>("OnePMTXDiffThreshold");
    _onepmt_pesum_threshold = pset.get<double>("OnePMTPESumThreshold");
    _onepmt_pefrac_threshold = pset.get<double>("OnePMTPEFracThreshold");

    _xpos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    _ypos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    _zpos_v.resize(DetectorSpecs::GetME().NOpDets(),0.);
    for(size_t ch=0; ch<DetectorSpecs::GetME().NOpDets(); ++ch) {
      auto const& pmt_pos = DetectorSpecs::GetME().PMTPosition(ch);
      _xpos_v[ch] = pmt_pos[0];
      _ypos_v[ch] = pmt_pos[1];
      _zpos_v[ch] = pmt_pos[2];
    }
    auto const& bbox = DetectorSpecs::GetME().ActiveVolume();
    _vol_xmax = bbox.Max()[0];
    _vol_xmin = bbox.Min()[0];
  }

  FlashMatch_t QLLMatch::Match(const QCluster_t &pt_v, const Flash_t &flash, const bool prohibit_touch_match) {

    //
    // Prepare TPC
    //
    _raw_trk.resize(pt_v.size());
    double min_x =  1e20;
    double max_x = -1e20;
    for (size_t i = 0; i < pt_v.size(); ++i) {
      auto const &pt = pt_v[i];
      _raw_trk[i] = pt;
      if (pt.x < min_x) { min_x = pt.x; _raw_xmin_pt = pt; }
      if (pt.x > max_x) { max_x = pt.x; _raw_xmax_pt = pt; }
    }
    for (auto &pt : _raw_trk) pt.x -= min_x;

    auto res1 = PESpectrumMatch(pt_v,flash,true, prohibit_touch_match);
    auto res2 = PESpectrumMatch(pt_v,flash,false, prohibit_touch_match);
    FLASH_INFO() << "Using   mid-x-init ... maximized 1/param Score=" << res1.score << " @ X=" << res1.tpc_point.x << " [cm]" << std::endl;
    FLASH_INFO() << "Without mid-x-init ... maximized 1/param Score=" << res2.score << " @ X=" << res2.tpc_point.x << " [cm]" << std::endl;

    auto res = (res1.score > res2.score ? res1 : res2);
    //auto res = res2;
    /*
    if(res.score < _onepmt_score_threshold) {

      FLASH_INFO() << "Resulting score below OnePMTScoreThreshold... calling OnePMTMatch" << std::endl;
      auto res_onepmt = OnePMTMatch(flash);

      if(res_onepmt.score >= 0.)
	return res_onepmt;
    }
    */
    return res;
  }

  FlashMatch_t QLLMatch::OnePMTMatch(const Flash_t& flash) {

    FlashMatch_t res;
    res.score=-1;
    res.num_steps = 1;
    res.minimizer_min_x = 0;
    res.minimizer_max_x = 0;
    // Check if pesum threshold condition is met to use this method
    double pesum = flash.TotalPE();
    if(pesum < _onepmt_pesum_threshold) {
      //std::cout <<"PESumThreshold not met (pesum=" << pesum << ")" << std::endl;
      return res;
    }

    // Check if pe max fraction condition is met to use this method
    size_t maxpmt = 0;
    double maxpe  = 0.;
    for(size_t pmt=0; pmt<flash.pe_v.size(); ++pmt) {
      if(flash.pe_v[pmt] < maxpe) continue;
      maxpe  = flash.pe_v[pmt];
      maxpmt = pmt;
    }
    if(maxpe / pesum < _onepmt_pefrac_threshold) {
      //std::cout << "PERatioThreshold not met (peratio=" << maxpe/pesum << ")" << std::endl;
      return res;
    }

    // Now see if Flash T0 can be consistent with an assumption MinX @ X=0.
    double xdiff = fabs(_raw_xmin_pt.x - flash.time * DetectorSpecs::GetME().DriftVelocity());
    if( xdiff > _onepmt_xdiff_threshold ) {
      //std::cout << "XDiffThreshold not met (xdiff=" << xdiff << ")" << std::endl;
      return res;
    }

    // Reaching this point means it is an acceptable match
    _reco_x_offset = 0.;
    _reco_x_offset_err = std::fabs(_xpos_v.at(maxpmt));

    // Compute hypothesis with MinX @ X=0 assumption.
    _hypothesis = flash;
    FillEstimate(_raw_trk,_hypothesis);
    res.hypothesis = _hypothesis.pe_v;

    // Compute TPC point
    res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;
    double weight = 0;
    for (size_t pmt_index = 0; pmt_index < DetectorSpecs::GetME().NOpDets(); ++pmt_index) {

      res.tpc_point.y += _ypos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];
      res.tpc_point.z += _zpos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];

      weight += _hypothesis.pe_v[pmt_index];
    }

    res.tpc_point.y /= weight;
    res.tpc_point.z /= weight;

    res.tpc_point.x = _reco_x_offset;
    res.tpc_point_err.x = _reco_x_offset_err;

    // Use MinX point YZ distance to the max PMT and X0 diff as weight
    res.score = 1.;
    // FIXME For now we do not have time distribution, so ignore this weighting
    //res.score *= 1. / xdiff;
    res.score *= 1. / (sqrt(pow(_raw_xmin_pt.y - _ypos_v.at(maxpmt),2) + pow(_raw_xmin_pt.z - _zpos_v.at(maxpmt),2)));

    return res;

  }

  FlashMatch_t QLLMatch::TouchingTrack(const QCluster_t &pt_v, const Flash_t & flash, double score, bool tpc0) {
      // Claim a match and don't call minuit.
      FlashMatch_t res;
      res.score = score;
      res.num_steps = 0;
      res.minimizer_min_x =  pt_v.min_x();//tpc0 ? pt_v.min_x() : pt_v.max_x();
      res.minimizer_max_x = pt_v.min_x();//tpc0 ? pt_v.min_x() : pt_v.max_x();
      res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = kINVALID_DOUBLE;
      // Find point with min x
      for(auto const& pt : pt_v) {
          if (pt.x == res.minimizer_min_x) {
              res.tpc_point.x = pt.x;
              res.tpc_point.y = pt.y;
              res.tpc_point.z = pt.z;
          }
      }
      //std::cout << flash.time << " " << pt_v.min_x() << " " << flash.time_true << " " << pt_v.time_true << " " << res.tpc_point.x << std::endl;
      res.tpc_point.x = res.tpc_point.x + (tpc0 ? -1.0 : 1.0) * (flash.time) * DetectorSpecs::GetME().DriftVelocity();
      //std::cout << res.tpc_point.x << std::endl;
      //res.hypothesis = flash.pe_v;
      Flash_t hypothesis;
      FillEstimate(pt_v, hypothesis);
      res.hypothesis = hypothesis.pe_v;
      res.touch_match = (score == 10.) ? 1 : 2;
      return res;
  }

  FlashMatch_t QLLMatch::PESpectrumMatch(const QCluster_t &pt_v, const Flash_t &flash, const bool init_x0, const bool prohibit_touch_match) {

    // before calling minuit, check if this is a touching track (score 10)
    // or a track crossing through both tpc planes (score 20)
    bool touching_both_planes = (pt_v.max_x() - pt_v.min_x() >= DetectorSpecs::GetME().ActiveVolume().Max()[0] - DetectorSpecs::GetME().ActiveVolume().Min()[0]);
    double time_tpc0 = (pt_v.min_x() - DetectorSpecs::GetME().ActiveVolume().Min()[0]) / DetectorSpecs::GetME().DriftVelocity();
    double time_tpc1 = (DetectorSpecs::GetME().ActiveVolume().Max()[0] - pt_v.max_x()) / DetectorSpecs::GetME().DriftVelocity();
    bool touching_tpc0 = std::fabs(flash.time - time_tpc0) < _touching_track_window;
    bool touching_tpc1 = std::fabs(flash.time - time_tpc1) < _touching_track_window;
    if (!prohibit_touch_match) {
        if (_check_touching_track || touching_both_planes) {
            double score = touching_both_planes ? 20.0 : 10.0;
            // if (std::fabs(pt_v.time_true - flash.time_true) < 1) {
            //     std::cout << flash.time << " " << time_tpc0 << " " << time_tpc1 << " " << pt_v.min_x() << " " << pt_v.max_x() << " " << flash.time_true << " " << pt_v.time_true << std::endl;
            // }
            if (touching_tpc0 || touching_tpc1) {
                return this->TouchingTrack(pt_v, flash, score, touching_tpc0);
            }
            else if (touching_both_planes) {
                // If the track is crossing both tpc planes, regardless of _check_touching_track
                // we do not want to go through minuit (might return nan).
                // std::cout << " touching both planes " << std::endl;
                FlashMatch_t res;
                res.num_steps = _num_steps;
                res.minimizer_min_x = _minimizer_min_x;
                res.minimizer_max_x = _minimizer_max_x;
                res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;
                res.score = 0;
                res.touch_match = 2;
                return res;
            }
        }
    }

    FlashMatch_t res;
    res.num_steps = _num_steps;
    res.minimizer_min_x = _minimizer_min_x;
    res.minimizer_max_x = _minimizer_max_x;
    res.tpc_point.x = res.tpc_point.y = res.tpc_point.z = 0;
    res.score = 0;
    if (touching_tpc0 || touching_tpc1) {
        res.touch_match = 1;
    }
    else if (touching_both_planes) {
        res.touch_match = 2;
    }
    else {
        res.touch_match = 0;
    }

    FLASH_INFO() << "true time " << pt_v.time_true << " " << flash.time_true << std::endl;
    auto return_value = this->CallMinuit(pt_v, flash, init_x0);
    // Shit happens line above in CallMinuit

    // Estimate position
    if (std::isnan(_qll)) return res;

    double weight = 0;

    for (size_t pmt_index = 0; pmt_index < DetectorSpecs::GetME().NOpDets(); ++pmt_index) {

      res.tpc_point.y += _ypos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];
      res.tpc_point.z += _zpos_v.at(pmt_index) * _hypothesis.pe_v[pmt_index];

      weight += _hypothesis.pe_v[pmt_index];
    }

    res.tpc_point.y /= weight;
    res.tpc_point.z /= weight;

    res.tpc_point.x = _reco_x_offset;
    res.tpc_point_err.x = _reco_x_offset_err;

    res.hypothesis  = _hypothesis.pe_v;

    //
    // Compute score
    //
    if(_mode == kSimpleLLHD)
      res.score = _qll * -1.;
    else
      res.score = 1. / _qll;

    // Compute X-weighting
    /*
    double x0 = _raw_xmin_pt.x - flash.time * DetectorSpecs::GetME().DriftVelocity();
    if( fabs(_reco_x_offset - x0) > _recox_penalty_threshold )
      res.score *= 1. / (1. + fabs(_reco_x_offset - x0) - _recox_penalty_threshold);
    // Compute Z-weighting
    double z0 = 0;
    weight = 0;
    for (size_t pmt_index = 0; pmt_index < DetectorSpecs::GetME().NOpDets(); ++pmt_index) {
      z0 += _zpos_v.at(pmt_index) * flash.pe_v[pmt_index];
      weight += flash.pe_v[pmt_index];
    }
    z0 /= weight;
    if( fabs(res.tpc_point.z - z0) > _recoz_penalty_threshold )
      res.score *= 1. / (1. + fabs(res.tpc_point.z - z0) - _recoz_penalty_threshold);
    */

    // FIXME
    // If this was re-running flash matching on a touching track, we still want to keep
    // the track min-x as best estimate instead of the minimizer estimate?
    if (_check_touching_track && prohibit_touch_match && res.touch_match == 1) {
        FlashMatch_t res2 = this->TouchingTrack(pt_v, flash, res.score, touching_tpc0);
        res.hypothesis = res2.hypothesis;
        res.tpc_point = res2.tpc_point;
    }
    return res;
  }

  const Flash_t &QLLMatch::ChargeHypothesis(const double xoffset) {
    //auto start = high_resolution_clock::now();
    if (_hypothesis.pe_v.empty()) _hypothesis.pe_v.resize(DetectorSpecs::GetME().NOpDets(), 0.);
    if (_hypothesis.pe_v.size() != DetectorSpecs::GetME().NOpDets()) {
      throw OpT0FinderException("Hypothesis vector length != PMT count");
    }

    for (auto &v : _hypothesis.pe_v) v = 0;

    // Apply xoffset
    _var_trk.resize(_raw_trk.size());
    for (size_t pt_index = 0; pt_index < _raw_trk.size(); ++pt_index) {
      //std::cout << "x point : " << _raw_trk[pt_index].x << "\t offset : " << xoffset << std::endl;
      _var_trk[pt_index].x = _raw_trk[pt_index].x + xoffset;
      _var_trk[pt_index].y = _raw_trk[pt_index].y;
      _var_trk[pt_index].z = _raw_trk[pt_index].z;
      _var_trk[pt_index].q = _raw_trk[pt_index].q;
    }
    //auto end = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 1 = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    FillEstimate(_var_trk, _hypothesis);
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 2 = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    if (_normalize) {
      double qsum = std::accumulate(std::begin(_hypothesis.pe_v),
				    std::end(_hypothesis.pe_v),
				    0.0);
      for (auto &v : _hypothesis.pe_v) v /= qsum;
    }
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis 3 = " << duration.count() << "us" << std::endl;

    return _hypothesis;
  }

  const Flash_t &QLLMatch::Measurement() const { return _measurement; }

  double QLLMatch::QLL(const Flash_t &hypothesis,
		       const Flash_t &measurement) {

    double nvalid_pmt = 0;

    double PEtot_Hyp = 0;
    for (auto const &pe : hypothesis.pe_v)
      PEtot_Hyp += pe;
    double PEtot_Obs = 0;
    for (auto const &pe : measurement.pe_v)
      PEtot_Obs += pe;

    _current_chi2 = _current_llhd = 0.;

    if (measurement.pe_v.size() != hypothesis.pe_v.size())
      throw OpT0FinderException("Cannot compute QLL for unmatched length!");

    double O, H, Error;
    const double epsilon = 1.e-6;
    for (size_t pmt_index = 0; pmt_index < hypothesis.pe_v.size(); ++pmt_index) {

      O = measurement.pe_v[pmt_index]; // observation
      H = hypothesis.pe_v[pmt_index];  // hypothesis

      if( H < 0 ) throw OpT0FinderException("Cannot have hypothesis value < 0!");

      if(O < _pe_observation_threshold) {
        if (!_penalty_value_v.empty()) {
          O = _penalty_value_v[pmt_index];
        }
        else {
          O = _pe_observation_threshold;
        }
      }
      if (H < _pe_hypothesis_threshold) {
        if(!_penalty_threshold_v.empty()) {
          H = _penalty_threshold_v[pmt_index];
        }
        else {
          H = _pe_hypothesis_threshold;
        }
      }

      if(_mode == kLLHD) {
	assert(H>0);
	double arg = TMath::Poisson(O,H) + epsilon;
	if(!std::isnan(arg) && !std::isinf(arg)) {
	  _current_llhd -= std::log10(arg);
	  nvalid_pmt += 1;
	  //if(_converged) FLASH_INFO() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg) << std::endl;
	}
      }
      else if(_mode == kWeightedLLHD) {
	assert(H>0);
	double arg = TMath::Poisson(O,H) + epsilon;
	if(!std::isnan(arg) && !std::isinf(arg)) {
	  _current_llhd -= std::log10(arg * sqrt(std::max(H,epsilon)));
	  nvalid_pmt += 1;
	  if(_converged) FLASH_INFO() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg * sqrt(std::max(H,epsilon))) << std::endl;
	}
      }else if(_mode == kIntegralLLHD) {
	double hmin = H-0.5;
	double hmax = H+0.5;
	double omin = O-sqrt(O);
	double omax = O+sqrt(O);
	if(hmin<0.) hmin = 0.;
	if(omin<0.) omin = 0.;
	if(hmax<hmin+1) hmax = hmin+1;
	if(omax<omin+1) omax = omin+1;
	double arg = _poisson.Integral(omin,omax,hmin,hmax) + epsilon;
	if(!std::isnan(arg) && !std::isinf(arg)) {
	  _current_llhd -= std::log10(arg);
	  nvalid_pmt += 1;
	  if(_converged) FLASH_INFO() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " LHD "<<arg << " -LLHD " << -1 * std::log10(arg) << std::endl;
	}
      }else if (_mode == kSimpleLLHD) {

	double arg = (H - O * std::log(H));
	_current_llhd += arg;
	if(_converged) FLASH_INFO() <<"PMT "<<pmt_index<<" O/H " << O << " / " << H << " ... -LLHD " << arg << std::endl;
	//nvalid_pmt += 1;

      } else if (_mode == kChi2) {

	Error = O;
	if( Error < 1.0 ) Error = 1.0;
	_current_chi2 += std::pow((O - H), 2) / (Error);
	nvalid_pmt += 1;

      } else {
	FLASH_ERROR() << "Unexpected mode" << std::endl;
	throw OpT0FinderException();
      }

    }
    //FLASH_DEBUG() <<"Mode " << (int)(_mode) << " Chi2 " << _current_chi2 << " LLHD " << _current_llhd << " nvalid " << nvalid_pmt << std::endl;

    _current_chi2 /= nvalid_pmt;
    _current_llhd /= (nvalid_pmt +1);
    if(_converged)
      FLASH_INFO() << "Combined LLHD: " << _current_llhd << " (divided by nvalid_pmt+1 = " << nvalid_pmt+1<<")"<<std::endl;

    return (_mode == kChi2 ? _current_chi2 : _current_llhd);
  }

  void MIN_vtx_qll(Int_t & /*Npar*/, // Number of parameters
		   Double_t * /*Grad*/, // Partial derivatives (return values)
		   Double_t &Fval, // Function value (return value)
		   Double_t *Xval, // Parameter values
		   Int_t) /*Flag*/{ // flag word
    //std::cout << "minuit offset : " << Fval << std::endl;
    //std::cout << "minuit Xval?? : " << *Xval << std::endl;

    //auto start = high_resolution_clock::now();
    auto const &hypothesis = QLLMatch::GetME()->ChargeHypothesis(*Xval);
    //auto end = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration ChargeHypothesis = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    auto const &measurement = QLLMatch::GetME()->Measurement();
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration Measurement = " << duration.count() << "us" << std::endl;

    //start = high_resolution_clock::now();
    Fval = QLLMatch::GetME()->QLL(hypothesis, measurement);
    //std::cout << "Fval " << Fval << std::endl;
    //end = high_resolution_clock::now();
    //duration = duration_cast<microseconds>(end - start);
    //std::cout << "Duration QLL = " << duration.count() << "us" << std::endl;

    QLLMatch::GetME()->Record(Xval[0]);
    QLLMatch::GetME()->OneStep(Xval[0]);

    return;
  }

  double QLLMatch::CallMinuit(const QCluster_t &tpc, const Flash_t &pmt, const bool init_x0) {

    if (_measurement.pe_v.empty()) {
      _measurement.pe_v.resize(DetectorSpecs::GetME().NOpDets(), 0.);
    }
    if (_measurement.pe_v.size() != pmt.pe_v.size()) {
      std::cout << _measurement.pe_v.size() << " " << pmt.pe_v.size() << std::endl;
      throw OpT0FinderException("PMT dimension has changed!");
    }

    if (!_penalty_threshold_v.empty() && _penalty_threshold_v.size() != pmt.pe_v.size()) {
      throw OpT0FinderException("Penalty threshold array has a different size than PMT array size!");
    }

    if (!_penalty_value_v.empty() && _penalty_value_v.size() != pmt.pe_v.size()) {
      throw OpT0FinderException("Penalty value array has a different size than PMT array size!");
    }

    _converged = false;

    //
    // Prepare PMT
    //
    double max_pe = 1.;

    // Debug: Print out expected PE spectrum
    //for(size_t i=0; i<pmt.pe_v.size(); ++i) {
    //std::cout << "PE meas: " << i << " " << pmt.pe_v[i] << std::endl;
    //}

    if (_normalize) {
      max_pe = 0;
      for (auto const &v : pmt.pe_v) if (v > max_pe) max_pe = v;
    }

    for (size_t i = 0; i < pmt.pe_v.size(); ++i)  _measurement.pe_v[i] = pmt.pe_v[i] / max_pe;

    _minimizer_record_chi2_v.clear();
    _minimizer_record_llhd_v.clear();
    _minimizer_record_x_v.clear();
    _num_steps = 0;
    _minimizer_min_x = std::numeric_limits<double>::max();
    _minimizer_max_x = -std::numeric_limits<double>::max();

    double reco_x = _vol_xmin + _offset; // FIXME
    if (!init_x0) {
      //reco_x = ((_vol_xmax - _vol_xmin) - (_raw_xmax_pt.x - _raw_xmin_pt.x)) / 2. + _vol_xmin;
      // Assume this is the right flash... in TPC0
      //bool tpc0 = (tpc.min_x() - DetectorSpecs::GetME().ActiveVolume().Min()[0]) <= (DetectorSpecs::GetME().ActiveVolume().Max()[0] - tpc.max_x());
      reco_x = _raw_xmin_pt.x - (pmt.time) * DetectorSpecs::GetME().DriftVelocity();
      if(reco_x < _vol_xmin || (reco_x + _raw_xmax_pt.x - _raw_xmin_pt.x) > _vol_xmax) {
          FLASH_INFO() << "returning invalid double" << std::endl;
          if (_minuit_ptr) delete _minuit_ptr;
          _minuit_ptr = 0;
	         return kINVALID_DOUBLE;
        }
    }
    double reco_x_err = ((_vol_xmax - _vol_xmin) - (_raw_xmax_pt.x - _raw_xmin_pt.x)) / 2.;
    double xmin = _vol_xmin;
    double xmax = (_vol_xmax - _vol_xmin) - (_raw_xmax_pt.x - _raw_xmin_pt.x) + _vol_xmin;
    FLASH_INFO() << _raw_xmax_pt.x << " " << _raw_xmin_pt.x << " " << (xmin < xmax) << " " << (xmin == xmax) << std::endl;
    FLASH_INFO() << "Running Minuit x: " << xmin << " => " << xmax
		 << " ... initial state x=" <<reco_x <<" x_err=" << reco_x_err << std::endl;


        if (!_minuit_ptr) _minuit_ptr = new TMinuit(4);
    double MinFval;
    int ierrflag, npari, nparx, istat;
    double arglist[4], Fmin, Fedm, Errdef;
    ierrflag = npari = nparx = istat = 0;

    assert(this == QLLMatch::GetME());

    _minuit_ptr->SetPrintLevel(-1);
    arglist[0] = 2.0;  // set strategy level
    _minuit_ptr->mnexcm("SET STR", arglist, 1, ierrflag);

    _minuit_ptr->SetFCN(MIN_vtx_qll);

    _minuit_ptr->DefineParameter(0, "X", reco_x, reco_x_err, xmin, xmax);

    _minuit_ptr->Command("SET NOW");

    // use Migrad minimizer

    arglist[0] = 5000;  // maxcalls
    arglist[1] = _migrad_tolerance; // tolerance*1e-3 = convergence condition
    _minuit_ptr->mnexcm("MIGRAD", arglist, 2, ierrflag);

    _converged = true;

    //arglist[0]   = 5.0e+2;
    //arglist[1]   = 1.0e-6;
    //_minuit_ptr->mnexcm ("simplex",arglist,2,ierrflag);
    FLASH_INFO() << " reco x before " << reco_x << std::endl;
    _minuit_ptr->GetParameter(0, reco_x, reco_x_err);
    FLASH_INFO() << " reco x after " << reco_x << std::endl;

    _minuit_ptr->mnstat(Fmin, Fedm, Errdef, npari, nparx, istat);

    // use this for debugging, maybe uncomment the actual minimzing function (MIGRAD / simplex calls)
    // scanning the parameter set
    //arglist[0] = 0;       // Parameter No (in our case x is the only and therefore the 0th parameter)
    //arglist[1] = 500;    // Number of points
    //arglist[2] = 0;       // Start point of scan
    //arglist[3] = 256;     // End point of scan
    //_minuit_ptr->mnexcm("scan", arglist,4, ierrflag);

    MinFval = Fmin;
    double *grad = 0;
    int nPar = 1;
    double fValue[1];
    fValue[0] = reco_x;
    // Transfer the minimization variables:
    MIN_vtx_qll(nPar, grad, Fmin, fValue, ierrflag);

    //static bool show = true;
    /*
      if(show){
      if(Fmin!=MinFval)std::cout<<"Fmin "<<Fmin<<" not equall to "<<MinFval<<std::endl;
      show=false;
      }
    */

    // Transfer the minimization variables:

    //FLASH_INFO() << " transfering variables " << reco_x << " " << MinFval << " " << Fmin << std::endl;
    _reco_x_offset = reco_x;
    _reco_x_offset_err = reco_x_err;
    _qll = MinFval;

    // Clear:
    _minuit_ptr->mnexcm("clear", arglist, 0, ierrflag);

    if (_minuit_ptr) delete _minuit_ptr;
    _minuit_ptr = 0;

    return _qll;
  }

}
#endif
