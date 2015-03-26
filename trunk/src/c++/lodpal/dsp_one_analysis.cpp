//****************************************************************************
//* File:      dsp_one_analysis.cpp                                          *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Apr. 01 *
//*                                                                          *
//* Notes:     This file implements DSP_one_analysis class.                  *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/dsp_one_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of DSP_one_analysis(non-Inline)                //
////////////////////////////////////////////////////////////////////////////

DSP_one_analysis::DSP_one_analysis(lodpal_pairs& p, cerrorstream& err)
                : ARP_one_analysis(p, err)
{}

////////////////////////////////////////////////////////////////////////////////////////

double
DSP_one_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
{
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << pos << "   ";
#endif

  rel_pair p = my_pairs.pairs_info()[pos].lodpal_pair;

  const double  f0 = p.prob_share(current_marker().marker, 0);
  const double  f2 = p.prob_share(current_marker().marker, 2);
  const double  f1 = 1. - f0 - f2;
  const double pf0 = p.prior_prob_share(0);
  const double pf2 = p.prior_prob_share(2);
  const double pf1 = 1. - pf0 - pf2;

  double y  = my_pairs.pairs_info()[pos].lodpal_cov[0].ad_pair_value;
  double y_ = 1. - y;

  double beta1y_ = theta[0] * y_;
  double delta1y = theta[1] * y;

  double alpha = parameters().autosomal_model().alpha;
  double ebeta1delta1y = ((alpha + 1)*exp(beta1y_ + delta1y) - alpha);

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << fp(y,7,4) << " " << fp(delta1y,7,4) << " " << fp(ebeta1delta1y,7,4) << " ";
#endif

  double lod_num   =  f0 +  f1*exp(beta1y_ + delta1y) +  f2*ebeta1delta1y;
  double lod_denom = pf0 + pf1*exp(beta1y_ + delta1y) + pf2*ebeta1delta1y;

  double w = my_pairs.pairs_info()[pos].lodpal_weight;
  double likelihood_ratio = w * (lod_num/lod_denom) + (1 - w);
  double current_log_lr   = 0.;

  if( !finite(likelihood_ratio) || likelihood_ratio < 0.0001 )
    current_log_lr = -4.0;
  else
    current_log_lr = log10(likelihood_ratio);

  if( !re_eval && current_log_lr > 0.45 )
    current_log_lr = 0.45; // my_pairs.pairs_map()[p];

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << fp(f0,7,4) << " " << fp(f2,7,4) << " " << fp(pf0,7,4) << " " << fp(pf2,7,4) << " " 
         << fp(lod_num,7,4) << " " << fp(lod_denom,7,4) << "   "
         << fp(likelihood_ratio,7,4) << " " << fp(current_log_lr,7,4) << endl;
  }
#endif

  return current_log_lr;
}

void
DSP_one_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
{
  size_t covariate_count = parameters().covariate_count();

  const double ne_inf = -std::numeric_limits<double>::infinity();
  const double inf    = std::numeric_limits<double>::infinity();

  switch(in)
  {
    case 0:
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, 0.0, ne_inf, inf, 0.02);

      break;

    case 1:
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, -0.5, ne_inf, inf, 0.02);

      break;

    case 2:
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, -1.5, ne_inf, inf, 0.02);

      break;

    case 3:
      double prev_b = previous_parameters().marker_parameters(0).beta1.value() * 0.9;
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, prev_b, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
      {
        double prev_d = previous_parameters().covariate_parameters(c).delta1.value() * 0.9;
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, prev_d, ne_inf, inf, 0.02);
      }

      break;
  }

  return;
}

double
DSP_one_analysis::compute_lod_with_covariate(const vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod = 0.0;
  size_t valid_pairs = 0;

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << "beta1 = " << theta[0] << ", delta1 = " << theta[1] << endl;
    cout << "index	y	d1*y	ebdy	f0	f2	pf0	pf2	lod_num	lod_denom lr	   cllr" << endl;
  }
#endif

  // Iterate over all pairs
  for( size_t pair = 0 ; pair < my_pairs.pairs_info().size(); ++pair )
  {
    ++valid_pairs;

    double current_log_lr = compute_a_log_lr(theta, pair, re_eval);

    if( !my_pairs.pairs_info()[pair].removed )
      lod  += current_log_lr;

    if( re_eval )
    {
      rel_pair p = my_pairs.pairs_info()[pair].lodpal_pair;
      double  f0 = p.prob_share(current_marker().marker, 0);
      double  f2 = p.prob_share(current_marker().marker, 2);

      my_pairs.pairs_info()[pair].f0 = f0;
      my_pairs.pairs_info()[pair].f2 = f2;
      my_pairs.pairs_info()[pair].likelihood = current_log_lr;
      my_pairs.pairs_map()[p] = current_log_lr;
    }
  }

  assert(valid_pairs >= my_pairs.pair_count());

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << "return lod = " << lod << endl;
  }
#endif

  return lod;
}

int
DSP_one_analysis::update_bounds_with_covariate(const vector<double>& theta)
{
  // one-parameter model with covariate
  double beta1  = theta[0];
  double delta1 = theta[1];

  my_good_bound = false;

//  if( beta1 < 0 || delta1 > 0 )
  if( beta1 < 0 || delta1 > -beta1 )
    return 1;

  return 0;
}

} // end of namespace LODPAL
} // end of namespace SAGE
