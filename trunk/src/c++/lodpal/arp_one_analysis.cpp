//****************************************************************************
//* File:      arp_one_analysis.cpp                                          *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Jul. 01 *
//*                                                                          *
//* Notes:     This file implements ARP_one_analysis class.                  *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_one_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_one_analysis(non-Inline)                //
////////////////////////////////////////////////////////////////////////////

ARP_one_analysis::ARP_one_analysis(lodpal_pairs& p, cerrorstream& err)
                : ARP_base_analysis(p, err)
{}

////////////////////////////////////////////////////////////////////////////////////////

double
ARP_one_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
{
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << pos << "   " << my_pairs.pairs_info()[pos].pair_type << " ";
#endif

  rel_pair& p = my_pairs.pairs_info()[pos].lodpal_pair;

  const double  f0 = p.prob_share(current_marker().marker, 0);
  const double  f2 = p.prob_share(current_marker().marker, 2);
  const double  f1 = 1. - f0 - f2;
        double pf0 = p.prior_prob_share(0);
        double pf2 = p.prior_prob_share(2);

  if( parameters().trait_parameters(0).pair_select == trait_parameter::contrast )
  {
    pf0 = my_pairs.get_drp_prob_share(my_pairs.pairs_info()[pos].pair_type, 0);
    pf2 = my_pairs.get_drp_prob_share(my_pairs.pairs_info()[pos].pair_type, 2);
  }

  const double pf1 = 1. - pf0 - pf2;

  const double beta1 = theta[0];

  double delta1y = 0.0;

  // Get covariate value of pair & adjust.
  for( size_t i = 0; i < parameters().covariate_count(); ++i )
  {
    double y = my_pairs.pairs_info()[pos].lodpal_cov[i].ad_pair_value;
    
    if( my_pairs.re_built_pairs_info() && i == my_pairs.re_built_covariate() )
      y = my_pairs.pairs_info()[pos].lodpal_cov[i].re_ad_pair_value;

    delta1y += (theta[i+1] * y);

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << fp(y,7,4) << " ";
#endif
  }

  double alpha = parameters().autosomal_model().alpha;
  double ebeta1delta1y = ((alpha + 1)*exp(beta1 + delta1y) - alpha);

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << fp(delta1y,7,4) << " " << fp(ebeta1delta1y,7,4) << " ";
#endif

  double lod_num   =  f0 +  f1*exp(beta1 + delta1y) +  f2*ebeta1delta1y;
  double lod_denom = pf0 + pf1*exp(beta1 + delta1y) + pf2*ebeta1delta1y;

  double w = my_pairs.pairs_info()[pos].lodpal_weight;
  double likelihood_ratio = w * (lod_num/lod_denom) + (1. - w);
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
ARP_one_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
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

      if( !covariate_count )
        pm.addParameter("marker", "lambda", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);

      break;

    case 1:
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, 0.5, ne_inf, inf, 0.02);

      if( !covariate_count )
        pm.addParameter("marker", "lambda", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);

      break;

    case 2:
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, -0.5, ne_inf, inf, 0.02);

      if( !covariate_count )
        pm.addParameter("marker", "lambda", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);

      break;

    case 3:
      double prev_b = previous_parameters().marker_parameters(0).beta1.value();
      pm.addParameter("marker", "beta", MAXFUN::Parameter::INDEPENDENT, prev_b, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
      {
        double prev_d = previous_parameters().covariate_parameters(c).delta1.value();
        pm.addParameter("covariate", "delta", MAXFUN::Parameter::INDEPENDENT, prev_d, ne_inf, inf, 0.02);
      }

      if( !covariate_count )
        pm.addParameter("marker", "lambda", MAXFUN::Parameter::DEPENDENT, exp(prev_b), exp(0.0), inf);

      break;
  }

  return;
}

void
ARP_one_analysis::decode_theta(const vector<double>& params_estimate,
                               const vector<double>& params_first_derivative)
{
  lodpal_parameters::marker_iterator mi = parameters().marker_begin();
  mi->beta1.set_value(params_estimate[0]);
  mi->beta1.set_first_derivative(params_first_derivative[0]);

  if( !parameters().covariate_count() )
  {
    mi->lambda1.set_value(params_estimate[1]);
    mi->lambda1.set_first_derivative(params_first_derivative[1]);
  }
  else
  {
    size_t pi = 1;
    lodpal_parameters::covariate_iterator ci;
    for( ci = parameters().covariate_begin(); ci != parameters().covariate_end(); ++ci, ++pi )
    {
      ci->delta1.set_value(params_estimate[pi]);
      ci->delta1.set_first_derivative(params_first_derivative[pi]);
    }
  }
}

double
ARP_one_analysis::compute_lod_no_covariate(vector<double>& theta, bool re_eval)
{
  return compute_lod_with_covariate(theta, re_eval);
}

double
ARP_one_analysis::compute_lod_with_covariate(vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod = 0.0;
  size_t valid_pairs = 0;

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << "beta1 = " << theta[0] << endl;
    for( size_t i = 0; i < parameters().covariate_count(); ++i )
      cout << "delta" << i+1 << " = " << theta[i+1] << ", ";
    cout << endl;
    cout << "index	";
    for( size_t i = 0; i < parameters().covariate_count(); ++i )
      cout << "y" << i+1 << "	";
    cout << "d*y	ebdy	f0	f2	pf0	pf2	lod_num	lod_denom lr	   cllr" << endl;
  }
  else if( current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << nfe << " " << "beta1 = " << theta[0];
    for( size_t i = 0; i < parameters().covariate_count(); ++i )
      cout << ", delta" << i+1 << " = " << theta[i+1];
    cout << endl;
  }
#endif

  // Iterate over all pairs
  for( size_t pair = 0 ; pair < my_pairs.pairs_info().size(); ++pair )
  {
    if( my_pairs.pairs_info()[pair].removed )
      continue;

    ++valid_pairs;

    double current_log_lr = compute_a_log_lr(theta, pair, re_eval);

    lod  += current_log_lr;

    if( re_eval )
    {
      rel_pair& p = my_pairs.pairs_info()[pair].lodpal_pair;
      double   f0 = p.prob_share(current_marker().marker, 0);
      double   f2 = p.prob_share(current_marker().marker, 2);

      my_pairs.pairs_info()[pair].f0         = f0;
      my_pairs.pairs_info()[pair].f2         = f2;
      my_pairs.pairs_info()[pair].likelihood = current_log_lr;
      my_pairs.pairs_map()[p]                = current_log_lr;
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

void
ARP_one_analysis::compute_lambda(vector<double>& theta)
{
  theta[1] = exp(theta[0]);
  return;
}

int
ARP_one_analysis::update_bounds_no_covariate(vector<double>& theta)
{
  return update_bounds_with_covariate(theta);
}

int
ARP_one_analysis::update_bounds_with_covariate(vector<double>& theta)
{
  if( !parameters().covariate_count() )
    compute_lambda(theta);

  if( parameters().autosomal_model().constraint == autosomal_model_type::unconstrained )
    return 0;

  // one-parameter model with covariate
  double beta1 = theta[0];

  if( beta1 < 0 )
    return 1;

  if( !parameters().covariate_count() )
    return 0;

  // If only mean-centered covariates, no constraints.
  bool min_cov_exist = false;
  vector<size_t> min_cov;

  for( size_t i = 0; i < parameters().covariate_count(); ++i )
  {
    if( parameters().covariate_parameters(i).adjust == covariate_type::minimum )
    {
      min_cov_exist = true;
      min_cov.push_back(i);
    }
  }

  if( !min_cov_exist )
    return 0;

  // Check for constraints on minimum-adjusted covariates.
  SampleInfo z1;

  // Iterate over all pairs
  for( size_t pair = 0 ; pair < my_pairs.pairs_info().size(); ++pair )
  {
    double delta1y = 0.0;

    // Get covariate value of pair & adjust.
    for( size_t i = 0; i < min_cov.size(); ++i )
    {
      double y = my_pairs.pairs_info()[pair].lodpal_cov[min_cov[i]].ad_pair_value;
      if( my_pairs.re_built_pairs_info() && min_cov[i] == my_pairs.re_built_covariate() )
        y = my_pairs.pairs_info()[pair].lodpal_cov[min_cov[i]].re_ad_pair_value;

      if( y > 0 )
        delta1y += (theta[min_cov[i]+1] * y);
    }
    z1.add(delta1y);
  }

  if( z1.min() < (-beta1) )
    return 1;

  if( beta1 > 0 && z1.min() == (-beta1) )
    my_good_bound = false;

  return 0;
}

} // end of namespace LODPAL
} // end of namespace SAGE
