//****************************************************************************
//* File:      arp_po_one_analysis.cpp                                       *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Feb. 03 *
//*                                                                          *
//* Notes:     This file implements ARP_po_one_analysis class.               *
//*                                                                          *
//* Copyright (c) 2003 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_po_one_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_po_one_analysis(non-Inline)             //
////////////////////////////////////////////////////////////////////////////

ARP_po_one_analysis::ARP_po_one_analysis(lodpal_pairs& p, cerrorstream& err)
                   : ARP_base_analysis(p, err)
{}

////////////////////////////////////////////////////////////////////////////////////////

double
ARP_po_one_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
{
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << pos << "   ";
#endif

  rel_pair& p = my_pairs.pairs_info()[pos].lodpal_pair;

  const double  f0   = p.prob_share(current_marker().marker, 0);
  const double  f1   = p.prob_share(current_marker().marker, 1);
  const double  f2   = p.prob_share(current_marker().marker, 2);
  const double  f1mp = p.prob_share(current_marker().marker, 3);
  const double  f1m  = (f1 + f1mp) / 2.0;
  const double  f1p  = (f1 - f1mp) / 2.0;
        double pf0   = p.prior_prob_share(0);
        double pf2   = p.prior_prob_share(2);

  if( parameters().trait_parameters(0).pair_select == trait_parameter::contrast )
  {
    pf0 = my_pairs.get_drp_prob_share(my_pairs.pairs_info()[pos].pair_type, 0);
    pf2 = my_pairs.get_drp_prob_share(my_pairs.pairs_info()[pos].pair_type, 2);
  }

  double pf1   = 1. - pf0 - pf2;
  double pf1m  = pf1 / 2.0;
  double pf1p  = pf1m;

  if( is_maternal_hsib(p.rels().pair) )
  {
    pf1m  = pf1;
    pf1p  = 0.;
  }
  else if( is_paternal_hsib(p.rels().pair) )
  {
    pf1p  = pf1;
    pf1m  = 0.;
  }

  size_t pi = 0;

  double beta1m = 0.;
  double beta1p = 0.;
  if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
  {
    beta1m = theta[pi];
    ++pi;
  }
  if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
  {
    beta1p = theta[pi];
    ++pi;
  }

  double delta1my = 0.0;
  double delta1py = 0.0;

  // Get covariate value of pair & adjust.
  for( size_t i = 0; i < parameters().covariate_count(); ++i )
  {
    double y = my_pairs.pairs_info()[pos].lodpal_cov[i].ad_pair_value;
    
    if( my_pairs.re_built_pairs_info() && i == my_pairs.re_built_covariate() )
      y = my_pairs.pairs_info()[pos].lodpal_cov[i].re_ad_pair_value;

    if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
    {
      delta1my += (theta[pi] * y);
      ++pi;
    }
    if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
    {
      delta1py += (theta[pi] * y);
      ++pi;
    }

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << fp(y,7,4) << " " << fp(delta1my,7,4) << " " << fp(delta1py,7,4) << " ";
#endif
  }

  double alpha    = parameters().autosomal_model().alpha;
  double lambda1m = exp(beta1m + delta1my);
  double lambda1p = exp(beta1p + delta1py);
  double lambda1  = (lambda1m + lambda1p) / 2.0;
  double lambda2  = (alpha + 1)*(lambda1) - alpha;

  double lod_num   =  f0 +  f1m*lambda1m +  f1p*lambda1p +  f2*lambda2;
  double lod_denom = pf0 + pf1m*lambda1m + pf1p*lambda1p + pf2*lambda2;

  if( !(p.is_fsib_pair() || p.is_hsib_pair()) )
  {
    lod_num   =  f0 +  f1*lambda1 +  f2*lambda2;
    lod_denom = pf0 + pf1*lambda1 + pf2*lambda2;
  }

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
    cout << fp(f0,7,4) << " " << fp(f1m,7,4) << " " << fp(f1p,7,4) << " " << fp(f2,7,4) << " "
         << fp(pf0,7,4) << " " << fp(pf1m,7,4) << " " << fp(pf1p,7,4) << " " << fp(pf2,7,4) << " "
         << fp(lod_num,7,4) << " " << fp(lod_denom,7,4) << "   "
         << fp(likelihood_ratio,7,4) << " " << fp(current_log_lr,7,4) << endl;
#endif

  return current_log_lr;
}

void
ARP_po_one_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
{
  size_t covariate_count = parameters().covariate_count();

  const double ne_inf = -std::numeric_limits<double>::infinity();
  const double inf    = std::numeric_limits<double>::infinity(); 

  switch(in)
  {
    case 0:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "delta1m", MAXFUN::Parameter::INDEPENDENT, 0.0, ne_inf, inf, 0.02);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "delta1p", MAXFUN::Parameter::INDEPENDENT, 0.0, ne_inf, inf, 0.02);
      }

      if( !covariate_count )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
      }

      break;

    case 1:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "delta1m", MAXFUN::Parameter::INDEPENDENT, 0.5, ne_inf, inf, 0.02);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "delta1p", MAXFUN::Parameter::INDEPENDENT, 0.5, ne_inf, inf, 0.02);
      }

      if( !covariate_count )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
      }

      break;

    case 2:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 1.1, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 1.1, 0.0, inf, 0.02);

      for( size_t c = 0; c < covariate_count; ++c )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "delta1m", MAXFUN::Parameter::INDEPENDENT, -0.5, ne_inf, inf, 0.02);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "delta1p", MAXFUN::Parameter::INDEPENDENT, -0.5, ne_inf, inf, 0.02);
      }

      if( !covariate_count )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(1.1), exp(0.0), inf);
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(1.1), exp(0.0), inf);
      }

      break;

    case 3:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      {
        double prev_b1m = previous_parameters().marker_parameters(0).beta1m.value();
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, prev_b1m, 0.0, inf, 0.02);
      }
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      {
        double prev_b1p = previous_parameters().marker_parameters(0).beta1p.value();
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, prev_b1p, 0.0, inf, 0.02);
      }

      for( size_t c = 0; c < covariate_count; ++c )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          double prev_d1m = previous_parameters().covariate_parameters(c).delta1m.value();
          pm.addParameter("marker", "delta1m", MAXFUN::Parameter::INDEPENDENT, prev_d1m, ne_inf, inf, 0.02);
        }
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          double prev_d1p = previous_parameters().covariate_parameters(c).delta1p.value();
          pm.addParameter("marker", "delta1p", MAXFUN::Parameter::INDEPENDENT, prev_d1p, ne_inf, inf, 0.02);
        }
      }

      if( !covariate_count )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          double prev_b1m = previous_parameters().marker_parameters(0).beta1m.value();
          pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(prev_b1m), exp(0.0), inf);
        }
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          double prev_b1p = previous_parameters().marker_parameters(0).beta1p.value();
          pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(prev_b1p), exp(0.0), inf);
        }
      }

      break;
  }  

  return;
}

void
ARP_po_one_analysis::decode_theta(const vector<double>& params_estimate,
                                  const vector<double>& params_first_derivative)
{
  size_t pi = 0;
  lodpal_parameters::marker_iterator mi = parameters().marker_begin();
  if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
  {
    mi->beta1m.set_value(params_estimate[pi]);
    mi->beta1m.set_first_derivative(params_first_derivative[pi]);
    ++pi;
  }
  if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
  {
    mi->beta1p.set_value(params_estimate[pi]);
    mi->beta1p.set_first_derivative(params_first_derivative[pi]);
    ++pi;
  }

  lodpal_parameters::covariate_iterator ci;
  for( ci = parameters().covariate_begin(); ci != parameters().covariate_end(); ++ci )
  {
    if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
    {
      ci->delta1m.set_value(params_estimate[pi]);
      ci->delta1m.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
    {
      ci->delta1p.set_value(params_estimate[pi]);
      ci->delta1p.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
  }

  if( !parameters().covariate_count() )
  {
    if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
    {
      mi->lambda1m.set_value(params_estimate[pi]);
      mi->lambda1m.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
    {
      mi->lambda1p.set_value(params_estimate[pi]);
      mi->lambda1p.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
  }
}

double
ARP_po_one_analysis::compute_lod_no_covariate(vector<double>& theta, bool re_eval)
{
  return compute_lod_with_covariate(theta, re_eval);
}

double
ARP_po_one_analysis::compute_lod_with_covariate(vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod = 0.0;
  size_t valid_pairs = 0;

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    size_t p = 0;
    if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      cout << "beta1m = " << theta[p++] << " ";
    if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      cout << "beta1p = " << theta[p++];
    cout << endl;

    for( size_t i = 0; i < parameters().covariate_count(); ++i )
    {
      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        cout << "delta1m = " << theta[p++] << " ";
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        cout << "delta1p = " << theta[p++] << " ";
      cout << endl;
    }
    cout << endl;

    cout << "index	";
    for( size_t i = 0; i < parameters().covariate_count(); ++i )
      cout << "y" << i+1 << "	";
    cout << "dm*y	dp*y	f0	f1m	f1p	f2	pf0	pf1m	pf1p	pf2	lod_num	lod_denom lr	   cllr" << endl;
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
      rel_pair& p   = my_pairs.pairs_info()[pair].lodpal_pair;
      double   f0   = p.prob_share(current_marker().marker, 0);
      double   f2   = p.prob_share(current_marker().marker, 2);
      double   f1mp = p.prob_share(current_marker().marker, 3);

      my_pairs.pairs_info()[pair].f0         = f0;
      my_pairs.pairs_info()[pair].f1mp       = f1mp;
      my_pairs.pairs_info()[pair].f2         = f2;
      my_pairs.pairs_info()[pair].likelihood = current_log_lr;
      my_pairs.pairs_map()[p]                = current_log_lr;
    }
  }

//  assert(valid_pairs >= my_pairs.pair_count());

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << "return lod = " << lod << endl;
#endif

  return lod;
}

void
ARP_po_one_analysis::compute_lambda(vector<double>& theta)
{
  size_t pi = 0;

  double beta1m = 0.;                                                         
  double beta1p = 0.;                                                         

  if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
  {
    beta1m = theta[pi];
    ++pi;
  }
  if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
  {
    beta1p = theta[pi];
    ++pi;
  }

  if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
  {
    theta[pi] = exp(beta1m);
    ++pi;
  }
  if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
  {
    theta[pi] = exp(beta1p);
    ++pi;
  }

  return;
}

int
ARP_po_one_analysis::update_bounds_no_covariate(vector<double>& theta)
{
  return update_bounds_with_covariate(theta);
}

int
ARP_po_one_analysis::update_bounds_with_covariate(vector<double>& theta)
{
  compute_lambda(theta);

  // one-parameter model with covariate
  size_t pi = 0;

  double beta1m = 0.;
  double beta1p = 0.;
  if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
  {
    beta1m = theta[pi];
    ++pi;
  }
  if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
  {
    beta1p = theta[pi];
    ++pi;
  }

  if( beta1m < 0 || beta1p < 0 )
    return 1;

  if( !parameters().covariate_count() )
    return 0;

  // If only mean-centered covariates, no constraints.
  bool min_cov_exist  = false;
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
  SampleInfo z1m;
  SampleInfo z1p;

  // Iterate over all pairs
  for( size_t pair = 0 ; pair < my_pairs.pairs_info().size(); ++pair )
  {
    double delta1my = 0.0;
    double delta1py = 0.0;

    // Get covariate value of pair & adjust.
    for( size_t i = 0; i < min_cov.size(); ++i )
    {
      double y = my_pairs.pairs_info()[pair].lodpal_cov[min_cov[i]].ad_pair_value;
      if( my_pairs.re_built_pairs_info() && min_cov[i] == my_pairs.re_built_covariate() )
        y = my_pairs.pairs_info()[pair].lodpal_cov[min_cov[i]].re_ad_pair_value;

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pi += min_cov[i];
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pi += min_cov[i];

      if( y > 0 )
      {
        if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          delta1my += (theta[pi] * y);
          ++pi;
        }
        if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          delta1py += (theta[pi] * y);
          ++pi;
        }
      }
    }
    z1m.add(delta1my);
    z1p.add(delta1py);
  }

  if(    (   parameters().autosomal_model().fixed != autosomal_model_type::maternal
          && beta1m > 0 && z1m.min() == (-beta1m))
      || (   parameters().autosomal_model().fixed != autosomal_model_type::paternal
          && beta1p > 0 && z1p.min() == (-beta1p)) )
    my_good_bound = false;

  return 0;
}

} // end of namespace LODPAL
} // end of namespace SAGE
