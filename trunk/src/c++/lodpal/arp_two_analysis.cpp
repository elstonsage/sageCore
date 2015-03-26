//****************************************************************************
//* File:      arp_two_analysis.cpp                                          *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Jul. 01 *
//*                                                                          *
//* Notes:     This file implements ARP_two_analysis class.                  *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_two_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_two_analysis(non-Inline)                //
////////////////////////////////////////////////////////////////////////////

ARP_two_analysis::ARP_two_analysis(lodpal_pairs& p, cerrorstream& err)
                : ARP_base_analysis(p, err)
{}

////////////////////////////////////////////////////////////////////////////////////////

double
ARP_two_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
{
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << pos << "   ";
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
  const double beta2 = theta[1];

  double lod_num   =  f0 +  f1*exp(beta1) +  f2*exp(beta2);
  double lod_denom = pf0 + pf1*exp(beta1) + pf2*exp(beta2);

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
ARP_two_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
{
  const double inf = std::numeric_limits<double>::infinity();
  
  switch(in)
  {
    case 0:
      pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);

      break;

    case 1:
      pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);

      break;

    case 2:
      pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 2.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);
      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(2.0), exp(0.0), inf);

      break;

    case 3:
      double prev_b1 = previous_parameters().marker_parameters(0).beta1.value();
      double prev_b2 = previous_parameters().marker_parameters(0).beta2.value();
      pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, prev_b1, 0.0, inf, 0.02);
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, prev_b2, 0.0, inf, 0.02);

      pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(prev_b1), exp(0.0), inf);
      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(prev_b2), exp(0.0), inf);

      break;
  }

  return;
}

void
ARP_two_analysis::decode_theta(const vector<double>& params_estimate,
                               const vector<double>& params_first_derivative)
{
  lodpal_parameters::marker_iterator mi = parameters().marker_begin();
  mi->beta1.set_value(params_estimate[0]);
  mi->beta2.set_value(params_estimate[1]);
  mi->beta1.set_first_derivative(params_first_derivative[0]);
  mi->beta2.set_first_derivative(params_first_derivative[1]);

  mi->lambda1.set_value(params_estimate[2]);
  mi->lambda2.set_value(params_estimate[3]);
  mi->lambda1.set_first_derivative(params_first_derivative[2]);
  mi->lambda2.set_first_derivative(params_first_derivative[3]);

  return;
}

double
ARP_two_analysis::compute_lod_no_covariate(vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod = 0.0;
  size_t valid_pairs = 0;

#if 0
  if( current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << "beta1 = " << theta[0] << ", beta2 = " << theta[1];

    if( re_eval )
    {
      cout << endl;
      cout << "index	f0	f2	pf0	pf2	lod_num	lod_denom lr	cllr" << endl;
    }
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

//  assert(valid_pairs >= my_pairs.pair_count());

#if 0
  if( current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << "   lod = " << lod << endl;
  }
#endif

  return lod;
}

double
ARP_two_analysis::compute_lod_with_covariate(vector<double>& theta, bool re_eval)
{
  return compute_lod_no_covariate(theta, re_eval);
}

void
ARP_two_analysis::compute_lambda(vector<double>& theta)
{
  double beta1 = theta[0];
  double beta2 = theta[1];

  theta[2] = exp(beta1);  
  theta[3] = exp(beta2);  

  return;
}

int
ARP_two_analysis::update_bounds_no_covariate(vector<double>& theta)
{
  compute_lambda(theta);

  if( parameters().autosomal_model().constraint == autosomal_model_type::unconstrained )
    return 0;

  // Two-parameter model with covariate
  double beta1 = theta[0];
  double beta2 = theta[1];

  const double beta2l = log(2*exp(beta1)-1);

  if( beta1 < 0 || beta2 < 0 || beta2 < beta2l )
    return 1;

  return 0;
}

int
ARP_two_analysis::update_bounds_with_covariate(vector<double>& theta)
{
  return update_bounds_no_covariate(theta);
}

} // end of namespace LODPAL
} // end of namespace SAGE
