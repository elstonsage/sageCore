//****************************************************************************
//* File:      arp_po_two_analysis.cpp                                       *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Feb. 03 *
//*                                                                          *
//* Notes:     This file implements ARP_po_two_analysis class.               *
//*                                                                          *
//* Copyright (c) 2003 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_po_two_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_po_two_analysis(non-Inline)             //
////////////////////////////////////////////////////////////////////////////

ARP_po_two_analysis::ARP_po_two_analysis(lodpal_pairs& p, cerrorstream& err)
                   : ARP_base_analysis(p, err)
{}

////////////////////////////////////////////////////////////////////////////////////////

double
ARP_po_two_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
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
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << "m ";
#endif
  }
  else if( is_paternal_hsib(p.rels().pair) )
  {
    pf1p  = pf1;
    pf1m  = 0.;
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << "p ";
#endif
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

  double beta2  = theta[pi];

  double lambda1m = exp(beta1m);
  double lambda1p = exp(beta1p);
  double lambda1  = (lambda1m + lambda1p) / 2.0;
  double lambda2  = exp(beta2);

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
    cout << fp(f0,18,15) << " " << fp(f2,18,15) << " " << fp(f1mp,18,15) << " " << fp(pf0,7,4) << " " << fp(pf2,7,4) << " "
         << fp(lod_num,7,4) << " " << fp(lod_denom,7,4) << "   "
         << fp(likelihood_ratio,7,4) << " " << fp(current_log_lr,7,4) << endl;
#endif

  return current_log_lr;
}

void
ARP_po_two_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
{
  const double inf = std::numeric_limits<double>::infinity();

  switch(in)
  {
    case 0:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf, 0.02);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf, 0.02);

      break;

    case 1:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 0.2, 0.0, inf, 0.02);

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf, 0.02);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(0.2), exp(0.0), inf, 0.02);

      break;

    case 2:

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "beta1m", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "beta1p", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 2.0, 0.0, inf, 0.02);

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf, 0.02);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(2.0), exp(0.0), inf, 0.02);

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

      double prev_b2 = previous_parameters().marker_parameters(0).beta2.value();
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, prev_b2, 0.0, inf, 0.02);

      if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      {
        double prev_b1m = previous_parameters().marker_parameters(0).beta1m.value();
        pm.addParameter("marker", "lambda1m", MAXFUN::Parameter::DEPENDENT, exp(prev_b1m), exp(0.0), inf);
      }
      if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      {
        double prev_b1p = previous_parameters().marker_parameters(0).beta1p.value();
        pm.addParameter("marker", "lambda1p", MAXFUN::Parameter::DEPENDENT, exp(prev_b1p), exp(0.0), inf, 0.02);
      }

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(prev_b2), exp(0.0), inf, 0.02);

      break;

  }

  return;
}

void
ARP_po_two_analysis::decode_theta(const vector<double>& params_estimate,
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

  mi->beta2.set_value(params_estimate[pi]);
  mi->beta2.set_first_derivative(params_first_derivative[pi]);
  ++pi;

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

  mi->lambda2.set_value(params_estimate[pi]);
  mi->lambda2.set_first_derivative(params_first_derivative[pi]);
}

double
ARP_po_two_analysis::compute_lod_no_covariate(vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod = 0.0;
  size_t valid_pairs = 0;

#if 0
  if( current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << nfe << " ";
    size_t p = 0;
    if( parameters().autosomal_model().fixed != autosomal_model_type::maternal )
    {
      cout << "beta1m = " << theta[p] << "	";
      ++p;
    }
    if( parameters().autosomal_model().fixed != autosomal_model_type::paternal )
    {
      cout << "beta1p = " << theta[p] << "	";
      ++p;
    }
    cout << "beta2 = " << theta[p];
  }
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    cout << endl;
    cout << "index	";
    cout << "f0	f2	f1mp	pf0	pf2	lod_num	lod_denom lr	   cllr" << endl;
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
      rel_pair& p  = my_pairs.pairs_info()[pair].lodpal_pair;
      double  f0   = p.prob_share(current_marker().marker, 0);
      double  f2   = p.prob_share(current_marker().marker, 2);
      double  f1mp = p.prob_share(current_marker().marker, 3);

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
  else if( current_marker().marker == parameters().diagnostic_marker() )
    cout << "	return lod = " << lod << endl;
#endif

  return lod;
}

double
ARP_po_two_analysis::compute_lod_with_covariate(vector<double>& theta, bool re_eval)
{
  return compute_lod_no_covariate(theta, re_eval);
}

void
ARP_po_two_analysis::compute_lambda(vector<double>& theta)
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

  double beta2 = theta[pi++];

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

  theta[pi] = exp(beta2);

  return;
}

int
ARP_po_two_analysis::update_bounds_no_covariate(vector<double>& theta)
{
  compute_lambda(theta);

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

  double beta2  = theta[pi];

  const double beta2l = log(exp(beta1m)+exp(beta1p)-1.);

  if( beta1m < 0 || beta1p < 0 || beta2 < beta2l )
    return 1;

  return 0;
}

int
ARP_po_two_analysis::update_bounds_with_covariate(vector<double>& theta)
{
  return update_bounds_no_covariate(theta);
}

} // end of namespace LODPAL
} // end of namespace SAGE
