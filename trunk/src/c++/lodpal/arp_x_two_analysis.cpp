//****************************************************************************
//* File:      arp_x_two_analysis.cpp                                        *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs May. 02 *
//*                                                                          *
//* Notes:     This file implements ARP_x_two_analysis class.                *
//*                                                                          *
//* Copyright (c) 2002 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_x_two_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_x_two_analysis(non-Inline)              //
////////////////////////////////////////////////////////////////////////////

ARP_x_two_analysis::ARP_x_two_analysis(lodpal_pairs& p, cerrorstream& err)
                  : ARP_base_analysis(p, err)
{
  my_ld_mm = std::numeric_limits<double>::quiet_NaN();
  my_ld_mf = std::numeric_limits<double>::quiet_NaN();
  my_ld_ff = std::numeric_limits<double>::quiet_NaN();
}

////////////////////////////////////////////////////////////////////////////////////////

double
ARP_x_two_analysis::compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval)
{
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << pos << "   ";
#endif

  lodpal_pairs::lodpal_pair_info p_info = my_pairs.pairs_info()[pos];

  rel_pair p = p_info.lodpal_pair;

  const double  f0 = p.prob_share(current_marker().marker, 0);
  const double  f2 = p.prob_share(current_marker().marker, 2);
  const double  f1 = 1. - f0 - f2;
        double pf0 = my_pairs.prior_x_ibds()[p_info.prior_x_ibd_index].pf0;
        double pf1 = my_pairs.prior_x_ibds()[p_info.prior_x_ibd_index].pf1;
        double pf2 = my_pairs.prior_x_ibds()[p_info.prior_x_ibd_index].pf2;

  if( parameters().trait_parameters(0).pair_select == trait_parameter::contrast )
  {
    pf0 = my_pairs.get_drp_x_prob_share(p_info.prior_x_ibd_index, 0);
    pf2 = my_pairs.get_drp_x_prob_share(p_info.prior_x_ibd_index, 2);
    pf1 = 1. - pf0 - pf2;
  }

  // Use the proper beta1 depends on sexes.
  //
  double beta1   = 0.;
  double beta2   = 0.;

  size_t pi = 0;

  if( !parameters().x_linkage_model().lambda1_equal )
  {
    double beta1mm = 0.;
    double beta1mf = 0.;
    double beta1ff = 0.;

    if( parameters().use_mm_pair() )
    {
      beta1mm = theta[pi];
      ++pi;
    }
    if( parameters().use_mf_pair() )
    {
      beta1mf = theta[pi];
      ++pi;
    }
    if( parameters().use_ff_pair() )
    {
      beta1ff = theta[pi];
      ++pi;
    }

    if( my_pairs.is_mm_pair(p_info.prior_x_ibd_index) )
      beta1 = beta1mm;
    else if( my_pairs.is_mf_pair(p_info.prior_x_ibd_index) )
      beta1 = beta1mf;
    else if( my_pairs.is_ff_pair(p_info.prior_x_ibd_index) )
      beta1 = beta1ff;
    else
      cout << "Impossible!!!" << endl;
  }
  else
  {
    beta1 = theta[pi];
    ++pi;
  }

  beta2 = theta[pi];

  double lod_num   =  f0 +  f1*exp(beta1) +  f2*exp(beta2);
  double lod_denom = pf0 + pf1*exp(beta1) + pf2*exp(beta2);

  double w = p_info.lodpal_weight;
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
    cout << fp(f0,7,4) << " " << fp(f2,7,4) << " " << fp(pf0,7,4) << " " << fp(pf2,7,4) << " "
         << fp(lod_num,7,4) << " " << fp(lod_denom,7,4) << "   "
         << fp(likelihood_ratio,7,4) << " " << fp(current_log_lr,7,4) << endl;
#endif

  return current_log_lr;
}

void
ARP_x_two_analysis::encode_params(size_t in, maxfun_param_mgr& pm)
{
  const double inf = std::numeric_limits<double>::infinity();

  switch(in)
  {
    case 0:

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "beta1mm", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "beta1mf", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "beta1ff", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);
      }
      else
        pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 0.0, 0.0, inf, 0.02);

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "lambda1mm", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "lambda1mf", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "lambda1ff", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);
      }
      else
        pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(0.0), exp(0.0), inf);

      break;

    case 1:

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "beta1mm", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "beta1mf", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "beta1ff", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);
      }
      else
        pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 0.1, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "lambda1mm", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "lambda1mf", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "lambda1ff", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);
      }
      else
        pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(0.1), exp(0.0), inf);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);

      break;

    case 2:

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "beta1mm", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "beta1mf", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "beta1ff", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);
      }
      else
        pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, 1.0, 0.0, inf, 0.02);

      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, 2.0, 0.0, inf, 0.02);

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
          pm.addParameter("marker", "lambda1mm", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);
        if( parameters().use_mf_pair() )
          pm.addParameter("marker", "lambda1mf", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);
        if( parameters().use_ff_pair() )
          pm.addParameter("marker", "lambda1ff", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);
      }
      else
        pm.addParameter("marker", "lambda1", MAXFUN::Parameter::DEPENDENT, exp(1.0), exp(0.0), inf);

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::DEPENDENT, exp(2.0), exp(0.0), inf);

      break;

    case 3:

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
        {
          double prev_b1mm = previous_parameters().marker_parameters(0).beta1mm.value();
          pm.addParameter("marker", "beta1mm", MAXFUN::Parameter::INDEPENDENT, prev_b1mm, 0.0, inf, 0.02);
        }
        if( parameters().use_mf_pair() )
        {
          double prev_b1mf = previous_parameters().marker_parameters(0).beta1mf.value();
          pm.addParameter("marker", "beta1mf", MAXFUN::Parameter::INDEPENDENT, prev_b1mf, 0.0, inf, 0.02);
        }
        if( parameters().use_ff_pair() )
        {
          double prev_b1ff = previous_parameters().marker_parameters(0).beta1ff.value();
          pm.addParameter("marker", "beta1ff", MAXFUN::Parameter::INDEPENDENT, prev_b1ff, 0.0, inf, 0.02);
        }
      }
      else
      {
        double prev_b1 = previous_parameters().marker_parameters(0).beta1mm.value();
        pm.addParameter("marker", "beta1", MAXFUN::Parameter::INDEPENDENT, prev_b1, 0.0, inf, 0.02);
      }

      double prev_b2 = previous_parameters().marker_parameters(0).beta2ff.value();
      pm.addParameter("marker", "beta2", MAXFUN::Parameter::INDEPENDENT, prev_b2, 0.0, inf, 0.02);

      if( !parameters().x_linkage_model().lambda1_equal )
      {
        if( parameters().use_mm_pair() )
        {
          double prev_b1mm = previous_parameters().marker_parameters(0).beta1mm.value();
          pm.addParameter("marker", "lambda1mm", MAXFUN::Parameter::DEPENDENT, exp(prev_b1mm), exp(0.0), inf);
        }
        if( parameters().use_mf_pair() )
        {
          double prev_b1mf = previous_parameters().marker_parameters(0).beta1mf.value();
          pm.addParameter("marker", "lambda1mf", MAXFUN::Parameter::INDEPENDENT, exp(prev_b1mf), exp(0.0), inf);
        }
        if( parameters().use_ff_pair() )
        {
          double prev_b1ff = previous_parameters().marker_parameters(0).beta1ff.value();
          pm.addParameter("marker", "lambda1ff", MAXFUN::Parameter::INDEPENDENT, exp(prev_b1ff), exp(0.0), inf);
        }
      }
      else
      {
        double prev_b1 = previous_parameters().marker_parameters(0).beta1mm.value();
        pm.addParameter("marker", "lambda1", MAXFUN::Parameter::INDEPENDENT, exp(prev_b1), exp(0.0), inf);
      }

      pm.addParameter("marker", "lambda2", MAXFUN::Parameter::INDEPENDENT, exp(prev_b2), exp(0.0), inf);

      break;
  }

  return;
}

void
ARP_x_two_analysis::decode_theta(const vector<double>& params_estimate,
                                 const vector<double>& params_first_derivative)
{
  size_t pi = 0;

  lodpal_parameters::marker_iterator mi = parameters().marker_begin();
  if( !parameters().x_linkage_model().lambda1_equal )
  {
    if( parameters().use_mm_pair() )
    {
      mi->beta1mm.set_value(params_estimate[pi]);
      mi->beta1mm.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().use_mm_pair() )
    {
      mi->beta1mf.set_value(params_estimate[pi]);
      mi->beta1mf.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().use_mm_pair() )
    {
      mi->beta1ff.set_value(params_estimate[pi]);
      mi->beta1ff.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
  }
  else
  {
    mi->beta1mm.set_value(params_estimate[pi]);
    mi->beta1mm.set_first_derivative(params_first_derivative[pi]);
    ++pi;
  }

  mi->beta2ff.set_value(params_estimate[pi]);
  mi->beta2ff.set_first_derivative(params_first_derivative[pi]);
  ++pi;

  if( !parameters().x_linkage_model().lambda1_equal )
  {
    if( parameters().use_mm_pair() )
    {
      mi->lambda1mm.set_value(params_estimate[pi]);
      mi->lambda1mm.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().use_mm_pair() )
    {
      mi->lambda1mf.set_value(params_estimate[pi]);
      mi->lambda1mf.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
    if( parameters().use_mm_pair() )
    {
      mi->lambda1ff.set_value(params_estimate[pi]);
      mi->lambda1ff.set_first_derivative(params_first_derivative[pi]);
      ++pi;
    }
  }
  else
  {
    mi->lambda1mm.set_value(params_estimate[pi]);
    mi->lambda1mm.set_first_derivative(params_first_derivative[pi]);
    ++pi;
  }

  mi->lambda2ff.set_value(params_estimate[pi]);
  mi->lambda2ff.set_first_derivative(params_first_derivative[pi]);

  return;
}

double
ARP_x_two_analysis::compute_lod_no_covariate(vector<double>& theta, bool re_eval)
{
  KahanAdder<double> lod    = 0.0;
  KahanAdder<double> lod_mm = 0.0;
  KahanAdder<double> lod_mf = 0.0;
  KahanAdder<double> lod_ff = 0.0;

  size_t valid_pairs = 0;

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
  {
    size_t p = 0;
    cout << endl;
    if( !parameters().x_linkage_model().lambda1_equal )
    {
      if( parameters().use_mm_pair() )
      {
        cout << "beta1mm = " << theta[p];
        ++p;
      }
      if( parameters().use_mf_pair() )
      {
        cout << ", beta1mf = " << theta[p];
        ++p;
      }
      if( parameters().use_ff_pair() )
      {
        cout << ", beta1ff = " << theta[p];
        ++p;
      }
    }
    else
    {
      cout << "beta1all = " << theta[p];
      ++p;
    }

    cout << ", beta2ff = " << theta[p] << endl;

    cout << "index	type	";
    cout << "d1*y	f0	f2	pf0	pf2	lod_num	lod_denom lr	   cllr" << endl;
  }
#endif

  // Iterate over all pairs
  for( size_t pair = 0 ; pair < my_pairs.pairs_info().size(); ++pair )
  {
    lodpal_pairs::lodpal_pair_info p_info = my_pairs.pairs_info()[pair];

    if( p_info.removed_x )
    {
#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << "skippling..." << endl;
#endif
      continue;
    }

#if 0
  if( re_eval && current_marker().marker == parameters().diagnostic_marker() )
    cout << setw(5) << my_pairs.prior_x_ibds()[p_info.prior_x_ibd_index].subtype << "   ";
#endif

    ++valid_pairs;

    double current_log_lr = compute_a_log_lr(theta, pair, re_eval);

    if( !p_info.removed )
      lod  += current_log_lr;

    if( re_eval )
    {
      rel_pair& p = p_info.lodpal_pair;
      double   f0 = p.prob_share(current_marker().marker, 0);
      double   f2 = p.prob_share(current_marker().marker, 2);

      my_pairs.pairs_info()[pair].f0 = f0;
      my_pairs.pairs_info()[pair].f2 = f2;
      my_pairs.pairs_info()[pair].likelihood = current_log_lr;
      my_pairs.pairs_map()[p] = current_log_lr;

      if( !p_info.removed )
        if( my_pairs.is_mm_pair(p_info.prior_x_ibd_index) )
          lod_mm  += current_log_lr;
        else if( my_pairs.is_mf_pair(p_info.prior_x_ibd_index) )
          lod_mf  += current_log_lr;
        else if( my_pairs.is_ff_pair(p_info.prior_x_ibd_index) )
          lod_ff  += current_log_lr;
        else
          cout << "Can't be!!" << endl;
    }
  }

//  assert(valid_pairs >= my_pairs.pair_count());

  if( re_eval )
  {
    my_ld_mm = lod_mm;
    my_ld_mf = lod_mf;
    my_ld_ff = lod_ff;
  }

  return lod;
}

double
ARP_x_two_analysis::compute_lod_with_covariate(vector<double>& theta, bool re_eval)
{
  return compute_lod_no_covariate(theta, re_eval);
}

void
ARP_x_two_analysis::compute_lambda(vector<double>& theta)
{
  size_t pi = 0;

  double beta1mm = 0.;
  double beta1mf = 0.;
  double beta1ff = 0.;
  double beta1   = 0.;

  if( !parameters().x_linkage_model().lambda1_equal )
  {
    if( parameters().use_mm_pair() )
    {
      beta1mm = theta[pi];
      ++pi;
    }
    if( parameters().use_mf_pair() )
    {
      beta1mf = theta[pi];
      ++pi;
    }
    if( parameters().use_ff_pair() )
    {
      beta1ff = theta[pi];
      ++pi;
    }
  }
  else
  {
    beta1 = theta[pi];
    ++pi;
  }

  double beta2 = theta[pi];
  ++pi;

  if( !parameters().x_linkage_model().lambda1_equal )
  {
    if( parameters().use_mm_pair() )
    {
      theta[pi] = beta1mm;
      ++pi;
    }
    if( parameters().use_mf_pair() )
    {
      theta[pi] = beta1mf;
      ++pi;
    }
    if( parameters().use_ff_pair() )
    {
      theta[pi] = beta1ff;
      ++pi;
    }
  }
  else
  {
    theta[pi] = beta1;
    ++pi;
  }

  theta[pi] = beta2;

  return;
}

int
ARP_x_two_analysis::update_bounds_no_covariate(vector<double>& theta)
{
  compute_lambda(theta);

  // two-parameter(lambda2 not fixed)
  double beta1ff = 0.;

  size_t pi = 0;
  if( !parameters().x_linkage_model().lambda1_equal )
  { 
    if( parameters().use_mm_pair() && theta[pi++] < 0 )
      return 1;

    if( parameters().use_mf_pair() && theta[pi++] < 0 )
      return 1;

    if( parameters().use_ff_pair() )
    {
      beta1ff = theta[pi++];
      if( beta1ff < 0 )
        return 1;
    }
  }
  else
  {
    beta1ff = theta[pi++];
    if( beta1ff < 0 )
      return 1;
  }
    
  double beta2ff = theta[pi++];

  double beta2l = log(2*exp(beta1ff)-1);

  if( beta2ff < 0 || beta2ff < beta2l )
    return 1;

  return 0;
}

int
ARP_x_two_analysis::update_bounds_with_covariate(vector<double>& theta)
{
  return update_bounds_no_covariate(theta);
}

} // end of namespace LODPAL
} // end of namespace SAGE
