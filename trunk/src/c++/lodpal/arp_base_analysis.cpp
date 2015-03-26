//****************************************************************************
//* File:      ARP_base_analysis.cpp                                         *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*                    1.0 one-parameter model added.            yjs Nov. 00 *
//*                    1.1 covariate added.                      yjs Nov. 00 *
//*                    1.2 bad_sib_pair storage & func added.    yjs Mar. 01 *
//*                    1.3 rel_pair_map storage & func added.    yjs Mar. 01 *
//*                    1.3 multipoint, singlepoint seperated.    yjs Apr. 01 *
//*                    1.4 dsp, re-parameterization added.       yjs Apr. 01 *
//*                    1.5 evaluate & update_bound for dsp added.yjs May. 01 *
//*                    1.6 diagnostic option added.              yjs Jun. 01 *
//*                    1.7 seperated from ARPTest.               yjs Jul. 01 *
//*                    1.8 one-parameter model as default.                   *
//*                        two-parameter with covariate removed. yjs Aug. 01 *
//*                                                                          *
//* Notes:     This file implements ARP_base_analysis class.                 *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_base_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of ARP_base_analysis(non-Inline)               //
////////////////////////////////////////////////////////////////////////////

ARP_base_analysis::ARP_base_analysis(lodpal_pairs& p, cerrorstream& err)
                 : my_pairs(p), errors(err)
{
  nfe = 0;
  invalidate();
}

void
ARP_base_analysis::invalidate()
{
  my_valid_parameter_count = 0;
  my_built = false;
  my_good_bound = true;
}

void
ARP_base_analysis::build()
{
  if( parameters().valid() && built() )
    return;

  // Clears counts
  invalidate();

  assert( parameters().marker_count() == 1);

  // Markers init
  lodpal_parameters::marker_iterator mi;
  for( mi = parameters().marker_begin(); mi != parameters().marker_end(); ++mi )
  {
    mi->valid = false;
    mi->info.clear();
    mi->info.sample_adjustment(1);
  }

  set_current_marker(*parameters().marker_begin());

  if( my_pairs.pair_count() == 0 )
  {
    invalidate();
    errors << priority(warning) <<  "No pairs found!" << endl;
    return;         // Error no pairs
  }

  for( mi = parameters().marker_begin(); mi != parameters().marker_end(); ++mi )
  {
    mi->valid = true;
    ++my_valid_parameter_count;
  }

  lodpal_parameters::covariate_iterator ci;
  for( ci = parameters().covariate_begin(); ci != parameters().covariate_end(); ++ci )
  {
    ci->valid = true;
    ++my_valid_parameter_count;
  }

  if( !my_valid_parameter_count )
  {
    cout << "No valid parameters?" << endl;
    return;
  }

  my_built = true;
  parameters().my_valid = true;

  return;
}

void
ARP_base_analysis::re_estimate_parameters()
{
  size_t re_built_covariate = pairs_info().re_built_covariate();

  double max_y = parameters().covariate_parameters(re_built_covariate).info_y.max();

  double beta1_  = parameters().marker_parameters(0).beta1.value();
  double beta2_  = parameters().marker_parameters(0).beta2.value();

  double delta1_ = parameters().covariate_parameters(re_built_covariate).delta1.value();
  double delta2_ = parameters().covariate_parameters(re_built_covariate).delta2.value();

  double beta1   = beta1_ + max_y * delta1_;
  double beta2   = beta2_ + max_y * delta2_;
  double delta1  = (beta1_ - beta1) / max_y;
  double delta2  = (beta2_ - beta2) / max_y;

  if( parameters().covariate_parameters(re_built_covariate).adjust == covariate_type::dsp )
  {
    beta1  = beta1_;
    delta1 = delta1_ - beta1;
  }

  parameters().marker_parameters(0).beta1.set_value(beta1);
  parameters().marker_parameters(0).beta2.set_value(beta2);
  parameters().covariate_parameters(re_built_covariate).delta1.set_value(delta1);
  parameters().covariate_parameters(re_built_covariate).delta2.set_value(delta2);

  return;
}

double
ARP_base_analysis::re_evaluate(const vector<double>& param_estimates)
{
  if(!built())
    return std::numeric_limits<double>::quiet_NaN();

  vector<double> theta;
  for( size_t i = 0; i < param_estimates.size(); ++i )
    theta.push_back(param_estimates[i]);              

  return compute_lod_with_covariate(theta, true);
}

void
ARP_base_analysis::update_result(const vector<double>& param_est, size_t removed_pos)
{
  rel_pair& p = my_pairs.pairs_info()[removed_pos].lodpal_pair;

  double  f0 = p.prob_share(current_marker().marker, 0);
  double  f2 = p.prob_share(current_marker().marker, 2);

  double current_log_lr = compute_a_log_lr(param_est, removed_pos, false);

  my_pairs.pairs_info()[removed_pos].f0         = f0;
  my_pairs.pairs_info()[removed_pos].f2         = f2;
  my_pairs.pairs_info()[removed_pos].likelihood = current_log_lr;
  my_pairs.pairs_map()[p]                       = current_log_lr;

  return;
}

double
ARP_base_analysis::evaluate(vector<double>& theta)
{
  if(!built())
    return std::numeric_limits<double>::quiet_NaN();

  nfe += 1;

  return compute_lod_with_covariate(theta, false);
}

int
ARP_base_analysis::update_bounds(vector<double>& theta)
{
  // Unconstrained model
  if( parameters().autosomal_model().constraint == autosomal_model_type::unconstrained )
    return 0;

  if( parameters().covariate_count() )
    return update_bounds_with_covariate(theta);  
  
  return update_bounds_no_covariate(theta);
}

} // end of namespace LODPAL
} // end of namespace SAGE
