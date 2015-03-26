//=============================================================================
// File:    regress_result.cpp
//
// Author:  Kevin Jacobs 
//
// History: Version 0.0 Initial implementation
//          Took it out from TraitRegression                        yjs  Jun.09
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress_result.h"

namespace SAGE   {
namespace SIBPAL {

parameter_estimate::parameter_estimate()
{
  clear();
}

void
parameter_estimate::clear()
{
  set_test_direction(dTWOSIDED);
  set_value(QNAN);
  set_variance(QNAN);
  set_expected_value(QNAN);
  set_degrees_of_freedom(0);
  set_adjusted_tvalue(QNAN);
  set_adjusted_degrees_of_freedom(0);
  set_adjusted_pvalue(QNAN);

  clear_simulation();
}

void
parameter_estimate::clear_simulation()
{
  my_significant_replicates = my_same_replicates = my_total_replicates = 0;
}

void
parameter_estimate::add_simulation_result(double v, double ss)
{
  if( !finite(v) || !finite(ss) || ss <= 0 )
    return;

  double obs_t = tvalue();

  if( !finite(obs_t) )
    return;

  ++my_total_replicates;

  double curr_t = v / sqrt(ss);

  bool sig = false;
  bool same = false;

  switch( test_direction() )
  {
    case dTWOSIDED:  if(      fabs(curr_t) >  fabs(obs_t)) sig  = true;
                     else if( fabs(curr_t) == fabs(obs_t)) same = true; break;
    case dRIGHTSIDE: if(           curr_t  >       obs_t ) sig  = true;
                     else if(      curr_t  ==      obs_t ) same = true; break;
    case dLEFTSIDE:  if(           curr_t  <       obs_t ) sig  = true;
                     else if(      curr_t  ==      obs_t ) same = true; break;
  }

  if( sig )
  {
    ++my_significant_replicates;
  }
  else if( same )
  {
    ++my_same_replicates;
  }

  return;
}

double
parameter_estimate::pvalue() const
{
  double p = adjusted_pvalue();

  if( !SAGE::isnan(p) )
    return p;

  return raw_pvalue();
}

double
parameter_estimate::raw_pvalue() const
{
  return compute_pvalue( raw_tvalue(), degrees_of_freedom(), test_direction() );
}

double
parameter_estimate::adjusted_tvalue_significance() const
{
  size_t df = adjusted_degrees_of_freedom();

  if(!df)
    df = degrees_of_freedom();

  return compute_pvalue( adjusted_tvalue(), df, test_direction() );
}

double
parameter_estimate::adjusted_pvalue() const
{
  double p = raw_adjusted_pvalue();

  if( !SAGE::isnan(p) )
    return p;

  return adjusted_tvalue_significance();
}

double
parameter_estimate::compute_pvalue(double t, size_t df, direction_type dir) const
{
  if(SAGE::isnan(t) || df == 0)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double p = T_cdf(df, -fabs(t));

  switch( dir )
  {
    case dTWOSIDED:  return 2.0*p;
    case dLEFTSIDE:  return (t<0.0)? p : 1.0 - p;
    case dRIGHTSIDE: return (t>0.0)? p : 1.0 - p;
    default:         return std::numeric_limits<double>::quiet_NaN();
  }
}

//
//-----------------------------------------------------------
//

regression_results::regression_results()
{
  clear();
  clear_full();
  clear_half();
}

void
regression_results::clear()
{
  my_variances.residual = my_variances.total = QNAN;
  my_variances.sum_residual = my_variances.diff_residual = QNAN;

  my_residual_info.clear();
  my_residual_square_info.clear();
  my_residual_square_info_reduced.clear();

  my_results.resize(0);
  my_full_results.resize(0);
  my_half_results.resize(0);

  my_intercept_count = 0;

  my_F_result.invalidate();

  return;
}

void
regression_results::clear_full()
{
  my_full_variances.residual = my_full_variances.total = QNAN;
  my_full_variances.sum_residual = my_full_variances.diff_residual = QNAN;

  my_full_residual_info.clear();
  my_full_results.resize(0);

  return;
}

void
regression_results::clear_half()
{
  my_half_variances.residual = my_half_variances.total = QNAN;
  my_half_variances.sum_residual = my_half_variances.diff_residual = QNAN;

  my_half_residual_info.clear();
  my_half_results.resize(0);

  return;
}

regression_results::~regression_results()
{}

void
regression_results::dump(ostream &out) const
{
  out << "regression_results:" << endl;

  for( size_t i = 0; i < my_results.size(); ++i )
  {
    out << "param type = " << my_results[i].type()
        << ", pair type = " << my_results[i].get_pair_type()
        << ", estimate = " << my_results[i].estimate.value()
        << ", variance = " << my_results[i].estimate.variance()
        << endl;
  }

  if( my_full_results.size() )
  {
    out << "regression_results with full sibs:" << endl;

    for( size_t i = 0; i < my_full_results.size(); ++i )
    {
      out << "param type = " << my_full_results[i].type()
          << ", pair type = " << my_full_results[i].get_pair_type()
          << ", estimate = " << my_full_results[i].estimate.value()
          << ", variance = " << my_full_results[i].estimate.variance()
          << endl;
    }
  }

  if( my_half_results.size() )
  {
    out << "regression_results with half sibs:" << endl;

    for( size_t i = 0; i < my_half_results.size(); ++i )
    {
      out << "param type = " << my_half_results[i].type()
          << ", pair type = " << my_half_results[i].get_pair_type()
          << ", estimate = " << my_half_results[i].estimate.value()
          << ", variance = " << my_half_results[i].estimate.variance()
          << endl;
    }
  }

  out << "intercept count = " << my_intercept_count << endl;

  return;
}

} //end of namespace SIBPAL
} //end of namespace SAGE
