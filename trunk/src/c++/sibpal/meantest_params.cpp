//=============================================================================
// File:    mean_test_params.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Took it out from sibmeantest.h                          yjs  Jun.05
//
// Notes:   This file contains implementation for following data structures.
//            class mean_estimate
//            struct marker_parameter
//            struct trait_parameter
//            class meantest_parameters
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/meantest_params.h"

using namespace std;

#define DEBUG_MEAN(x)

namespace SAGE   {
namespace SIBPAL {

//////////////////////////////////////////////////////////////////////////

mean_estimate::mean_estimate()
{
  clear();
}

double mean_estimate::f(size_t n) const
{
  if( !finite(my_pi) || !finite(my_f1))
    return numeric_limits<double>::quiet_NaN();

  switch(n)
  {
    case 0: return 1.0 - my_pi - (1.-my_w)*my_f1;
    case 1: return my_f1;
    case 2: return my_pi - my_w*my_f1;
  }
  return numeric_limits<double>::quiet_NaN();
}

double mean_estimate::residual_variance() const
{
  return my_ss(0,0);
}

double mean_estimate::standard_error() const
{
  double ss = residual_variance();
  if(finite(ss) && ss >= 0.0)
    return sqrt(ss);
  return numeric_limits<double>::quiet_NaN();
}

double mean_estimate::residual_variance(size_t n)      const
{
  if( !finite(my_ss(0,0)) || !finite(my_ss(0,1)) || !finite(my_ss(1,1)) )
    return numeric_limits<double>::quiet_NaN();

  switch(n) 
  {
    case 0: return my_ss(0,0)+(my_w-1.)*(my_w-1.)*my_ss(1,1)-(2.)*(my_w-1.)*my_ss(0,1);
    case 1: return my_ss(1,1);
    case 2: return my_ss(0,0)+(my_w)*(my_w)*my_ss(1,1)      -(2.)*(my_w)*my_ss(0,1);
  }
  return numeric_limits<double>::quiet_NaN();
}

double mean_estimate::standard_error(size_t n) const
{
  double ss = residual_variance(n);
  if(finite(ss) && ss >= 0.0)
    return sqrt(ss);
  return numeric_limits<double>::quiet_NaN();
}

void mean_estimate::clear()
{
  const double qNaN = std::numeric_limits<double>::quiet_NaN();
  my_ss.resize_fill(2, 2, qNaN);
  my_pi = my_f1 = qNaN;
  my_w = 0.5;
}

////////////////////////////////////////////////////////////////////////////


marker_parameter::marker_parameter() 
{ 
  clear(); 
}

marker_parameter::marker_parameter(size_t m) 
{ 
  clear(); 
  marker = m; 
}

void marker_parameter::clear()
{
  estimate.clear();
  marker = (size_t) -1;
  affecteds = -1;
  pair_count = 0;
}

double marker_parameter::t_value(double mean) const
{
  double param = estimate.pi();
  double s     = estimate.standard_error();
//  double mean  = 0.5;

  if(finite(param) && finite(s) && s >= 0.0)
    return ( param-mean ) / s;

  return numeric_limits<double>::quiet_NaN();
}

double marker_parameter::t_value(size_t n, double mean) const
{
//  if(n == 1)
//    return numeric_limits<double>::quiet_NaN();
  
  double param = estimate.f(n);
  double s     = estimate.standard_error(n);
//  double mean  = 0.25;

//  if(n==1)
//    mean = 0.5;

  if(finite(param) && finite(s) && s >= 0.0)
    return ( param-mean ) / s;

  return numeric_limits<double>::quiet_NaN();
}

double marker_parameter::p_value(long df, double mean) const
{
  // P( T_df > |mean| )

  // when testing all pairs we do a two-sided test
  if( affecteds < 0 || affecteds > 2)
  {
    double t = t_value(mean);

    if( !finite(t) )
      return numeric_limits<double>::quiet_NaN();
    return T_mean_test(df,t);
  }

  double t = t_value(mean);

  if( !finite(t) )
    return numeric_limits<double>::quiet_NaN();

  // return one sided test in the correct direction
  double dir = -1;
  if( affecteds == 1 )
    dir *= -1;

  return T_cdf(df, dir*t);
}


double marker_parameter::p_value(size_t n, long df, double mean) const
{
  // P( T_df > |f_n| )
  // when testing all pairs we do a two-sided test
  if( affecteds < 0 || affecteds > 2 )
  {
    double t = t_value(n, mean);

    if( !finite(t) )
      return numeric_limits<double>::quiet_NaN();

    return T_mean_test(df,t);
  }

  double t = t_value(n, mean);

  if( !finite(t) )
    return numeric_limits<double>::quiet_NaN();

  // return one sided test in the correct direction
  double dir = 1;
  if( affecteds == 1 )
    dir *= -1;
  if( n == 2 )
    dir *= -1;
 
  return T_cdf(df, dir*t);
}

//////////////////////////////////////////////////////////////////////////

trait_parameter::trait_parameter() 
{ 
  clear(); 
}

trait_parameter::trait_parameter(size_t t, int a)
{
  clear();
  trait = t;
  if(a >= 0)
    set_affection(a);
}

void trait_parameter::set_affection(int a)
{
  if(a >= 0 && a < 3)
  {
    affection[0] = affection[1] = affection[2] = false;
    affection[a] = true;
  }
  else
    affection[0] = affection[1] = affection[2] = true;
}
      
void trait_parameter::clear()
{
  trait = (size_t) -1;
  affection[0] = affection[1] = affection[2] = true;
}

//////////////////////////////////////////////////////////////////////////

meantest_parameters::meantest_parameters()
{
  invalidate();

  my_pval_s_notation = false;
  my_export_output   = false;

  my_use_full_sibs = false;
  my_use_half_sibs = false;

  my_w = 0.5;
}

void
meantest_parameters::clear_marker_parameters()
{
  my_betas.clear();
  invalidate();
}

std::pair<meantest_parameters::marker_parameter_const_iterator, bool> 
meantest_parameters::add_marker_parameter(const marker_parameter& p)
{
  marker_parameter_iterator i = std::find(my_betas.begin(),
                                          my_betas.end(), p);

  bool inserted = true;
  if( i == my_betas.end() )
  {
    my_betas.push_back(p);
    i = my_betas.end();
    --i;
  }
  else
  {
    inserted = false;
    *i = p;
  }
  invalidate();

  return std::pair<marker_parameter_const_iterator,bool>(i,inserted);
}

std::pair<meantest_parameters::marker_parameter_const_iterator, bool>
meantest_parameters::set_marker_parameter(const marker_parameter& p)
{
  clear_marker_parameters();
  return add_marker_parameter(p);
}

std::pair<meantest_parameters::marker_parameter_const_iterator, bool>
meantest_parameters::add_marker(size_t m)
{
  return add_marker_parameter( marker_parameter(m) );
}

std::pair<meantest_parameters::marker_parameter_const_iterator, bool>
meantest_parameters::set_marker(size_t m)
{
  return set_marker_parameter( marker_parameter(m) );
}

std::pair<meantest_parameters::trait_parameter_const_iterator, bool>
meantest_parameters::add_subset(const trait_parameter& t)
{
  trait_parameter_iterator i = std::find(my_subsets.begin(),
                                         my_subsets.end(), t);

  bool inserted = true;
  if( i == my_subsets.end() )
  {
    my_subsets.push_back(t);
    i = my_subsets.end();
    --i;
  }
  else
  {
    inserted = false;
    *i = t;
  }
  invalidate();

  return std::pair<trait_parameter_const_iterator,bool>(i,inserted);
}

std::pair<meantest_parameters::trait_parameter_const_iterator, bool>
meantest_parameters::add_subset(size_t t)
{
  return add_subset( trait_parameter(t) );
}

std::pair<meantest_parameters::trait_parameter_const_iterator, bool>
meantest_parameters::add_subset(size_t t, int a)
{
  return add_subset( trait_parameter(t, a) );
}

} // end of namespace SIBPAL
} // end of namespace SAGE
