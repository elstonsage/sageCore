//=======================================================================
//  File:	member_covariate_data.ipp
//
//  Purpose:	This class stores information about variate and covariate
//              traits.
//
//  Author:	Stephen Gross
//
//  History:	0.1  sag  Initial Implementation	Jul 11 01
//
//  Copyright (c) 2001 R.C. Elston
//  All Rights Reserved
//=======================================================================

#ifndef MEMBER_COVARIATE_DATA_H
#include "segreg/member_covariate_data.h"
#endif

namespace SAGE   {
namespace SEGREG {

// Defines to make life easier with the enumeration.

#define missing    ascertainment_sub_model::not_specified
#define actual     ascertainment_sub_model::actual
#define gte_thresh ascertainment_sub_model::gte_thresh
#define lte_thresh ascertainment_sub_model::lte_thresh

inline
member_covariate_data::member_covariate_data()
  : Analysis_Trait(QNAN),
    my_member_type(missing)
{}

inline void
member_covariate_data::set_number_of_covariates(
  size_t Number_of_Composite_Covariates, 
  size_t Number_of_Susceptibility_Covariates, 
  size_t Number_of_Covariate_Covariates,
  size_t Number_of_Variance_Covariates)
{
  composite_covariates      .resize(Number_of_Composite_Covariates,      QNAN);
  susceptibility_covariates .resize(Number_of_Susceptibility_Covariates, QNAN);
  mean_covariates           .resize(Number_of_Covariate_Covariates,      QNAN);
  variance_covariates       .resize(Number_of_Variance_Covariates,       QNAN);
}

inline size_t member_covariate_data::get_composite_cov_count      () const { return composite_covariates      .size(); }
inline size_t member_covariate_data::get_susceptibility_cov_count () const { return susceptibility_covariates .size(); }
inline size_t member_covariate_data::get_mean_cov_count           () const { return mean_covariates           .size(); }    
inline size_t member_covariate_data::get_variance_cov_count       () const { return variance_covariates       .size(); }    

inline double
member_covariate_data::get_analysis_trait() const
{ return Analysis_Trait; }

inline double
member_covariate_data::set_analysis_trait(double value)
{
  Analysis_Trait = value;
  return Analysis_Trait;
}

inline double 
member_covariate_data::get_composite_cov(size_t TraitIndex) const
{ return composite_covariates[TraitIndex]; }

inline bool
member_covariate_data::set_composite_cov(size_t TraitIndex, double value)
{
  composite_covariates[TraitIndex] = value;
  return true;
}
    
inline double 
member_covariate_data::get_susceptibility_cov(size_t TraitIndex) const
{ return susceptibility_covariates[TraitIndex]; }

inline bool
member_covariate_data::set_susceptibility_cov(size_t TraitIndex, double value)
{
  susceptibility_covariates[TraitIndex] = value;
  return true;
}
    
inline double
member_covariate_data::get_mean_cov(size_t TraitIndex) const
{ return mean_covariates[TraitIndex]; }

inline bool
member_covariate_data::set_mean_cov(size_t TraitIndex, double value)
{
  mean_covariates[TraitIndex] = value;
  return true;
}

inline double
member_covariate_data::get_variance_cov(size_t TraitIndex) const
{ return variance_covariates[TraitIndex]; }

inline bool
member_covariate_data::set_variance_cov(size_t TraitIndex, double value)
{
  variance_covariates[TraitIndex] = value;
  return true;
}

inline void member_covariate_data::set_member_type(member_type t)
{
  my_member_type = t;
}

inline member_covariate_data::member_type
     member_covariate_data::get_member_type() const
{
  return my_member_type;
}

inline bool
member_covariate_data::validate()
{
  // If the member type is missing already, then there's nothing to do.
  if(my_member_type == missing) return false;

  // Otherwise, we check for missing values

  if(SAGE::isnan(Analysis_Trait))
  { 
    my_member_type = missing;
    return is_valid();
  }

  for(vector<double>::iterator i = composite_covariates.begin(); i != composite_covariates.end(); ++i)
    if(SAGE::isnan(*i))
    {
      my_member_type = missing;
      return is_valid();
    }

  for(vector<double>::iterator i = mean_covariates.begin(); i != mean_covariates.end(); ++i)
    if(SAGE::isnan(*i))
    {
      my_member_type = missing;
      return is_valid();
    }

  for(vector<double>::iterator i = variance_covariates.begin(); i != variance_covariates.end(); ++i)
    if(SAGE::isnan(*i))
    {
      my_member_type = missing;
      return is_valid();
    }

  for(vector<double>::iterator i  = susceptibility_covariates.begin();
                               i != susceptibility_covariates.end(); ++i)
    if(SAGE::isnan(*i))
    {
      my_member_type = missing;
      return is_valid();
    }

  return is_valid();
}

inline bool
member_covariate_data::is_valid() const
{
  return my_member_type != missing;
}

// Undefine the ennumeration

#undef missing
#undef actual
#undef gte_thresh
#undef lte_thresh

}}

