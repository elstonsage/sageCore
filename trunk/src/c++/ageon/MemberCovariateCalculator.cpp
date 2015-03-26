//======================================================================
//
//  File:  MemberCovariateCalculator.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================

#include "globals/SAGEConstants.h"
#include "ageon/MemberCovariateCalculator.h"

namespace SAGE {
namespace AO {

//======================================================================
//
//  MemberCovariateCalculator(...) CONSTRUCTOR
//
//======================================================================
MemberCovariateCalculator::MemberCovariateCalculator(
  const Model      & mod,
        const SAMPLING::PartitionedMemberDataSample & sample,
        int             t)
  :
  my_model           (mod),
  my_sample          (sample),
  my_analysis_type   (t),
  my_cache           (mod.GetParameterMgr(), sample)
        
{
  reset();
}

//======================================================================
//
//  MemberCovariateCalculator(...) COPY CONSTRUCTOR
//
//======================================================================
MemberCovariateCalculator::MemberCovariateCalculator(
  const MemberCovariateCalculator & other) 
  :
  my_model  (other.my_model),
  my_sample (other.my_sample),
  my_cache  (other.my_cache)
{
  copy(other);
}

//======================================================================
//
//  copy(...)
//
//======================================================================
void
MemberCovariateCalculator::copy(const MemberCovariateCalculator & other)
{
  my_analysis_type      = other.my_analysis_type;
  my_AO_mean            = other.my_AO_mean;
  my_AO_stdev           = other.my_AO_stdev;
  my_AO_var             = other.my_AO_var;
  my_genetic_suscept    = other.my_genetic_suscept;
  my_num_of_inds        = other.my_num_of_inds;
  my_suscept_intercepts = other.my_suscept_intercepts;
  my_transformed_AO     = other.my_transformed_AO;
  my_transformed_AE     = other.my_transformed_AE;
  my_transformed_vals   = other.my_transformed_vals;
  my_zero_transf        = other.my_zero_transf;
}

//======================================================================
//
//  reset()
//
//======================================================================
void
MemberCovariateCalculator::reset()
{
  my_num_of_inds = my_sample.getTotalIndividualCount();

  my_genetic_suscept .resize(my_num_of_inds,QNAN);
  my_AO_mean         .resize(my_num_of_inds,QNAN);
  my_AO_var          .resize(my_num_of_inds,QNAN);
  my_AO_stdev        .resize(my_num_of_inds,QNAN);
  my_transformed_AO  .resize(my_num_of_inds,QNAN);
  my_transformed_AE  .resize(my_num_of_inds,QNAN);

  my_zero_transf = 0.0;

  // Performance features:
  my_suscept_intercepts.resize(my_model.GetParameterMgr().getParamCount("Susceptibility intercepts"));
  my_transformed_vals.resize(150, 0.0);
  precalculate_transformed_vals();
}

//======================================================================
//
//  int update()
//
//======================================================================
int
MemberCovariateCalculator::update()
{
  // 1. Update values:

  if(calculate_AO_means            ()) return 1;
  if(calculate_AO_vars             ()) return 1;
  if(calculate_genetic_suscepts    ()) return 1;
  if(calculate_AO_stdevs           ()) return 1;
  if(precalculate_transformed_vals ()) return 1;
  if(transform_AOs                 ()) return 1;
  if(transform_AEs                 ()) return 1;
  if(calculate_zero_transf         ()) return 1;

  // 3. Return success:

  return 0;
}

//===============================================================
//
//  dumpContents()
//
//===============================================================
void 
MemberCovariateCalculator::dumpContents()
{
  cout << setw(5) << " "
       << setw(15) << left << "Genetic suscept"
       << setw(15) << left << "AO mean"
       << setw(15) << left << "AO var"
       << setw(15) << left << "AO stdev"
       << endl
       << endl;

  for(size_t i = 0; i < my_genetic_suscept.size(); ++i)
  {
    cout << setw(5) << i
         << setw(15) << my_genetic_suscept[i]
         << setw(15) << my_AO_mean[i]
         << setw(15) << my_AO_var[i]
         << setw(15) << my_AO_stdev[i]
         << endl;
  }
}


//======================================================================
//
//   calculate_genetic_suscepts()
//
//======================================================================
int 
MemberCovariateCalculator::calculate_genetic_suscepts()
{
  // 1. Intercept:

  // Copy intercepts from maxfunapi to local storage:

  for(MAXFUN::ParameterConstIterator p = my_cache.suscept_int_begin; p != my_cache.suscept_int_end; ++p)
    my_suscept_intercepts[p->getGroupIndex()] = p->getCurrentEstimate();

  // Copy intercepts from local storage to individual storage:

  for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
    my_genetic_suscept[i] = my_suscept_intercepts[my_sample.getPartition(i)];

  // 2. Covariates:

  MAXFUN::ParameterConstIterator p(my_cache.suscept_cov_param_begin);
  SAMPLING::FieldConstIterator   f(my_cache.suscept_cov_field_begin);

  for(; p != my_cache.suscept_cov_param_end; ++p, ++f)
    for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
      my_genetic_suscept[i] += p->getCurrentEstimate() * f->getAdjValue(i);

  // 3. Exponentiate:

  double e_term;

  for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
  {
    e_term = exp(my_genetic_suscept[i]);

    my_genetic_suscept[i] = e_term / (1.0 + e_term);
  }

  // 4. Return success:

  return 0;
}

//======================================================================
//
//   calculate_AO_means()
//
//======================================================================
int
MemberCovariateCalculator::calculate_AO_means()
{
  // 1. Intercept:

  double mean_intercept = my_cache.mean_int->getCurrentEstimate();

  for(size_t i = 0; i < my_AO_mean.size(); i++)
    my_AO_mean[i] = mean_intercept;

  // 2. Covariates:

  MAXFUN::ParameterConstIterator p(my_cache.mean_cov_param_begin);
  SAMPLING::FieldConstIterator         f(my_cache.mean_cov_field_begin);

  for(; p != my_cache.mean_cov_param_end; ++p, ++f)
    for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
      my_AO_mean[i] += p->getCurrentEstimate() * f->getAdjValue(i);

  // 3. Return success:

  return 0;
}

//======================================================================
//
//  calculate_AO_vars()
//
//======================================================================
int
MemberCovariateCalculator::calculate_AO_vars()
{
  // 1. Intercept:

  double var_intercept = my_cache.var_int->getCurrentEstimate();

  for(size_t i = 0; i < my_AO_mean.size(); i++)
    my_AO_var[i] = var_intercept;

  // 2. Covariates:

  MAXFUN::ParameterConstIterator p(my_cache.var_cov_param_begin);
  SAMPLING::FieldConstIterator         f(my_cache.var_cov_field_begin);

  for(; p != my_cache.var_cov_param_end; ++p, ++f)
    for(size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i)
      my_AO_var[i] += p->getCurrentEstimate() * f->getAdjValue(i);

  // 3. Check for errors:

  for(size_t i = 0; i < my_AO_var.size(); i++)
    if(my_AO_var[i] <= 0)
      return 1;

  // 4. Return success:

  return 0;
}

//======================================================================
//
//  calculate_AO_stdevs()
//
//======================================================================
int
MemberCovariateCalculator::calculate_AO_stdevs()
{
  // 1. Calculate square roots:

  for(size_t i = 0; i < my_AO_var.size(); i++)
  {
    if( isnan(my_AO_var[i]) )
      my_AO_stdev[i] = QNAN;
    else

      my_AO_stdev[i] = sqrt(my_AO_var[i]);
  }

  // 2. Return success:

  return 0;
}

//======================================================================
//
//  calculate_zero_transf()
//
//======================================================================
int 
MemberCovariateCalculator::calculate_zero_transf()
{
  // 1. Are we using truncation?

  if(!NoTruncation(my_analysis_type))
    return 0;

  // 2. I guess we are!

  else
  {
    my_zero_transf = my_transformed_vals[0];

    return 0;
  }
}

//======================================================================
//
//  precalculate_transformed_vals()
//
//======================================================================
int 
MemberCovariateCalculator::precalculate_transformed_vals()
{
  // 1. Precalculate transformed vals:

  for(size_t i = 0; i < 150; i++)
  {
    my_transformed_vals[i] = i;

    if(transform(my_transformed_vals[i]))
      return 1;
  }

  // X. Return success:

  return 0;
}

//======================================================================
//
//  transform_AOs()
//
//======================================================================
int 
MemberCovariateCalculator::transform_AOs()
{
  // 1. Loop through all individuals and transform those with AO known:

  for( size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i )
  {
    if( my_sample.getField(my_model.AO_id).isAdjValuePresent(i) )
      my_transformed_AO[i] = my_transformed_vals[(size_t)my_sample.getAdjValue(i, my_model.AO_id)];
  }

  // 2. Return success:

  return 0;
}

//======================================================================
//
//  transform_AEs()
//
//======================================================================
int
MemberCovariateCalculator::transform_AEs()
{
  // 1. Loop through all individuals and transform those with AO known:

  for( size_t i = 0; i < my_sample.getTotalIndividualCount(); ++i )
  {
    if( my_sample.getField(my_model.AE_id).isAdjValuePresent(i) )
      my_transformed_AE[i] = my_transformed_vals[(size_t)my_sample.getAdjValue(i, my_model.AE_id)];
  }

  // 2. Return success:

  return 0;
}

} // End namespace AO
} // End namespace SAGE
