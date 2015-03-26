//============================================================================
// File:      model.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   5/15/01 - created.                         djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              model
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/model.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{


//============================================================================
// IMPLEMENTATION:  model
//============================================================================
//

set<FPED::MemberConstPointer> model::cond_mem_set; // due to JA for member of conditioned subset

// - Set parameters to their default values.
//
void
model::reset(cerrorstream& errors)
{
  // - Set model members to their default values.
  //
  mean_sub_model        = genotype_specific_mean_sub_model(errors);
  susc_sub_model        = genotype_specific_susceptibility_sub_model(errors);
  var_sub_model         = genotype_specific_variance_sub_model(&(this->mean_sub_model), errors);
  freq_sub_model        = genotype_frequency_sub_model(errors);
  resid_sub_model       = residual_correlation_sub_model(errors);
  transm_sub_model      = transmission_sub_model(&(this->freq_sub_model), errors);
  transf_sub_model      = MAXFUN::TransformationSubmodel(errors);
  mean_cov_sub_model    = MeanCovariateSubmodel(&(this->mean_sub_model), errors);
  var_cov_sub_model     = VarianceCovariateSubmodel(&(this->mean_sub_model), errors);
  susc_cov_sub_model    = SusceptibilityCovariateSubmodel(&(this->susc_sub_model), errors);
  comp_trait_sub_model  = CompositeTraitSubmodel(&(this->mean_sub_model), errors);
  fpmm_sub_model        = finite_polygenic_mixed_model_sub_model(errors);
  ons_sub_model         = onset_sub_model(errors);
  ascer_sub_model       = ascertainment_sub_model(errors);
  prev_sub_model        = prevalence_sub_model(&mean_sub_model,
                                               &susc_sub_model,
                                               &var_sub_model,
                                               &mean_cov_sub_model,
                                               &susc_cov_sub_model,
                                               &var_cov_sub_model,
                                               &freq_sub_model,
                                               &fpmm_sub_model,
                                               &ons_sub_model,
                                               &transf_sub_model,
                                               false,
                                               false);

  file_name_root  = "segreg";
  
  m_class          = model_D; 
  each_pedigree    = false;
  pen_func_output  = false;
  type_prob        = false;
  mean_missing     = true;
  susc_missing     = true;
  trans_missing    = true;
  
  primary_trait  = string();
  primary_trait_type = pt_NONE;
}

// --- GENERAL INTER-SUB-MODEL INCONSISTENCIES ---
//
// - Two or three types must be specified in the genotype specific mean
//   sub-model if transmission sub-model is specified (transmission sub-
//   model missing replaced by transmission sub-model not no_trans or homog_
//   no_trans per gcw 1-17-02).
//
bool  
model::transm_one_type() const
{
  // - If mean sub-model missing, commingling analysis is being requested
  //   and transmission sub-block has been ignored so there is no need to 
  //   invalidate the model.
  //
  if(get_type_missing()) return false;

  // Otherwise, if the transmission model is no_trans or homog_no_trans, the
  // check is irrelevant.
  //
  if(transm_sub_model.option() == transmission_sub_model::no_trans      ||
     transm_sub_model.option() == transmission_sub_model::homog_no_trans  ) 
    return false;

  // Finally, check to see if there's only one type.
  //
  return (type_dependent_sub_model().option() == genotype_specific_mean_sub_model::one);
}

// - For "homogeneous" options in transmission sub-model, hwe must be
//   the genotype frequency sub-model option.
//
bool  
model::transm_nhwe() const
{
  return ((transm_sub_model.option() == transmission_sub_model::homog_no_trans   || 
           transm_sub_model.option() == transmission_sub_model::homog_mendelian  ||
           transm_sub_model.option() == transmission_sub_model::homog_general      ) &&
           freq_sub_model.option()   == genotype_frequency_sub_model::nhwe             );  
}

/// Returns \c true if the primary trait is listed in any of the covariate
/// lists, \c false otherwise.
bool model::prim_cov() const
{
  return (comp_trait_sub_model.has_covariate(primary_trait) ||
          mean_cov_sub_model.has_covariate(primary_trait)   ||
          var_cov_sub_model.has_covariate(primary_trait)    ||
          susc_cov_sub_model.has_covariate(primary_trait)     );
}

// - Penetrance function output file only relevant if transmission sub-model
//   option is homogeneous mendelian.
//
bool  
model::transm_pen_out() const
{
  return (pen_func_output &&
          transm_sub_model.option() != transmission_sub_model::homog_mendelian);
}

// - Option to calculate type probabilities only relevant if two or three means
//   are specified in genotype specific mean sub-model.
//
bool  
model::type_prob_one_type() const
{
  return (type_prob &&
          type_dependent_sub_model().option() == genotype_specific_mean_sub_model::one);
}

// - genotype frequency sub-model option must be NONE if there is only one mean.
//
bool  
model::freq_not_none_one_type() const
{
  // If type is missing, then it's a commingling analysis, and freq is HWE
  // or NONE
  //
  if(get_type_missing()) return false;

  // If more than one type, no problem.
  //
  if(type_dependent_sub_model().option() != genotype_specific_mean_sub_model::one)
    return false;

  // Check for freq being NONE
  //
  return (freq_sub_model.option() != genotype_frequency_sub_model::NONE);
}

// - SCR 190. Geno_freq sm option may not be nhwe if there are two means.
//
bool
model::freq_nhwe_two_types() const
{
  return (freq_sub_model.option() == genotype_frequency_sub_model::nhwe &&
          type_dependent_sub_model().is_two_option());
}

// - SCR 190.  Transmission sm option may not be no trans if there are two types.
//
bool
model::transm_no_trans_two_types() const
{
  return (transm_sub_model.option() == transmission_sub_model::no_trans &&
          type_dependent_sub_model().is_two_option());
}

/// This function determines if we have a binary model with more than one
/// susceptibility, no transmission and no residuals.  Such models are 
/// equvalent to the one susceptibility model, due to bernoulli theory.
bool
model::binary_nt_no_residuals()    const
{
  // We can ignore anything that's not MLM

  if(m_class != model_MLM) return false;

  // Test for commingling analysis

  if(get_type_missing())
  {
    // A comingling analysis is always homog_no_trans, so we only need to
    // test for residuals

    return !resid_sub_model.has_residuals();
  }

  // Test for a one mean model.  One mean models are valid.

  bool one_mean_model = susc_sub_model.option() == genotype_specific_mean_sub_model::one;

  if(one_mean_model) return false;

  // Test the transmission model.  If it's not a no trans model, we're ok.

  if(trans_missing) return false;

  bool transm_nt = 
        transm_sub_model.option() == transmission_sub_model::homog_no_trans  ||
        transm_sub_model.option() == transmission_sub_model::no_trans;

  if(!transm_nt) return false;

  // Determine if there are residuals

  return !resid_sub_model.has_residuals();
}

/// Tests for covariate exclusivity.  There are four kinds of covariates,
/// mean, susceptibility, variance and composite trait.  Depending on the model,
/// a trait must not appear in multiple lists as follows:
///
/// * In continuous models, mean and composite trait are exclusive.
/// * In binary models, there are no exclusivity requirements (only
///   susceptibility is used)
/// * In age of onset models, mean is exclusive with all other covariates,
///   susceptibility is also exclusive with composite trait.  Only the
///   pairs variance/susceptibility and variance/composite trait are
///   non-exclusive.
///
/// This function tests each appropriate CovariateSubmodel pair and
/// sets the appropriate failure bits in the incon bitset as required.
void model::test_covariate_exclusivity(bitset<i_num_i>& incon) const
{
  switch(get_primary_trait_type())
  {
    case pt_CONTINUOUS :
      if(!are_covariates_exclusive(mean_cov_sub_model, comp_trait_sub_model))
        incon.set(i_mc_cov_non_exclusive);
      break;

    case pt_ONSET      :
    {
      // Check all exclusivities.  We do the checks separately to make
      // it possible for all errors to be reported.
      if(!are_covariates_exclusive(mean_cov_sub_model, comp_trait_sub_model))
        incon.set(i_mc_cov_non_exclusive);
      if(!are_covariates_exclusive(mean_cov_sub_model, susc_cov_sub_model))
        incon.set(i_ms_cov_non_exclusive);
      if(!are_covariates_exclusive(mean_cov_sub_model, var_cov_sub_model))
        incon.set(i_mv_cov_non_exclusive);
      if(!are_covariates_exclusive(susc_cov_sub_model, comp_trait_sub_model))
        incon.set(i_sc_cov_non_exclusive);
    }

    case pt_BINARY     :
    case pt_NONE       :
      break;
  }
}

// - Return a bit field w. information about each inconsistency
//   checked.
//
bitset<i_num_i>
model::check_consistency()
{
  //lint --e{534}

  bitset<i_num_i>  incon(0);
  
  mv_types()                  ? incon.set(i_mv_types)                  : incon;
  transm_one_type()           ? incon.set(i_transm_one_type)           : incon;
  transm_nhwe()               ? incon.set(i_transm_nhwe)               : incon;
  prim_cov()                  ? incon.set(i_prim_cov)                  : incon;
  transm_pen_out()            ? incon.set(i_transm_pen_out)            : incon;
  type_prob_one_type()        ? incon.set(i_type_prob_one_type)        : incon;
  freq_not_none_one_type()    ? incon.set(i_freq_not_none_one_type)    : incon;
  freq_nhwe_two_types()       ? incon.set(i_freq_nhwe_two_types)       : incon;
  transm_no_trans_two_types() ? incon.set(i_transm_no_trans_two_types) : incon;
  binary_nt_no_residuals()    ? incon.set(i_binary_nt_no_residuals)    : incon;
  
  test_covariate_exclusivity(incon);
  
  return incon; 
}

// --- MEAN-VARIANCE SUB-MODEL INCONSISTENCIES ---
//
// - Is there an inconsistency between mean and variance sub-models in terms of number
//   of types?
//
bool  
model::mv_types() const
{
  return var_sub_model.mv_types(var_sub_model.option());
}

// --- GENERAL INTER-SUB-MODEL INCONSISTENCY FIXES ---
//
void
model::fix_freq_not_none_one_type()
{
  //lint -e{534}
  freq_sub_model.set_none();
}


}
}

