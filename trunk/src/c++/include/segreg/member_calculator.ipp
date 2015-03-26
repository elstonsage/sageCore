//=========================================================================
//  File:	member_calculator.ipp
//
//  Author:	Stephen Gross
//
//  History:	0.1  sag  Initial implementation	Jul 11 01
//
//  Notes:	Calculates analysis traits, expected means, and expected
//		variates.
//
//  Copyright (c) 2001 R. C. Elston
//=========================================================================

#ifndef MEMBER_CALCULATOR_H
#include "segreg/member_calculator.h"
#endif

namespace SAGE {
namespace SEGREG {

inline size_t 
member_calculator_base::get_abs_mem_ref(const member_type& mem) const
{ 
  if(mem.pedigree() != last_pedigree_accessed)
  {
    last_pedigree_accessed = mem.pedigree();
    last_pedigree_location = pedigree_index_map.find(mem.pedigree())->second;
  }

  return last_pedigree_location + mem.index();
}

inline bool member_calculator_base::is_member_valid
    (const member_type& mem)  const
{
  return is_member_valid(get_abs_mem_ref(mem));
}

inline bool member_calculator_base::is_member_valid
    (size_t i)  const
{
  return get_member_class(i) != missing;
}

inline member_calculator_base::member_class
    member_calculator_base::get_member_class
      (const member_type& mem)  const
{
  return get_member_class(get_abs_mem_ref(mem));
}

inline member_calculator_base::member_class
    member_calculator_base::get_member_class
      (size_t i)  const
{
  return member_classes[i];
}

inline
bool member_calculator_base::import_trait
    (size_t                 trait_index,
     pedigree_const_pointer ped,
     size_t                 member_index,
     size_t                 abs_mem_ref,
     double&                val)
{
  // Get the status for the member

  bool b = ped->info().trait_missing(member_index,trait_index);

  if(b)  // If trait is missing.
  {
    member_classes[abs_mem_ref] = missing;
  }
  else
  {
    val = ped->info().trait(member_index,trait_index);
  }

  return !b;
}

//lint -e{601} <- spurious errors caused by vector<bool>
inline
bool member_calculator_base::import_trait
    (size_t                  trait_index,
     pedigree_const_pointer  ped,
     size_t                  member_index,
     size_t                  abs_mem_ref,
     vector<bool>::reference val)
{
  // Get the status for the member

  bool b = ped->info().trait_missing(member_index,trait_index);

  if(b)  // If trait is missing.
  {
    member_classes[abs_mem_ref] = missing;
  }
  else
  {
    val = (bool) ped->info().trait(member_index,trait_index);
  }

  //lint --e{550} <- Lint thinks val not used.

  return !b;
}

/// a mini-function for centering a vector of doubles around a mean.
/// This is used in all the member calculators.

inline void
member_calculator_base::center_trait(vector<double>& v, double mean) const
{
  for(size_t i = 0; i < v.size(); ++i)
    v[i] = v[i] - mean;
}


//=====================================================================================
//
//  update()
//
//=====================================================================================

inline int continuous_member_calculator::update()
{
  // First perform all the calculations on those that might produce errors
  // Note that while currently some of these don't ever produce errors,
  // the error check is left in as a safety measure.

  int err_code = 0;

  err_code = calculate_composite_traits  (); if(err_code) return err_code;
  err_code = calculate_expected_means    (); if(err_code) return err_code;
  err_code = calculate_expected_variances(); if(err_code) return err_code;
  err_code = calculate_expected_sds      (); if(err_code) return err_code;
  err_code = calculate_standardizations  (); if(err_code) return err_code;

  // If we're using ascertainment, we calculate the ascertained
  // standardizations.  This uses several possible values for ti. 
  // Otherwise, the ascertained standardizations are the same as the normal
  // standardizations.

  if(use_ascertainment)
  {
    err_code = calculate_ascertained_standardizations();
    if(err_code) return err_code;
  }
  else
  {
    ascertained_standardizations = standardizations;
  }

  // Calculate model specific elements.

  if(mod.get_model_class() == model_D)
  {
    err_code = calculate_estimated_standardizations();
    if(err_code) return err_code;
  }
  else if(mod.get_model_class() == model_FPMM)
  {
    err_code = calculate_expected_polygenic_means();
    if(err_code) return err_code;
  }

  return segreg_errors::EVAL_OK;
}

//=====================================================================================
//
//  get_composite_trait(...)
//
//=====================================================================================

inline double 
continuous_member_calculator::get_composite_trait(const member_type& mem) const
{
  return get_composite_trait(get_abs_mem_ref(mem));
}

inline double
continuous_member_calculator::get_composite_trait(size_t abs_mem_ref) const
{
  return composite_traits[abs_mem_ref];
}

//=====================================================================================
//
//  get_expected_mean(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_expected_mean(const member_type& mem,
                                                genotype_index         genotype) const
{
  return get_expected_mean(get_abs_mem_ref(mem),genotype);
}

inline double 
continuous_member_calculator::get_expected_mean(
	size_t         abs_mem_ref,
	genotype_index genotype) const
{
  //lint -e{732}
  return expected_means[genotype][abs_mem_ref];
}

//=====================================================================================
//
//  get_expected_mean(...) POLYGENIC
//
//=====================================================================================

inline double
continuous_member_calculator::get_expected_mean(const member_type& mem,
          	                               genotype_index         genotype,
                                               size_t                 polygenotype) const
{
  return get_expected_mean(get_abs_mem_ref(mem),genotype,polygenotype);
}

inline double 
continuous_member_calculator::get_expected_mean(size_t         abs_mem_ref,
     	                                       genotype_index genotype,
					       size_t         polygenotype) const
{
  //lint -e{732}
  return expected_polygenic_means[polygenotype][genotype][abs_mem_ref];
}

//=====================================================================================
//
//  get_expected_variance(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_expected_variance(
        const member_type& mem,
	genotype_index         genotype) const
{       
  return get_expected_variance(get_abs_mem_ref(mem),genotype);
} 

inline double
continuous_member_calculator::get_expected_variance(
        size_t         abs_mem_ref,
	genotype_index genotype) const
{
  //lint -e{732}
  return expected_variances[genotype][abs_mem_ref];
}

//=====================================================================================
//
//  get_expected_sd(...)
//
//=====================================================================================

inline double 
continuous_member_calculator::get_expected_sd(
        const member_type& mem,
	genotype_index         genotype) const
{
  return get_expected_sd(get_abs_mem_ref(mem),genotype);
}

inline double 
continuous_member_calculator::get_expected_sd(
        size_t                 abs_mem_ref,
	genotype_index         genotype) const
{ 
  //lint -e{732}
  return expected_sds[genotype][abs_mem_ref]; 
}

//=====================================================================================
//
//  get_standardization(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_standardization(
        const member_type& mem,
        genotype_index         genotype) const
{
  return get_standardization(get_abs_mem_ref(mem),genotype);
}

inline double
continuous_member_calculator::get_standardization(
        size_t         abs_mem_ref,
        genotype_index genotype) const
{
  //lint -e{732}
  return standardizations[genotype][abs_mem_ref]; 
}

//=====================================================================================
//
//  get_ascertained_standardization(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_ascertained_standardization(
        const member_type& mem,
        genotype_index         genotype) const
{
  return get_ascertained_standardization(get_abs_mem_ref(mem),genotype);
}

inline double
continuous_member_calculator::get_ascertained_standardization(
        size_t         abs_mem_ref,
        genotype_index genotype) const
{
  //lint -e{732}
  return ascertained_standardizations[genotype][abs_mem_ref]; 
}

//=====================================================================================
//
//  get_estimated_standardization(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_estimated_standardization(
        const member_type& mem,
        genotype_index         genotype_mother,
        genotype_index         genotype_father) const
{
  return get_estimated_standardization(get_abs_mem_ref(mem),genotype_mother,genotype_father);
}

inline double
continuous_member_calculator::get_estimated_standardization(
        size_t         abs_mem_ref,
        genotype_index genotype_mother,
        genotype_index genotype_father) const
{
  //lint -e{732}
  return estimated_standardizations[genotype_mother][genotype_father][abs_mem_ref]; 
}

//=====================================================================================
//
//  get_polygenic_standardization(...)
//
//=====================================================================================

inline double
continuous_member_calculator::get_polygenic_standardization(
        const member_type& mem,
        genotype_index         genotype,
        size_t                 polygenotype) const
{
  return get_polygenic_standardization(get_abs_mem_ref(mem),genotype,polygenotype);
}

inline double
continuous_member_calculator::get_polygenic_standardization(
        size_t         abs_mem_ref,
        genotype_index genotype,
        size_t         polygenotype) const
{
  //lint -e{732}
  return polygenic_standardizations[polygenotype][genotype][abs_mem_ref]; 
}

//=====================================================================================
//
//  update()
//
//=====================================================================================

inline int binary_member_calculator::update()
{
  // First perform all the calculations on those that might produce errors

  int err_code = 0;

  err_code = calculate_susceptibilities(); if(err_code) return err_code;
  err_code = calculate_penetrances ();     if(err_code) return err_code;

  if(mod.get_model_class() == model_FPMM)
  {
    err_code = calculate_polygenic_penetrances (); if(err_code) return err_code;
  }

  return segreg_errors::EVAL_OK;
}

//=====================================================================================
//
//  get_penetrance(...)
//
//=====================================================================================

inline double
binary_member_calculator::get_penetrance(
    const member_type& mem,
    genotype_index     genotype) const
{
  return int_get_penetrance(get_abs_mem_ref(mem),genotype);
}

inline double 
binary_member_calculator::int_get_penetrance(
	size_t         abs_mem_ref,
	genotype_index genotype) const
{
  //lint -e{732}
  return penetrances[genotype][abs_mem_ref];
}

inline bool
binary_member_calculator::get_aff_status(const member_type& mem) const
{
  return int_get_aff_status(get_abs_mem_ref(mem));
}

inline bool 
binary_member_calculator::int_get_aff_status(size_t abs_mem_ref) const
{
  //lint -e{56,48,734} <- Problems with vector<bool> again 
  return affections[abs_mem_ref];
}

//=====================================================================================
//
//  get_expected_susceptibility(...)
//
//=====================================================================================

inline double
binary_member_calculator::get_expected_susc(
    const member_type& mem,
    genotype_index     genotype) const
{
  return int_get_expected_susc(get_abs_mem_ref(mem),genotype);
}

inline double 
binary_member_calculator::int_get_expected_susc(
	size_t         abs_mem_ref,
	genotype_index genotype) const
{
  //lint -e{732}
  return expected_susceptibilities[genotype][abs_mem_ref];
}

//=====================================================================================
//
//  get_penetrance(...) POLYGENIC
//
//=====================================================================================

inline double
binary_member_calculator::get_penetrance(
    const member_type& mem,
    genotype_index     genotype,
    size_t             polygenotype) const
{
  return int_get_penetrance(get_abs_mem_ref(mem),genotype,polygenotype);
}

inline double 
binary_member_calculator::int_get_penetrance(
    size_t         abs_mem_ref,
    genotype_index genotype,
    size_t         polygenotype) const
{
  //lint -e{732}
  return polygenic_penetrances[polygenotype][genotype][abs_mem_ref];
}

//===========================================================================
//
//  unconnected individual likelihood in MLM model(...)
//
//  Equation 84 in page 52 - Appendix F
//===========================================================================
inline double
binary_member_calculator::unconnected_likelihood
  (const member_type& mem, genotype_index geno) const
{
  if(!is_member_valid(mem)) return 1.0;

  double      psi  = mod.freq_sub_model.prob(geno);
  double      pen  = get_penetrance         (mem,geno);

  return psi * pen;
}

/// logistic form of penetrance function.

/// Calculates:
///
/// \f[\frac{e^{my}}{1 + e^{m}}\f]
///
/// where
///
/// \f$m\f$ is the mean and
/// \f$y\f$ is the affection status.
///
/// This is used throughout the mlm (ie, binary model) equations. 

inline double
binary_member_calculator::calculate_penetrance(double mean, double affection) const 
{
  return exp(mean*affection) / (1.0 + exp(mean) );
}                   

//=====================================================================================
//=====================================================================================
//
//  onset_member_calculator
//
//=====================================================================================
//=====================================================================================

inline int onset_member_calculator::update()
{
  // If the geometric mean couldn't be set, we know we can't do a
  // transformation on the ages.

  if(!valid_geom_mean) return segreg_errors::MCC_FAILED_TRANSFORM;

  // First perform all the calculations on those that might produce errors

  int err_code = 0;

  err_code = calculate_transf_age_onset  (); if(err_code) return err_code;
  err_code = calculate_transf_age_exam   (); if(err_code) return err_code;

  err_code = calculate_expected_means    (); if(err_code) return err_code;
  err_code = calculate_expected_suscs    (); if(err_code) return err_code;
  err_code = calculate_expected_alphas   (); if(err_code) return err_code;

  return segreg_errors::EVAL_OK;
}


inline
bool onset_member_calculator::get_aff_status
    (const member_type& mem)  const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_aff_status(mem_index);
}

inline
double onset_member_calculator::get_age_onset
    (const member_type& mem)  const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_age_onset(mem_index);
}


inline
double onset_member_calculator::get_age_exam
    (const member_type& mem)  const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_age_exam(mem_index);
}


inline
double onset_member_calculator::get_expected_age_onset 
    (const member_type& mem,
     genotype_index         genotype,
     size_t                 polygenotype)         const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_expected_age_onset(mem_index, genotype, polygenotype);
}


inline
double onset_member_calculator::get_expected_susc
    (const member_type& mem,
     genotype_index         genotype,
     size_t                 polygenotype)         const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_expected_susc(mem_index, genotype, polygenotype);
}


inline
double onset_member_calculator::get_alpha_i
    (const member_type& mem,
     genotype_index         genotype)             const
{
  size_t mem_index = get_abs_mem_ref(mem);

  return get_alpha_i(mem_index, genotype);
}

// ================
// Internal inlines
// ================

inline
bool onset_member_calculator::get_aff_status (size_t mem_index)  const
{
  //lint -e{56,48,734} <- Problems with vector<bool> again 
  return affections[mem_index];
}

inline
double onset_member_calculator::get_age_onset (size_t mem_index)  const
{
  return transf_age_onsets[mem_index];
}


inline
double onset_member_calculator::get_age_exam (size_t mem_index)  const
{
  return transf_age_exams[mem_index];
}

inline
double onset_member_calculator::get_expected_age_onset
    (size_t mem_index,genotype_index genotype,size_t polygenotype) const     
{
  if(mod.ons_sub_model.m_option() == onset_sub_model::m_A)
    //lint -e{732}
    return expected_means[polygenotype][genotype][mem_index];
  else
    //lint -e{732}
    return expected_means[0][genotype][mem_index];
}


inline
double onset_member_calculator::get_expected_susc 
    (size_t mem_index,genotype_index genotype,size_t polygenotype) const     
{
  if(mod.ons_sub_model.m_option() == onset_sub_model::m_S)
    //lint -e{732}
    return expected_susceptibilities[polygenotype][genotype][mem_index];
  else
    //lint -e{732}
    return expected_susceptibilities[0][genotype][mem_index];
}


inline
double onset_member_calculator::get_alpha_i
    (size_t mem_index,genotype_index genotype) const
{
  //lint -e{732}
  return expected_alphas[genotype][mem_index];
}

}}
