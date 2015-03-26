#ifndef SEGREG_TRANS_SUB_MODEL_H
#include "segreg/transmission_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{

// ================================
//    psm_info inlines
// ================================

inline psm_info::prev_cov::prev_cov()
  : my_name  (),
    my_value (QNAN),
    my_index ((size_t) -1)
{ }

inline bool psm_info::prev_cov::operator==(const prev_cov& rhs) const
{
  return my_name  == rhs.my_name &&
         my_value == rhs.my_value;
}

inline bool psm_info::prev_constraint::operator==
    (const prev_constraint& rhs) const
{
  return my_number_affected == rhs.my_number_affected &&
         my_sample_size     == rhs.my_sample_size     &&
         (SAGE::isnan(my_age) || my_age == rhs.my_age)      &&
         my_susc_covs       == rhs.my_susc_covs       &&
         my_mean_covs       == rhs.my_mean_covs       &&
         my_var_covs        == rhs.my_var_covs;
}

inline bool psm_info::prev_estimate::operator==
    (const prev_estimate& rhs) const
{
  return (SAGE::isnan(my_age) || my_age == rhs.my_age) &&
         my_susc_covs       == rhs.my_susc_covs  &&
         my_mean_covs       == rhs.my_mean_covs  &&
         my_var_covs        == rhs.my_var_covs;

}
                      
// ================================
//    psm_builder inlines
// ================================

inline
psm_builder::psm_builder
    (const MeanCovariateSubmodel*            mean_covs,
     const SusceptibilityCovariateSubmodel*  susc_covs,
     const VarianceCovariateSubmodel*        var_covs)
  : my_mean_covs       (mean_covs),
    my_susc_covs       (susc_covs),
    my_var_covs        (var_covs),
    my_mean_cov_values (mean_covs->get_covariate_count(), QNAN),
    my_susc_cov_values (susc_covs->get_covariate_count(), QNAN),
    my_var_cov_values  (var_covs->get_covariate_count(), QNAN),
    my_age             (QNAN),
    my_sample_size     (QNAN),
    my_number_affected (QNAN)
{
    // Strictly speaking, this shouldn't be needed, but we do it anyway.

    clear_prev_data();
}

inline
size_t psm_builder::set_age     (double age)
{
  if(SAGE::isnan(age) || age <= 0) return 1;

  my_age = age;

  return 0;
}

inline
size_t psm_builder::set_sample_size     (double sample_size)
{
  if(!finite(sample_size) || sample_size <= 0) return 1;

  my_sample_size = sample_size;

  return 0;
}

inline
size_t psm_builder::set_number_affected (double number_affected)
{
  if(!finite(number_affected) || number_affected <= 0) return 1;

  my_number_affected = number_affected;

  return 0;
}

inline
void psm_builder::clear_prev_data()
{
    fill(my_mean_cov_values.begin(),
         my_mean_cov_values.end(),
         QNAN);
    fill(my_susc_cov_values.begin(),
         my_susc_cov_values.end(),
         QNAN);
    fill(my_var_cov_values.begin(),
         my_var_cov_values.end(),
         QNAN);

    my_age             = QNAN;
    my_sample_size     = QNAN;
    my_number_affected = QNAN;
}

/// If value valid, replaces lhs with value and returns either valid or
/// duplicate (depending on if lhs was finite previously).  If value invalid
/// returns value_err.
inline
psm_builder::add_cov_ret_type
    psm_builder::modify_value(double& lhs, double value)
{
  if(finite(value))
  {
    bool previously_valid = finite(lhs);
    
    lhs = value;
    
    if(previously_valid) return psm_builder::duplicate;
    else                 return psm_builder::valid;
  }
  else
  {
    return psm_builder::value_err;
  }
}
                            
inline psm_builder::add_cov_ret_type
    psm_builder::add_mean_cov(const string& cov_name, double value)
{
    size_t index = my_mean_covs->covariate_index(cov_name);

    return modify_value(my_mean_cov_values[index], value);
}

inline psm_builder::add_cov_ret_type
    psm_builder::add_susc_cov(const string& cov_name, double value)
{
    size_t index = my_susc_covs->covariate_index(cov_name);

    return modify_value(my_susc_cov_values[index], value);
}

inline psm_builder::add_cov_ret_type
    psm_builder::add_var_cov(const string& cov_name, double value)
{
    size_t index = my_var_covs->covariate_index(cov_name);

    return modify_value(my_var_cov_values[index], value);
}

// ================================
//    prevalence_sub_model inlines
// ================================

inline
prevalence_sub_model::prevalence_sub_model
   (const genotype_specific_mean_sub_model*           means,
    const genotype_specific_susceptibility_sub_model* suscs,
    const genotype_specific_variance_sub_model*       vars,
    const MeanCovariateSubmodel*                   mean_covs,
    const SusceptibilityCovariateSubmodel*         susc_covs,
    const VarianceCovariateSubmodel*               var_covs,
    const genotype_frequency_sub_model*               freqs,
    const finite_polygenic_mixed_model_sub_model*     fpmms,
    const onset_sub_model*                            onsets,
    const MAXFUN::TransformationSubmodel*             transforms,
    bool                                              fpmm_option,
    bool                                              onset_option)
  : my_means      (means),
    my_suscs      (suscs),
    my_vars       (vars),
    my_mean_covs  (mean_covs),
    my_susc_covs  (susc_covs),
    my_var_covs   (var_covs),
    my_freqs      (freqs),
    my_fpmms      (fpmms),
    my_onsets     (onsets),
    my_transforms (transforms),

    my_fpmm_option  (fpmm_option),
    my_onset_option (onset_option),

    my_estimate_data   (),
    my_constraint_data ()

{ }

inline
prevalence_sub_model::prevalence_sub_model
   (const prevalence_sub_model& psm)
  : SegregSubmodel(psm),
    my_means      (psm.my_means),
    my_suscs      (psm.my_suscs),
    my_vars       (psm.my_vars),
    my_mean_covs  (psm.my_mean_covs),
    my_susc_covs  (psm.my_susc_covs),
    my_var_covs   (psm.my_var_covs),
    my_freqs      (psm.my_freqs),
    my_fpmms      (psm.my_fpmms),
    my_onsets     (psm.my_onsets),
    my_transforms (psm.my_transforms),

    my_fpmm_option  (psm.my_fpmm_option),
    my_onset_option (psm.my_onset_option),

    my_estimate_data   (psm.my_estimate_data),
    my_constraint_data (psm.my_constraint_data)

{ }

inline prevalence_sub_model&
prevalence_sub_model::operator=
   (const prevalence_sub_model& psm)
{
    if(this == &psm) return *this;

    SegregSubmodel::operator=(psm);

    my_means      = psm.my_means;
    my_suscs      = psm.my_suscs;
    my_vars       = psm.my_vars;
    my_mean_covs  = psm.my_mean_covs;
    my_susc_covs  = psm.my_susc_covs;
    my_var_covs   = psm.my_var_covs;
    my_freqs      = psm.my_freqs;
    my_fpmms      = psm.my_fpmms;
    my_onsets     = psm.my_onsets;
    my_transforms = psm.my_transforms;

    my_fpmm_option  = psm.my_fpmm_option;
    my_onset_option = psm.my_onset_option;

    my_estimate_data   = psm.my_estimate_data;
    my_constraint_data = psm.my_constraint_data;

    return *this;
}

inline
prevalence_sub_model::~prevalence_sub_model()
{ }

inline string
prevalence_sub_model::name() const
{
  return "";
}
inline string
prevalence_sub_model::option_description() const
{
  return "";
}

inline bool prevalence_sub_model::set_fpmm_option(bool using_fpmm)
{
  // only fail when onset is true and fpmm is to be set off, because that's
  // not a valid combination.
 
  if(my_onset_option && !using_fpmm) return false;

  my_fpmm_option = using_fpmm;

  return true;
}

inline bool prevalence_sub_model::set_onset_option(bool using_age_onset)
{
  // If the new value is the same as the old one, we do nothing

  if(using_age_onset == my_onset_option) return true;

  // If we're turning onset on, verify that fpmm is already on.

  if(using_age_onset && !my_fpmm_option) return false;

  // Regardless of onset state, the data must be empty to switch.

  if(!my_estimate_data.empty() || ! my_constraint_data.empty())
    return false;

  my_onset_option = using_age_onset;

  return true;
}


inline 
prevalence_sub_model::add_elt_ret_type
    prevalence_sub_model::add_constraint(psm_builder& builder)
{
  add_elt_ret_type ret = add_constraint_internal(builder);
  
  builder.clear_prev_data();
  
  return ret;
}

inline 
prevalence_sub_model::add_elt_ret_type
    prevalence_sub_model::add_estimate(psm_builder& builder)
{
  add_elt_ret_type ret = add_estimate_internal(builder);
  
  builder.clear_prev_data();
  
  return ret;
}

inline
size_t prevalence_sub_model::get_estimate_count() const
{
  return my_estimate_data.size();
}

inline
double prevalence_sub_model::get_estimate_age(size_t i) const
{
  return my_estimate_data[i].my_age;
}

inline
size_t prevalence_sub_model::get_estimate_susc_covariate_count(size_t i) const
{
  return my_estimate_data[i].my_susc_covs.size();
}

inline
const string& prevalence_sub_model::get_estimate_susc_covariate_name
    (size_t est, size_t cov) const
{
  return my_estimate_data[est].my_susc_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_estimate_susc_covariate_value(size_t est, size_t cov) const
{
  return my_estimate_data[est].my_susc_covs[cov].my_value;
}

inline
size_t prevalence_sub_model::get_estimate_mean_covariate_count(size_t i) const
{
  return my_estimate_data[i].my_mean_covs.size();
}

inline
const string& prevalence_sub_model::get_estimate_mean_covariate_name
    (size_t est, size_t cov) const
{
  return my_estimate_data[est].my_mean_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_estimate_mean_covariate_value(size_t est, size_t cov) const
{
  return my_estimate_data[est].my_mean_covs[cov].my_value;
}

inline
size_t prevalence_sub_model::get_estimate_var_covariate_count(size_t i) const
{
  return my_estimate_data[i].my_var_covs.size();
}

inline
const string& prevalence_sub_model::get_estimate_var_covariate_name
    (size_t est, size_t cov) const
{
  return my_estimate_data[est].my_var_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_estimate_var_covariate_value(size_t est, size_t cov) const
{
  return my_estimate_data[est].my_var_covs[cov].my_value;
}

inline
size_t prevalence_sub_model::get_constraint_count() const
{
  return my_constraint_data.size();
}

inline
double prevalence_sub_model::get_constraint_age(size_t i) const
{
  return my_constraint_data[i].my_age;
}

inline
double prevalence_sub_model::get_constraint_number_affected(size_t i) const
{
  return my_constraint_data[i].my_number_affected;  
}

inline
double prevalence_sub_model::get_constraint_sample_size    (size_t i) const
{
  return my_constraint_data[i].my_sample_size;  
}

inline
size_t prevalence_sub_model::get_constraint_susc_covariate_count(size_t i) const
{
  return my_constraint_data[i].my_susc_covs.size();
}

inline
const string& prevalence_sub_model::get_constraint_susc_covariate_name
    (size_t est, size_t cov) const
{
  return my_constraint_data[est].my_susc_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_constraint_susc_covariate_value(size_t est, size_t cov) const
{
  return my_constraint_data[est].my_susc_covs[cov].my_value;
}

inline
size_t prevalence_sub_model::get_constraint_mean_covariate_count(size_t i) const
{
  return my_constraint_data[i].my_mean_covs.size();
}

inline
const string& prevalence_sub_model::get_constraint_mean_covariate_name
    (size_t est, size_t cov) const
{
  return my_constraint_data[est].my_mean_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_constraint_mean_covariate_value(size_t est, size_t cov) const
{
  return my_constraint_data[est].my_mean_covs[cov].my_value;
}

inline
size_t prevalence_sub_model::get_constraint_var_covariate_count(size_t i) const
{
  return my_constraint_data[i].my_var_covs.size();
}

inline
const string& prevalence_sub_model::get_constraint_var_covariate_name
    (size_t est, size_t cov) const
{
  return my_constraint_data[est].my_var_covs[cov].my_name;
}

inline
double prevalence_sub_model::get_constraint_var_covariate_value(size_t est, size_t cov) const
{
  return my_constraint_data[est].my_var_covs[cov].my_value;
}

// Internal inlines

inline double
  prevalence_sub_model::calculate_susc_adj(const psm_info::prev_cov_vector& covs) const
{
  return calculate_adj(covs, my_susc_covs);
}

inline double
  prevalence_sub_model::calculate_mean_adj(const psm_info::prev_cov_vector& covs) const
{
  return calculate_adj(covs, my_mean_covs);
}

inline double
  prevalence_sub_model::calculate_var_adj(const psm_info::prev_cov_vector& covs) const
{
  return calculate_adj(covs, my_var_covs);
}

inline double
prevalence_sub_model::calculate_susc
    (genotype_index u, const psm_info::prev_cov_vector& covs) const
{
  double susc = my_suscs->parameter(u);

  double adj  = calculate_susc_adj(covs);

  return susc + adj;
}

inline double
prevalence_sub_model::calculate_mean
    (genotype_index u, const psm_info::prev_cov_vector& covs) const
{
  double mean = my_means->parameter(u);

  double adj  = calculate_mean_adj(covs);

  return mean + adj;
}

inline double
prevalence_sub_model::calculate_var
    (genotype_index u, const psm_info::prev_cov_vector& covs) const
{
  double var = my_vars->parameter(u);

  double adj  = calculate_var_adj(covs);

  return var + adj;
}

inline double
prevalence_sub_model::calculate_alpha
    (genotype_index u, const psm_info::prev_cov_vector& var_covs) const
{
  double var = calculate_var(u, var_covs);

  double alpha = M_PI / ( sqrt(3.0 * var));

  return alpha;
}

inline double
prevalence_sub_model::calculate_age
    (const psm_info::prev_estimate& est) const
{
  double age = est.my_age;

  bool valid = my_transforms->transform(age);

  if(!valid) return QNAN;
  else       return age;
}

}
}
