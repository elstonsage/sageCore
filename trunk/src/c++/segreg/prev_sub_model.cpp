#include "segreg/prev_sub_model.h"

namespace SAGE
{
namespace SEGREG
{

psm_builder::add_cov_ret_type
    psm_builder::add_covariate(const string& cov_name, double value)
{
    string covariate = toUpper(cov_name);

    // There are three ways for a covariate to be valid.  It must be
    //   1.  in the mean covariates, 
    //   2.  in the susceptibility covariates, or
    //   3.  in the variance covariates.
    //
    // Note that these are mutually exclusive.  We favor mean and susc over
    // var.
    
    // Case 1: Mean Covariates

    if(my_mean_covs->has_covariate(covariate))
    {
      return add_mean_cov(covariate, value);
    }

    // Case 2: Susceptibility Covariates
    
    else if(my_susc_covs->has_covariate(covariate))
    {
      return add_susc_cov(covariate, value);
    }

    // Case 3: Variance Covariates
    
    else if(my_var_covs->has_covariate(covariate))
    {
      return add_var_cov(covariate, value);
    }

    // Otherwise invalid
    
    else
    {
      return unknown;
    }
}

log_double prevalence_sub_model::get_prevalence_penalty() const
{
  // Return the product of the constraints
  
  log_double prev_penalty_product(1.0);
  
  for(size_t i = 0; i < get_constraint_count(); ++i)
  {
    double Prev = calculate_prev(my_constraint_data[i]);

    // Now we must calculate based on n and r as shown in Appendix D of the
    // Segreg Forumla Doc.

    double N = my_constraint_data[i].my_sample_size;
    double R = my_constraint_data[i].my_number_affected;

    double single_penalty = R * log(Prev) + (N - R) * log(1.0 - Prev);

    log_double sp(exp(1.0));

    sp = sp.pow(single_penalty);

    // Store in Product

    prev_penalty_product *= sp;
  }

  return prev_penalty_product;
}

void prevalence_sub_model::copy_prev_covs
    (psm_info::prev_estimate& t, const psm_builder& builder)
{
  // Clear out the vectors, just in case
  
  t.my_susc_covs.clear();
  t.my_mean_covs.clear();
  t.my_var_covs.clear();
  
  // For each finite susc covariate in the psm_builder, add it to the
  // covariate vector

  for(size_t i = 0; i != builder.my_susc_cov_values.size(); ++i)
  {
    if(finite(builder.my_susc_cov_values[i]))
    {
      psm_info::prev_cov pc;
      
      pc.my_name  = my_susc_covs->covariates()[i].trait_name;
      pc.my_value = builder.my_susc_cov_values[i];
      pc.my_index = i;
      
      t.my_susc_covs.push_back(pc);
    }
  }
  
  // For each finite mean covariate in the psm_builder, add it to the
  // covariate vector

  for(size_t i = 0; i != builder.my_mean_cov_values.size(); ++i)
  {
    if(finite(builder.my_mean_cov_values[i]))
    {
      psm_info::prev_cov pc;
      
      pc.my_name  = my_mean_covs->covariates()[i].trait_name;
      pc.my_value = builder.my_mean_cov_values[i];
      pc.my_index = i;
      
      t.my_mean_covs.push_back(pc);
    }
  }
  
  // For each finite var covariate in the psm_builder, add it to the
  // covariate vector

  for(size_t i = 0; i != builder.my_var_cov_values.size(); ++i)
  {
    if(finite(builder.my_var_cov_values[i]))
    {
      psm_info::prev_cov pc;
      
      pc.my_name  = my_var_covs->covariates()[i].trait_name;
      pc.my_value = builder.my_var_cov_values[i];
      pc.my_index = i;
      
      t.my_var_covs.push_back(pc);
    }
  }
}


prevalence_sub_model::add_elt_ret_type
    prevalence_sub_model::add_constraint_internal(const psm_builder& builder)
{
  // Check for simple errors

  if(!finite(builder.my_number_affected))
    return na_not_finite;

  if(!finite(builder.my_sample_size))
    return ss_not_finite;

  if(builder.my_number_affected >= builder.my_sample_size)
    return na_gt_ss;

  if(my_onset_option)
  {
    if(SAGE::isnan(builder.my_age))  return age_required;
  }
  else
  {
    if(!SAGE::isnan(builder.my_age)) return age_present;
  }

  // Create the new constraint

  psm_info::prev_constraint new_constraint;

  new_constraint.my_age             = builder.my_age;

  new_constraint.my_sample_size     = builder.my_sample_size;
  new_constraint.my_number_affected = builder.my_number_affected;

  copy_prev_covs(new_constraint, builder);

  // Check for constraints that are identical to this one.  A constraint is
  // identical if the sample size, number affected, covariates and covariate
  // values are all identical

  prev_constraint_vector::const_iterator i = find(my_constraint_data.begin(),
                                                  my_constraint_data.end(),
                                                  new_constraint);
  
  if(i != my_constraint_data.end())
    return duplicate;
  
  // Add the new constraint to our constraint set

  my_constraint_data.push_back(new_constraint);

  return valid;
}

prevalence_sub_model::add_elt_ret_type
    prevalence_sub_model::add_estimate_internal(const psm_builder& builder)
{
  // Check for simple errors

  if(my_onset_option)
  {
    if(SAGE::isnan(builder.my_age))  return age_required;
  }
  else
  {
    if(!SAGE::isnan(builder.my_age)) return age_present;
  }

  // Create the new estimate

  psm_info::prev_estimate new_estimate;

  new_estimate.my_age             = builder.my_age;

  copy_prev_covs(new_estimate, builder);

  // Check for estimates that are identical to this one.  A estimate is
  // identical if the sample size, number affected, covariates and covariate
  // values are all identical

  prev_estimate_vector::const_iterator i = find(my_estimate_data.begin(),
                                                my_estimate_data.end(),
                                                new_estimate);
  
  if(i != my_estimate_data.end())
    return duplicate;
  
  // Add the new estimate to our estimate set

  my_estimate_data.push_back(new_estimate);

  // Estimates create maxfun parameters.  We need to handle this now that we
  // know the estimate is good.

  MAXFUN::ParameterInput new_est_param
    ("PREVALENCES",
     "prev. est",
     MAXFUN::Parameter::DEPENDENT,
     0,
     NEGATIVE_INF,
     POSITIVE_INF);

  my_parameters.push_back(new_est_param);

  return valid;
}

int prevalence_sub_model::update()
{
  for(size_t i = 0; i < get_estimate_count(); ++i)
  {
    getParam(i) = calculate_prev(my_estimate_data[i]);
  }

  return 0;
}

double prevalence_sub_model::calculate_adj
    (const psm_info::prev_cov_vector& covs,
     const CovariateSubmodel*     csm) const
{
  double adjustment = 0.0;

  psm_info::prev_cov_vector::const_iterator cov = covs.begin();

  for( ; cov != covs.end(); ++cov)
  {
    double mean  = csm->get_covariate_mean(cov->my_index);
    double value = cov->my_value;
    double coeff = csm->covariates()[cov->my_index].coefficient;

    adjustment += coeff * (value - mean);
    continue;
  }
     
  return adjustment;
}

double 
prevalence_sub_model::penetrance(const psm_info::prev_estimate& pr, genotype_index u) const
{
  // Since there is no polygenic component, we skip all that silliness (see
  // below) and simply calculate the susceptibility

  double susc = calculate_susc(u, pr.my_susc_covs);

  double pen = exp(susc) / (1.0 + exp(susc));

  return pen;
}

double 
prevalence_sub_model::penetrance 
    (const psm_info::prev_estimate& pr, genotype_index u, size_t v) const
{
  double susc = calculate_susc(u, pr.my_susc_covs);

  // Determine if the polygenic component is for susceptibility.
  // It is for susceptibility if we're working on a binary trait (not
  // onset) or the susceptibility is multi-dependent

  if(!my_onset_option || my_onsets->m_option() == onset_sub_model::m_S)
  {
    susc += my_fpmms->mean(v);
  }

  double pen = exp(susc) / (1.0 + exp(susc));

  // If we're onsetting

  if(my_onset_option)
  {
    double alpha = calculate_alpha(u, pr.my_var_covs);

    double age = calculate_age(pr);

    double mean = calculate_mean(u, pr.my_mean_covs);

    // Determine if the polygenic component is for mean.
    // It is for mean if the mean is multi-dependent

    if(my_onsets->m_option() == onset_sub_model::m_A)
    {
      mean += my_fpmms->mean(v);
    }

    // Transform the mean

    bool valid = my_transforms->transform(mean);

    if(!valid) return QNAN;
    
    // Calculate the age exponent (this and following from the b1, b4, etc
    // equations in the SEGREG Formula Document, p 30 in the jan 31, 2004
    // edition
    
    double age_exp = alpha * (age - mean);
    
    // Calculate the age penetrance (logistic function of the age_exponent)

    double age_pen = exp(age_exp) / (1.0 + exp(age_exp));

    // Adjust the penetrance based on the above
    
    pen *= age_pen;
  }

  return pen;
}

double prevalence_sub_model::calculate_prev
    (const psm_info::prev_estimate& est) const
{
  double Prev = 0.0;

  // To calculate P_rev, we iterate over main types, and polygenotypes (if
  // FPMM) and calculate the frequency for each.

  for(genotype_index mtype = index_AA; mtype != index_INVALID; ++mtype)
  {
    double pop_frequency = my_freqs->prob(mtype);

    if(!my_fpmm_option)
    {
      double pen = penetrance(est, mtype);
      
      Prev += pop_frequency * pen;
    }
    else
    {
      // We must iterate over polygenotypes
      for(size_t poly = 0; poly != my_fpmms->max_pgt(); ++poly)
      {
        double poly_frequency = my_fpmms->pop_freq(poly);
        
        double pen = penetrance(est, mtype, poly);
        
        Prev += poly_frequency * pop_frequency * pen;
      }
    }
  }

  return Prev;
}



}
}
