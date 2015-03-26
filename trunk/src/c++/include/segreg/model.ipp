//============================================================================
// File:      model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/15/01 created        -djb
//                                                                          
// Notes:     Inline implementation of struct, model, and class, parser.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_MODEL_H
#include "segreg/model.h"
#endif

namespace SAGE
{

namespace SEGREG
{

inline
const genotype_specific_mean_susc_sub_model&
    model::type_dependent_sub_model() const
{
  switch(get_primary_trait_type())
  {
    case pt_CONTINUOUS :
      return mean_sub_model;

    case pt_BINARY     :
      return susc_sub_model;

    case pt_ONSET      :
      if(ons_sub_model.t_option() == onset_sub_model::t_A)
      {
        return mean_sub_model;
      }
      else
      {
        return susc_sub_model;
      }

    case pt_NONE       :
    default            :
    
      SAGE_internal_error();  // Should never happen!
  }

  SAGE_internal_error(); // Should never happen

  return mean_sub_model;
}

inline
genotype_specific_mean_susc_sub_model& model::type_dependent_sub_model()
{
  switch(get_primary_trait_type())
  {
    case pt_CONTINUOUS :
      return mean_sub_model;

    case pt_BINARY     :
      return susc_sub_model;

    case pt_ONSET      :
      if(ons_sub_model.t_option() == onset_sub_model::t_A)
      {
        return mean_sub_model;
      }
      else
      {
        return susc_sub_model;
      }

    case pt_NONE       :
    default            :
    
      cout << get_primary_trait_type() << endl;

      SAGE_internal_error();  // Should never happen!
  }

  SAGE_internal_error(); // Should never happen

  return mean_sub_model;
}

inline std::string
model_class_2_string(model_class mc)
{
  switch(mc)
  {
    case model_A:
      return "A";
    case model_D:
      return "D";
    case model_FPMM:
      return "mixed";
    case model_MLM:
      return "mlm";
    case model_INVALID:
      return "invalid";
    default:
      return ""; 
  }
}

inline std::string
primary_type_2_string(primary_type pt)
{
  switch(pt)
  {
    case pt_NONE       : return "none";
    case pt_CONTINUOUS : return "continuous";
    case pt_BINARY     : return "binary";
    case pt_ONSET      : return "onset";
    default            : SAGE_internal_error();
  }
  return "";
}

//============================================================================
// IMPLEMENTATION:  model
//============================================================================
//
inline
model::model(cerrorstream& errors)
    : mean_sub_model(),
      var_sub_model(&(this->mean_sub_model), errors),
      susc_sub_model(),
      comp_trait_sub_model(&(this->mean_sub_model), errors),
      mean_cov_sub_model(&(this->mean_sub_model), errors),
      var_cov_sub_model(&(this->mean_sub_model), errors),
      susc_cov_sub_model(&(this->susc_sub_model), errors),
      fpmm_sub_model(),
      ons_sub_model(),
      resid_sub_model(),
      transf_sub_model(),
      freq_sub_model(),
      transm_sub_model(&(this->freq_sub_model), errors),
      ascer_sub_model(),
      prev_sub_model(&mean_sub_model,
                     &susc_sub_model,
                     &var_sub_model,
                     &mean_cov_sub_model,
                     &susc_cov_sub_model,
                     &var_cov_sub_model,
                     &freq_sub_model,
                     &fpmm_sub_model,
                     &ons_sub_model,
                     &transf_sub_model,
                     (m_class == model_FPMM),
                     (primary_trait_type == pt_ONSET))
{
  reset(errors);
}

inline
model::model(const model& other)
   : mean_sub_model(other.mean_sub_model),
     var_sub_model(other.var_sub_model),
     susc_sub_model(other.susc_sub_model),
     comp_trait_sub_model(other.comp_trait_sub_model),
     mean_cov_sub_model(other.mean_cov_sub_model),
     var_cov_sub_model(other.var_cov_sub_model),
     susc_cov_sub_model(other.susc_cov_sub_model),
     fpmm_sub_model(other.fpmm_sub_model),
     ons_sub_model(other.ons_sub_model),
     resid_sub_model(other.resid_sub_model),
     transf_sub_model(other.transf_sub_model),
     freq_sub_model(other.freq_sub_model),
     transm_sub_model(other.transm_sub_model),
     ascer_sub_model(other.ascer_sub_model),
     prev_sub_model(other.prev_sub_model)
{
  file_name_root = other.file_name_root;
  title = other.title;
  m_class = other.m_class;
  each_pedigree = other.each_pedigree;
  pen_func_output = other.pen_func_output;
  type_prob = other.type_prob;
  primary_trait = other.primary_trait;
  primary_trait_type = other.primary_trait_type;
  mean_missing = other.mean_missing;
  susc_missing = other.susc_missing;
  trans_missing = other.trans_missing;
  transm_sub_model.my_gf_ptr = &(this->freq_sub_model);
  var_sub_model.my_m_ptr = &(this->mean_sub_model);

  prev_sub_model.my_means      = &(mean_sub_model);
  prev_sub_model.my_suscs      = &(susc_sub_model);
  prev_sub_model.my_vars       = &(var_sub_model);
  prev_sub_model.my_mean_covs  = &(mean_cov_sub_model);
  prev_sub_model.my_susc_covs  = &(susc_cov_sub_model);
  prev_sub_model.my_var_covs   = &(var_cov_sub_model);
  prev_sub_model.my_freqs      = &(freq_sub_model);
  prev_sub_model.my_fpmms      = &(fpmm_sub_model);
  prev_sub_model.my_onsets     = &(ons_sub_model);
  prev_sub_model.my_transforms = &(transf_sub_model);
  
  my_maxfun_debug            = other.my_maxfun_debug;

  rebuild_covariates();
}

inline model&
model::operator=(const model& other)
{
  if(this != &other)
  {
    mean_sub_model        = other.mean_sub_model;
    susc_sub_model        = other.susc_sub_model;
    var_sub_model         = other.var_sub_model;
    freq_sub_model        = other.freq_sub_model;
    resid_sub_model       = other.resid_sub_model;
    transm_sub_model      = other.transm_sub_model;
    transf_sub_model      = other.transf_sub_model;
    mean_cov_sub_model    = other.mean_cov_sub_model;
    var_cov_sub_model     = other.var_cov_sub_model;
    susc_cov_sub_model    = other.susc_cov_sub_model;
    comp_trait_sub_model  = other.comp_trait_sub_model;
    fpmm_sub_model        = other.fpmm_sub_model;
    ons_sub_model         = other.ons_sub_model;
    ascer_sub_model       = other.ascer_sub_model;
    prev_sub_model        = other.prev_sub_model;
    
    transm_sub_model.my_gf_ptr = &(this->freq_sub_model);
    var_sub_model.my_m_ptr        = &(this->mean_sub_model);

    prev_sub_model.my_means      = &(mean_sub_model);
    prev_sub_model.my_suscs      = &(susc_sub_model);
    prev_sub_model.my_vars       = &(var_sub_model);
    prev_sub_model.my_mean_covs  = &(mean_cov_sub_model);
    prev_sub_model.my_susc_covs  = &(susc_cov_sub_model);
    prev_sub_model.my_var_covs   = &(var_cov_sub_model);
    prev_sub_model.my_freqs      = &(freq_sub_model);
    prev_sub_model.my_fpmms      = &(fpmm_sub_model);
    prev_sub_model.my_onsets     = &(ons_sub_model);
    prev_sub_model.my_transforms = &(transf_sub_model);

    rebuild_covariates();

    file_name_root = other.file_name_root;
    title = other.title;
    m_class = other.m_class;
    each_pedigree = other.each_pedigree;
    pen_func_output = other.pen_func_output;
    type_prob = other.type_prob;
    primary_trait = other.primary_trait;
    mean_missing = other.mean_missing;
    susc_missing = other.susc_missing;
    trans_missing = other.trans_missing;
    primary_trait_type = other.primary_trait_type;
    
    my_maxfun_debug            = other.my_maxfun_debug;
  }
  
  return *this;
}

inline std::string  
model::get_file_name_root() const
{
  return file_name_root;
}

inline std::string  
model::get_title() const
{
  return title;
}

inline model_class  
model::get_model_class() const
{
  return m_class;
}

inline bool         
model::get_each_pedigree() const
{
  return each_pedigree;
}
    
inline bool         
model::get_pen_func_output() const
{
  return pen_func_output;
}

inline bool         
model::get_type_prob() const
{
  return type_prob;
}
    
inline std::string  
model::get_primary_trait() const
{
  return primary_trait;
}

inline primary_type
model::get_primary_trait_type() const
{
  return primary_trait_type;
}

inline bool
model::get_type_missing() const
{
  switch(get_primary_trait_type())
  {
    case pt_CONTINUOUS :
      return mean_missing;

    case pt_BINARY     :
      return susc_missing;

    case pt_ONSET      :
      if(ons_sub_model.t_option() == onset_sub_model::t_A)
      {
        return mean_missing;
      }
      else
      {
        return susc_missing;
      }

    case pt_NONE       :

      // Should never be none if model is invalid, but if it is, we can
      // always think of it as missing.
      if(get_model_class() == model_INVALID)
        return true;

      //lint -fallthrough

    default            :

      SAGE_internal_error();  // Should never happen!

  }

  SAGE_internal_error(); // Should never happen

  return mean_missing;
}

inline bool
model::get_trans_missing() const
{
  return trans_missing;
}

/// Returns \c true if the model likelihood is affected by the sex of
/// members in it, \c false otherwise.  Note that this only includes
/// things where sex changes behaviors in calculating likelihoods.  
/// Sex as a covariate is not relevant to this check.
inline bool
model::has_sex_effect() const
{
  bool has_residual_sex_effect     = resid_sub_model.has_sex_effect();
  bool has_transmission_sex_effect = transm_sub_model.has_sex_effect();;
  
  return has_residual_sex_effect || has_transmission_sex_effect;
}

inline void model::rebuild_covariates()
{
  mean_cov_sub_model.my_m_ptr   = &(this->mean_sub_model);
  var_cov_sub_model.my_m_ptr    = &(this->mean_sub_model);
  susc_cov_sub_model.my_m_ptr   = &(this->susc_sub_model);
  comp_trait_sub_model.my_m_ptr = &(this->mean_sub_model);
}    

}}
