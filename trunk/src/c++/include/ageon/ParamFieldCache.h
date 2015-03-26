#ifndef PARAM_FIELD_CACHE_H
#define PARAM_FIELD_CACHE_H

#include "maxfunapi/maxfunapi.h"
#include "sampling/sampling.h"

namespace SAGE {
namespace AO   {

class ParamFieldCache
{
  public:
  
  /// @name Constructors
  //@{
  
    ParamFieldCache(const MAXFUN::ParameterMgr & parameter_mgr, const SAMPLING::MemberDataSample & analysis_sample);
    ParamFieldCache(const ParamFieldCache &);
    
  //@}

  /// @name Mean intercept & covariates
  //@{
  
    MAXFUN::ParameterConstIterator mean_int;
    
    MAXFUN::ParameterConstIterator mean_cov_param_begin;
    MAXFUN::ParameterConstIterator mean_cov_param_end;

    SAMPLING::FieldConstIterator mean_cov_field_begin;
    SAMPLING::FieldConstIterator mean_cov_field_end;

  //@}  

  /// @name Variance intercept & covariates
  //@{
  
    MAXFUN::ParameterConstIterator var_int;
    
    MAXFUN::ParameterConstIterator var_cov_param_begin;
    MAXFUN::ParameterConstIterator var_cov_param_end;

    SAMPLING::FieldConstIterator var_cov_field_begin;
    SAMPLING::FieldConstIterator var_cov_field_end;

  //@}

  /// @name Susceptibility intercepts & covariates
  //@{
  
    MAXFUN::ParameterConstIterator suscept_int_begin;
    MAXFUN::ParameterConstIterator suscept_int_end;
    
    MAXFUN::ParameterConstIterator suscept_cov_param_begin;
    MAXFUN::ParameterConstIterator suscept_cov_param_end;
    
    SAMPLING::FieldConstIterator suscept_cov_field_begin;
    SAMPLING::FieldConstIterator suscept_cov_field_end;
    
  //@}

  private:
  
    ParamFieldCache& operator=(const ParamFieldCache &);
};

//=================================================
//  INLINES
//=================================================

inline
ParamFieldCache::ParamFieldCache(
  const MAXFUN::ParameterMgr & mgr, 
  const SAMPLING::MemberDataSample & sample)
  :
  mean_int                 (mgr.getParamBegin    ("Mean intercept")),
  mean_cov_param_begin     (mgr.getParamBegin    ("Mean covariates")),
  mean_cov_param_end       (mgr.getParamEnd      ("Mean covariates")),
  mean_cov_field_begin     (sample.getFieldBegin ("Mean covariates")),
  mean_cov_field_end       (sample.getFieldEnd   ("Mean covariates")),
  var_int                  (mgr.getParamBegin    ("Variance intercept")),
  var_cov_param_begin      (mgr.getParamBegin    ("Variance covariates")),
  var_cov_param_end        (mgr.getParamEnd      ("Variance covariates")),
  var_cov_field_begin      (sample.getFieldBegin ("Variance covariates")),
  var_cov_field_end        (sample.getFieldEnd   ("Variance covariates")),
  suscept_int_begin        (mgr.getParamBegin    ("Susceptibility intercepts")),
  suscept_int_end          (mgr.getParamEnd      ("Susceptibility intercepts")),
  suscept_cov_param_begin  (mgr.getParamBegin    ("Susceptibility covariates")),
  suscept_cov_param_end    (mgr.getParamEnd      ("Susceptibility covariates")),
  suscept_cov_field_begin  (sample.getFieldBegin ("Susceptibility covariates")),
  suscept_cov_field_end    (sample.getFieldEnd   ("Susceptibility covariates"))
{ }

inline
ParamFieldCache::ParamFieldCache(const ParamFieldCache & other) :
  mean_int                 (other.mean_int),
  mean_cov_param_begin     (other.mean_cov_param_begin),
  mean_cov_param_end       (other.mean_cov_param_end),
  mean_cov_field_begin     (other.mean_cov_field_begin),
  mean_cov_field_end       (other.mean_cov_field_end),
  var_int                  (other.var_int),
  var_cov_param_begin      (other.var_cov_param_begin),
  var_cov_param_end        (other.var_cov_param_end),
  var_cov_field_begin      (other.var_cov_field_begin),
  var_cov_field_end        (other.var_cov_field_end),
  suscept_int_begin        (other.suscept_int_begin),
  suscept_int_end          (other.suscept_int_end),
  suscept_cov_param_begin  (other.suscept_cov_param_begin),
  suscept_cov_param_end    (other.suscept_cov_param_end),
  suscept_cov_field_begin  (other.suscept_cov_field_begin),
  suscept_cov_field_end    (other.suscept_cov_field_end)
{ }




} // End namespace AO
} // End namespace SAGE

#endif
