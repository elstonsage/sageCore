//============================================================================
// File:      linkage_tests.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   12/19/2 - created.                         djb
//                                                                          
// Notes:     Implementation of linkage test classes.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/linkage_tests.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  non_ss_lod_ratio_test
//============================================================================
//
void
non_ss_lod_ratio_test::calculate()
{
  do_task_calculations<non_ss_lod_ratio_test, non_ss_lod_ratio_result>(*this);
}

// - Alternate hypothesis.
//
void
non_ss_lod_ratio_test::calculate_alt(size_t trait_index, size_t marker_index, non_ss_lod_ratio_result& result)
{
  mle_sub_model  alt_mle(false, false);
  alt_mle.set_strict_limits();
  
  result.alt_ln_like = QNAN;
  
  if(likelihood_finite(alt_mle, trait_index, marker_index, "lod ratio test"))
  {
    MAXFUN::Results  data = maximize(alt_mle, trait_index, marker_index);
    result.alt_ln_like = max_checked_value(data);
  
    if(! SAGE::isnan(result.alt_ln_like))
    {
      result.restricted_alt_theta = alt_mle.average_theta();
    
      result.set_var_cov(data);
    
      result.set_perform_test(data);
      calculate_relaxed_theta(trait_index, marker_index, result, est_at_point_five(data));  
    }
  }
  
  if(SAGE::isnan(result.alt_ln_like))
  {
    result.perform_test = false;
  }
  
  result.restricted_alt_theta_ub = alt_mle.average_theta_ub();
}

// - Calculate likelihood value at complement to estimate of theta in the range
//   [0, 5] and maximize in the range [0, 1] using complement as a starting point
//   if the resulting likelihood is greater than that for the estimate of theta 
//   or the estimate of theta equals .5.
//
void  
non_ss_lod_ratio_test::calculate_relaxed_theta(size_t trait_index, size_t marker_index, 
                                               non_ss_lod_ratio_result& result, bool point_five)
{
  mle_sub_model  relaxed_mle(false, false);
  relaxed_mle.set_relaxed_limits();
  relaxed_mle.set_average_theta(1 - result.restricted_alt_theta);
  mped_calculator  relaxed_mp_calc(my_mped, relaxed_mle, trait_index, marker_index);
  
  double  relaxed_ln_like = relaxed_mp_calc.likelihood().get_log();
  
  if(relaxed_ln_like > result.alt_ln_like || point_five)
  {
    if(likelihood_finite(relaxed_mle, trait_index, marker_index, "lod ratio test"))
    {
      MAXFUN::Results  data = maximize(relaxed_mle, trait_index, marker_index);
      if(! SAGE::isnan(max_checked_value(data)))
      {
        result.relaxed_alt_theta = relaxed_mle.average_theta();
      
        result.set_var_cov_relaxed(data);
      }
    }
  }
}

// - Is theta fixed at .5 ?
//
bool
non_ss_lod_ratio_test::est_at_point_five(const MAXFUN::Results& data)
{
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  const MAXFUN::Parameter      parameter = parameter_manager.getParameter(AVERAGE);

  double  theta_hat = parameter.getFinalEstimate();
  MAXFUN::Parameter::ParamTypeEnum  status = parameter.getFinalType();
  
  return  (status == MAXFUN::Parameter::IND_FIXED_AT_BOUND   ||
           status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND   ) &&
           (.5 - theta_hat) < BOUND_EPS;
}                     

// - Null hypothesis.
//
void
non_ss_lod_ratio_test::calculate_null(size_t trait_index, size_t marker_index, non_ss_lod_ratio_result& result)
{
  mle_sub_model  null_mle(false, false);
  null_mle.set_average_theta(NULL_THETA);
  mped_calculator  null_mp_calc(my_mped, null_mle, trait_index, marker_index);
  
  result.null_ln_like = null_mp_calc.likelihood().get_log();
}

void
non_ss_lod_ratio_test::write_summary(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << non_ss_lod_ratio_result::meta
                  << non_ss_lod_ratio_result::columns;
              
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
non_ss_lod_ratio_test::write_detail(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << non_ss_lod_ratio_result::vc_meta << endl; 
  
  ios::fmtflags old_flags = out.flags();
  string  label_small = ESTIMATES_IN + SMALL_RANGE;
  string  label_large = ESTIMATES_IN + LARGE_RANGE;
  const size_t  label_large_lw = label_large.size();
  
  const size_t  label_small_lw = label_small.size();  
  const size_t  label_small_sw = 10;
  const size_t  label_small_tw = label_small_lw + label_small_sw;

  out << left;
  
  out << setw(_LOCUS.tw()) << _LOCUS[0] 
      << setw(label_small_tw) << label_small
      << setw(label_large_lw) << label_large << endl;
      
  out << setfill(U_CHR) << setw(_LOCUS.lw())    << "" 
      << setfill(' ')   << setw(_LOCUS.sw())    << ""
      << setfill(U_CHR) << setw(label_small_lw) << ""
      << setfill(' ')   << setw(label_small_sw) << ""
      << setfill(U_CHR) << setw(label_large_lw) << "" << endl;

  out << setfill(' ');
  out.flags(old_flags);
    
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_vc_matrix(out);
  }
}

//============================================================================
// IMPLEMENTATION:  ss_lod_ratio_test
//============================================================================
//
void
ss_lod_ratio_test::calculate()
{
  do_task_calculations<ss_lod_ratio_test, ss_lod_ratio_result>(*this);
}

// - Alternate hypothesis.
//
void
ss_lod_ratio_test::calculate_alt(size_t trait_index, size_t marker_index, ss_lod_ratio_result& result)
{
  mle_sub_model  alt_mle(true, false);
  alt_mle.set_strict_limits();
  
  result.alt_ln_like = QNAN;
  
  if(likelihood_finite(alt_mle, trait_index, marker_index, "lod ratio test"))
  {
    MAXFUN::Results  data = maximize(alt_mle, trait_index, marker_index);
    result.alt_ln_like = max_checked_value(data);
  
    if(! SAGE::isnan(result.alt_ln_like))
    {
      result.restricted_alt_thetas.male_theta = alt_mle.male_theta();
      result.restricted_alt_thetas.female_theta = alt_mle.female_theta();
    
      result.set_var_cov(data);
            
      result.set_perform_test(data);  
      calculate_relaxed_theta(trait_index, marker_index, result, 
                          male_est_at_point_five(data), female_est_at_point_five(data));
    }
  }
  
  if(SAGE::isnan(result.alt_ln_like))
  {
    result.perform_test = false;
  }  
  
  result.restricted_alt_theta_ubs.male_theta = alt_mle.male_theta_ub();
  result.restricted_alt_theta_ubs.female_theta = alt_mle.female_theta_ub();
}

// - If one or more recombinations fixed at .5 by Maxfun, do a maximization in [0, 1].
//   If not see if likelihoods of complements have a larger likelihood and if so do a maximization
//   in [0, 1].
//
void  
ss_lod_ratio_test::calculate_relaxed_theta(size_t trait_index, size_t marker_index, 
                        ss_lod_ratio_result& result, bool male_point_five, bool female_point_five)
{
  if(male_point_five || female_point_five)
  {
    calculate_relaxed_theta_bound(trait_index, marker_index, result);
  }
  else
  {
    calculate_relaxed_theta_compl(trait_index, marker_index, result);
  }
}

// - per rce on 4-16-3, use .5, .5 as starting values for maximization in [0, 1]
//   in this case.
//
void  
ss_lod_ratio_test::calculate_relaxed_theta_bound(size_t trait_index, size_t marker_index, 
                                                  ss_lod_ratio_result& result)
{
  mle_sub_model  relaxed_mle(true, false);
  relaxed_mle.set_relaxed_limits();
  relaxed_mle.set_male_theta(.5);
  relaxed_mle.set_female_theta(.5);
    
  if(likelihood_finite(relaxed_mle, trait_index, marker_index, "lod ratio test"))
  {
    MAXFUN::Results  data = maximize(relaxed_mle, trait_index, marker_index);   
    if(! SAGE::isnan(max_checked_value(data)))
    {
      result.relaxed_alt_thetas.male_theta = relaxed_mle.male_theta();
      result.relaxed_alt_thetas.female_theta = relaxed_mle.female_theta();
      result.set_var_cov_relaxed(data);        
    }
  }
}                        

// - Calculate likelihood values at the complements to estimate of thetas in the range
//   [0, 5] and maximize in the range [0, 1] using the complement as a starting point
//   if the resulting likelihood is greater than that for the estimate of thetas.
//
void  
ss_lod_ratio_test::calculate_relaxed_theta_compl(size_t trait_index, size_t marker_index, 
                                           ss_lod_ratio_result& result)
{
  mle_sub_model  relaxed_mle(true, false);
  relaxed_mle.set_relaxed_limits();
  mped_calculator  relaxed_mp_calc(my_mped, relaxed_mle, trait_index, marker_index);

  double  max_ln_like = result.alt_ln_like;
  double  max_male_theta = result.restricted_alt_thetas.male_theta;
  double  max_female_theta = result.restricted_alt_thetas.female_theta;
  
  // - Likelihood of 1st complement.
  //
  relaxed_mle.set_male_theta(result.restricted_alt_thetas.male_theta);
  relaxed_mle.set_female_theta(1 - result.restricted_alt_thetas.female_theta);
  double  relaxed_ln_like1 = relaxed_mp_calc.likelihood().get_log();

  if(relaxed_ln_like1 > max_ln_like)
  {
    max_ln_like = relaxed_ln_like1;
    max_male_theta = result.restricted_alt_thetas.male_theta;
    max_female_theta = 1 - result.restricted_alt_thetas.female_theta;
  }
  
  // - Likelihood of 2nd complement.
  //
  relaxed_mle.set_male_theta(1 - result.restricted_alt_thetas.male_theta);
  relaxed_mle.set_female_theta(result.restricted_alt_thetas.female_theta);
  double  relaxed_ln_like2 = relaxed_mp_calc.likelihood().get_log();
  
  if(relaxed_ln_like2 > max_ln_like)
  {
    max_ln_like = relaxed_ln_like2;
    max_male_theta = 1 - result.restricted_alt_thetas.male_theta;
    max_female_theta = result.restricted_alt_thetas.female_theta;
  }
  
  // - Likelihood of 3rd complement.
  //
  relaxed_mle.set_male_theta(1 - result.restricted_alt_thetas.male_theta);
  relaxed_mle.set_female_theta(1 - result.restricted_alt_thetas.female_theta);
  double  relaxed_ln_like3 = relaxed_mp_calc.likelihood().get_log();
  
  if(relaxed_ln_like3 > max_ln_like)
  {
    max_ln_like = relaxed_ln_like3;
    max_male_theta = 1 - result.restricted_alt_thetas.male_theta;
    max_female_theta = 1 - result.restricted_alt_thetas.female_theta;
  }
  
  // - Maximization.
  //
  if(max_ln_like > result.alt_ln_like)
  {
    relaxed_mle.set_male_theta(max_male_theta);
    relaxed_mle.set_female_theta(max_female_theta);
    
    if(likelihood_finite(relaxed_mle, trait_index, marker_index, "lod ratio test"))
    {
      MAXFUN::Results  data = maximize(relaxed_mle, trait_index, marker_index);
      if(! SAGE::isnan(max_checked_value(data)))
      { 
        result.relaxed_alt_thetas.male_theta = max_male_theta;
        result.relaxed_alt_thetas.female_theta = max_female_theta;
      
        result.set_var_cov_relaxed(data);          
      }
    }
  }
}

// - Is male theta fixed at .5 ?
//
bool
ss_lod_ratio_test::male_est_at_point_five(const MAXFUN::Results& data)
{
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  const MAXFUN::Parameter      parameter = parameter_manager.getParameter(MALE);

  double  theta_hat = parameter.getFinalEstimate();
  MAXFUN::Parameter::ParamTypeEnum  status = parameter.getFinalType();
  
  return  (status == MAXFUN::Parameter::IND_FIXED_AT_BOUND   ||
           status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND   ) &&
           (.5 - theta_hat) < BOUND_EPS;  
}

// - Is female theta fixed at .5 ?
//
bool
ss_lod_ratio_test::female_est_at_point_five(const MAXFUN::Results& data)
{
  const MAXFUN::ParameterMgr&  parameter_manager = data.getParameterMgr();
  const MAXFUN::Parameter      parameter = parameter_manager.getParameter(FEMALE);

  double  theta_hat = parameter.getFinalEstimate();
  MAXFUN::Parameter::ParamTypeEnum  status = parameter.getFinalType();
  
  return  (status == MAXFUN::Parameter::IND_FIXED_AT_BOUND   ||
           status == MAXFUN::Parameter::IND_FIXED_NEAR_BOUND   ) &&
           (.5 - theta_hat) < BOUND_EPS;  
}


// - Null hypothesis.
//
void
ss_lod_ratio_test::calculate_null(size_t trait_index, size_t marker_index, ss_lod_ratio_result& result)
{
  mle_sub_model  null_mle(true, false);
  null_mle.set_male_theta(NULL_THETA);
  null_mle.set_female_theta(NULL_THETA);
  mped_calculator  null_mp_calc(my_mped, null_mle, trait_index, marker_index);
  
  result.null_ln_like = null_mp_calc.likelihood().get_log();
}

void
ss_lod_ratio_test::write_summary(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << ss_lod_ratio_result::meta
                  << ss_lod_ratio_result::columns;
              
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
ss_lod_ratio_test::write_detail(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << ss_lod_ratio_result::vc_meta << endl; 
  
  ios::fmtflags old_flags = out.flags();
  string  label_small = ESTIMATES_IN + SMALL_RANGE;
  string  label_large = ESTIMATES_IN + LARGE_RANGE;
  const size_t  label_large_lw = label_large.size();
  
  const size_t  label_small_lw = label_small.size();  
  const size_t  label_small_sw = 20;
  const size_t  label_small_tw = label_small_lw + label_small_sw;

  out << left;
  
  out << setw(_LOCUS.tw()) << _LOCUS[0] 
      << setw(label_small_tw) << label_small
      << setw(label_large_lw) << label_large << endl;
      
  out << setfill(U_CHR) << setw(_LOCUS.lw())    << "" 
      << setfill(' ')   << setw(_LOCUS.sw())    << ""
      << setfill(U_CHR) << setw(label_small_lw) << ""
      << setfill(' ')   << setw(label_small_sw) << ""
      << setfill(U_CHR) << setw(label_large_lw) << "" << endl;

  out << setfill(' ');
  out.flags(old_flags);
    
  for(size_t r = 0; r < my_results.size(); ++r)

  {
    my_results[r]->write_vc_matrix(out);
  }
}

//============================================================================
// IMPLEMENTATION:  cleves_elston_test
//============================================================================
//
void
cleves_elston_test::calculate()
{
  do_task_calculations<cleves_elston_test, cleves_elston_result>(*this);
}

// - Alternate hypothesis.
//
void
cleves_elston_test::calculate_alt(size_t trait_index, size_t marker_index, cleves_elston_result& result)
{
  mle_sub_model  alt_mle(true, false);
  alt_mle.set_relaxed_limits();
  
  result.alt_ln_like = QNAN;
  
  if(likelihood_finite(alt_mle, trait_index, marker_index, "Cleves Elston test"))
  {
    MAXFUN::Results  data = maximize(alt_mle, trait_index, marker_index);
    result.alt_ln_like = max_checked_value(data);
  
    if(! SAGE::isnan(result.alt_ln_like))
    {
      result.alt_thetas.male_theta = alt_mle.male_theta();
      result.alt_thetas.female_theta = alt_mle.female_theta();
      result.set_var_cov(data);    
    }
  }
  
  result.alt_theta_ubs.male_theta = alt_mle.male_theta_ub();
  result.alt_theta_ubs.female_theta = alt_mle.female_theta_ub();
}

// - Null hypothesis.
//
void
cleves_elston_test::calculate_null(size_t trait_index, size_t marker_index, 
                                                cleves_elston_result& result)
{
  mle_sub_model  null_mle(true, false);
  null_mle.set_relaxed_limits();
  null_mle.constrain_thetas();
  
  result.null_ln_like = QNAN;
  
  if(likelihood_finite(null_mle, trait_index, marker_index, "Cleves Elston test"))
  {
    MAXFUN::Results  data = maximize(null_mle, trait_index, marker_index);
    result.null_ln_like = max_checked_value(data);
    if(! SAGE::isnan(result.null_ln_like))
    {
      result.null_thetas.male_theta = null_mle.male_theta();
      result.null_thetas.female_theta = null_mle.female_theta();
    }
  }
  
  result.null_theta_ubs.female_theta = null_mle.female_theta_ub();
  result.null_theta_ubs.male_theta = null_mle.male_theta_ub();
}

void
cleves_elston_test::write_summary(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << cleves_elston_result::meta
              << cleves_elston_result::columns;
              
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
cleves_elston_test::write_detail(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << cleves_elston_result::vc_meta << endl; 
  
  ios::fmtflags old_flags = out.flags();

  out << left;
  out << setw(_LOCUS.tw()) << _LOCUS[0] << endl;
  out << setfill(U_CHR) << setw(_LOCUS.lw()) << "" << endl;

  out.flags(old_flags);
    
  for(size_t r = 0; r < my_results.size(); ++r)

  {
    my_results[r]->write_vc_matrix(out);
  }
}

}
}



