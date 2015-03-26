//============================================================================
// File:      tasks.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   12/19/2 - created.                         djb
//                                                                          
// Notes:     Implementation of misc classes.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/tasks.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  task
//============================================================================
//
// - Find the maximum likelihood for all of the data using the given
//   mle_sub_model.
//
MAXFUN::Results
task::maximize(mle_sub_model& mle, size_t trait_index, size_t marker_index) const
{
  mped_calculator  mp_calc(my_mped, mle, trait_index, marker_index);
  MAXFUN::ParameterMgr  param_mgr;
  MAXFUN::DebugCfg      debug_cfg;
  MAXFUN::SequenceCfg   sequence_cfg;
  
  int  failure = param_mgr.addSubModel(&mle);
  assert(! failure);
  
  set_maxfun_configuration(sequence_cfg);
  
  return  MAXFUN::maximize(param_mgr, mp_calc, sequence_cfg, debug_cfg);  
}

// - Find the maximum likelihood for data in the given group using the given
//   mle_sub_model.
//
MAXFUN::Results
task::maximize(const group& g, mle_sub_model& mle, size_t trait_index, size_t marker_index)
{
  group_calculator      group_calc(g, my_mped, mle, trait_index, marker_index);
  MAXFUN::ParameterMgr  param_mgr;
  MAXFUN::DebugCfg      debug_cfg;
  MAXFUN::SequenceCfg   sequence_cfg;
  
  int  failure = param_mgr.addSubModel(&mle);
  assert(! failure);
  
  set_maxfun_configuration(sequence_cfg);
  
  return  MAXFUN::maximize(param_mgr, group_calc, sequence_cfg, debug_cfg);  
}

// - Evaluate likelihood function for entire multipedigree.  If result is finite, 
//   return true; otherwise, write an error message and return false.
//
bool
task::likelihood_finite(mle_sub_model& mle, size_t trait, size_t marker, const string& test)
{
  mped_calculator  mp_calc(my_mped, mle, trait, marker);
  
  if(finite(mp_calc.likelihood().get_log()))
  {
    return  true;
  }
  else
  {
    my_errors << priority(error) << "Likelihood could not be evaluated or had a value "
              << "of zero for marker '" << marker_name(marker) << "' while performing "
              << test << ".  Skipping " << test << " for marker'" << marker_name(marker)
              << "' ..." << endl;
    return  false;
  }
}

// - Evaluate likelihood function for the given group.  If result is finite, 
//   return true; otherwise, write an error message and return false.
//
bool
task::likelihood_finite(const string& group_name, const group& g, mle_sub_model& mle, 
                        size_t trait, size_t marker, const string& test)
{
  group_calculator  group_calc(g, my_mped, mle, trait, marker);
  
  if(finite(group_calc.likelihood().get_log()))
  {
    return  true;
  }
  else
  {
    my_errors << priority(error) << "Likelihood for pedigree '" << group_name 
              << "' could not be evaluated or had a value "
              << "of zero for marker '" << marker_name(marker) << "' while performing "
              << test << ".  Skipping " << test << " for marker'" << marker_name(marker) 
              << "' ..." << endl;
    return  false;
  }
}

string
task::marker_name(size_t marker_index)
{
  const RPED::RefMPedInfo&  mp_info = my_mped.info();
  
  assert(marker_index < mp_info.marker_count());
  return  mp_info.marker_info(marker_index).name();
}


//============================================================================
// IMPLEMENTATION:  non_ss_smiths_faraways_test
//============================================================================
//
bool  non_ss_smiths_faraways_test::alt_calculated = false;
vector<non_ss_alt_result_ptr>  non_ss_smiths_faraways_test::alt_results;

// - Clear static members.
//
void
non_ss_smiths_faraways_test::clear()
{
  alt_calculated = false;
  alt_results.clear();
}

void
non_ss_smiths_faraways_test::calculate()
{
  if(my_instructions.valid)
  {
    const RPED::RefMPedInfo&  mp_info = my_mped.info();
    size_t  trait_index = mp_info.marker_find(my_instructions.trait);
    
    for(size_t marker_index = 0; marker_index < mp_info.marker_count(); ++marker_index)
    {
      if(marker_index == trait_index)
      {
        continue;
      }
      
      non_ss_smiths_faraways_result*  result = 0;
      
      if(my_type == sf_LINKAGE)
      {
        result = new non_ss_faraways_result;
      }
      else
      {
        result = new non_ss_smiths_result;
      }
      
      result->trait = my_instructions.trait;
      result->marker = mp_info.marker_info(marker_index).name();
      
      calculate_alt(trait_index, marker_index, *result);
      calculate_null(trait_index, marker_index, *result);      
      if(my_type == sf_LINKAGE)
      {
        calculate_posteriors(trait_index, marker_index, *result);
      }
      
      my_results.push_back(result_ptr(result)); 
    }

    alt_calculated = true;    
    completed = true;
  }
}


// - Maximize likelihood while estimating the proportion of families w.
//   linkage (allowing for heterogeneity).
//
void
non_ss_smiths_faraways_test::calculate_alt(size_t trait_index, size_t marker_index, 
                                              non_ss_smiths_faraways_result& result)
{
  if(! alt_calculated)    // Keep a static copy of results.
  {
    non_ss_alt_result*  alt_result = new non_ss_alt_result;
    alt_result->marker = result.marker;
    
    mle_sub_model  alt_mle(false, true);
    alt_mle.set_strict_limits();
  
    result.alt_ln_like = QNAN;
    
    if(likelihood_finite(alt_mle, trait_index, marker_index, "likelihood calculation assuming heterogeneity"))
    {
      MAXFUN::Results  data = maximize(alt_mle, trait_index, marker_index);
      result.alt_ln_like = max_checked_value(data);
    
      if(! SAGE::isnan(result.alt_ln_like))
      {
        result.alt_theta = alt_mle.average_theta();
        alt_result->alt_theta = result.alt_theta;
        result.alpha = alt_mle.alpha();
        alt_result->alpha = result.alpha;
        result.set_var_cov(data);
        alt_result->var_cov = result.var_cov;          
      }
    }
  
    alt_result->alt_ln_like = result.alt_ln_like;  
  
    result.alt_theta_ub = alt_mle.average_theta_ub();
    alt_result->alt_theta_ub = result.alt_theta_ub;
    
    alt_results.push_back(non_ss_alt_result_ptr(alt_result));
  }
  else   // Get results from static copy.
  {
    vector<non_ss_alt_result_ptr>::const_iterator  r_iter = find_if(alt_results.begin(), alt_results.end(), 
                                                            bind2nd(has_marker<non_ss_alt_result_ptr>(), result.marker));
    assert(r_iter != alt_results.end());
    
    result.alt_ln_like = (*r_iter)->alt_ln_like;
    result.alt_theta = (*r_iter)->alt_theta;
    result.alpha = (*r_iter)->alpha;
    result.alt_theta_ub = (*r_iter)->alt_theta_ub;
    result.var_cov = (*r_iter)->var_cov;
  }
}

// - Maximize likelihood for the homogeneous case.
//
void
non_ss_smiths_faraways_test::calculate_null(size_t trait_index, size_t marker_index, 
                                                non_ss_smiths_faraways_result& result)
{
  mle_sub_model  null_mle(false, false);
  null_mle.set_strict_limits();
  
  if(my_type == sf_LINKAGE)
  {
    null_mle.set_average_theta(NULL_THETA);
    mped_calculator  null_mp_calc(my_mped, null_mle, trait_index, marker_index);
  
    result.null_ln_like = null_mp_calc.likelihood().get_log();
  }
  else
  {
    result.null_ln_like = QNAN; 
    
    if(likelihood_finite(null_mle, trait_index, marker_index, "Smith's test"))
    {
      MAXFUN::Results  data = maximize(null_mle, trait_index, marker_index);
      result.null_ln_like = max_checked_value(data);
    }
      
    if(! SAGE::isnan(result.null_ln_like))
    {
      result.null_theta = null_mle.average_theta();
    }
  }
  
  result.null_theta_ub = null_mle.average_theta_ub();
}

// - Calculate posterior likelihoods for each pedigree as a whole AND for each
//   subpedigree w/i a given pedigree.  Subpedigrees may be identified
//   by a combination of pedigree id and name of an arbitrary member.
//   All calculations based on the alpha and theta(s) that was estimated for
//   the entire multipedigree.
//                                            -per discussion w. RCE 1/23/3.
//
void
non_ss_smiths_faraways_test::calculate_posteriors(size_t trait_index, size_t marker_index, 
                                                     non_ss_smiths_faraways_result& result)
{
  if(SAGE::isnan(result.alt_theta) || SAGE::isnan(result.alpha)) 
  {
    return;
  }
  
  mle_sub_model  mle(false, false);
  mle.set_average_theta(result.alt_theta);

  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    string  ped_name = ped_iter->name();
    pedigree_posterior  ped_p(ped_name);
    double              ped_alt_ln_like = 0;
    double              ped_null_ln_like = 0;
  
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      peeler  p(*subped_iter, mle, trait_index, marker_index);
      subped_calculator  sp_calc(p);
      double  alt_ln_like  = sp_calc.likelihood().get_log();
      double  null_ln_like = sp_calc.unlinked_likelihood().get_log();
      
      double  term1 = result.alpha * exp(alt_ln_like);
      double  term2 = (1 - result.alpha) * exp(null_ln_like);
      double  posterior = term1 / (term1 + term2);
      
      string  member_name = subped_iter->member_begin()->name();
      
      ped_alt_ln_like += alt_ln_like;
      ped_null_ln_like += null_ln_like;
      ped_p.sub_posteriors.push_back(subpedigree_posterior(member_name, posterior)); 
    }
    
    double  ped_term1 = result.alpha * exp(ped_alt_ln_like);
    double  ped_term2 = (1 - result.alpha) * exp(ped_null_ln_like);
    ped_p.posterior = ped_term1 / (ped_term1 + ped_term2);
      
    result.posteriors.push_back(ped_p);
  }
}

void
non_ss_smiths_faraways_test::write_summary(ostream& out) const
{
  assert(completed);

  if(my_type == sf_HOMOGENEITY)
  {
    out << "\n\n\n" << non_ss_smiths_result::meta
                    << non_ss_smiths_result::columns;
  }
  else 
  {
    out << "\n\n\n" << non_ss_faraways_result::meta
                    << non_ss_faraways_result::columns;  
  }
  
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
non_ss_smiths_faraways_test::write_detail(ostream& out) const
{
  assert(completed);

  if(my_type == sf_LINKAGE)
  {
    ios::fmtflags old_flags = out.flags();
  
    out << setfill(' ');
    out << left;
    
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << POSTERIORS_FAM_ << "\n"  
        << setw(_LOCUS.tw()) << ""
        << EVIDENCE_ << "\n\n" << endl;
  
    pedigree_const_iterator  ped_iter;
    for(ped_iter = my_mped.pedigree_begin(); ped_iter != my_mped.pedigree_end(); ++ped_iter)
    {
      subpedigree_const_iterator  subped_iter;
      for(subped_iter = ped_iter->subpedigree_begin(); 
          subped_iter != ped_iter->subpedigree_end();
          ++subped_iter)
      {
        out << setw(_LOCUS.tw()) << ""
            << "Constituent Pedigree in Pedigree " << ped_iter->name() 
            << " Containing Member " << subped_iter->member_begin()->name()
            << endl;
            
        out << non_ss_faraways_result::detail_columns;
        for(size_t i = 0; i < my_results.size(); ++i)
        {
          non_ss_faraways_result* ptr = static_cast<non_ss_faraways_result*>(my_results[i].get());
          ptr->write_family_detail(out, ped_iter->name(), subped_iter->member_begin()->name());
        }
        
        out << "\n" << endl;
      }
    }    
  
    out.flags(old_flags); 
  }
  
  write_vc_matrix(out);
}

void
non_ss_smiths_faraways_test::write_vc_matrix(ostream& out) const
{
  assert(completed);
  
  if(my_type == sf_LINKAGE)
  {
    out << "\n\n\n" << non_ss_faraways_result::vc_meta << endl; 
  }
  else
  {
    out << "\n\n\n" << non_ss_smiths_result::vc_meta << endl;     
  }
  
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


//============================================================================
// IMPLEMENTATION:  ss_smiths_faraways_test
//============================================================================
//
bool  ss_smiths_faraways_test::alt_calculated = false;
vector<ss_alt_result_ptr>  ss_smiths_faraways_test::alt_results;

// - Clear static members.
//
void
ss_smiths_faraways_test::clear()
{
  alt_calculated = false;
  alt_results.clear();
}

void
ss_smiths_faraways_test::calculate()
{
  if(my_instructions.valid)
  {
    const RPED::RefMPedInfo&  mp_info = my_mped.info();
    size_t  trait_index = mp_info.marker_find(my_instructions.trait);
    
    for(size_t marker_index = 0; marker_index < mp_info.marker_count(); ++marker_index)
    {
      if(marker_index == trait_index)
      {
        continue;
      }
      
      ss_smiths_faraways_result*  result = 0;
      
      if(my_type == sf_LINKAGE)
      {
        result = new ss_faraways_result;
      }
      else
      {
        result = new ss_smiths_result;
      }
      
      result->trait = my_instructions.trait;
      result->marker = mp_info.marker_info(marker_index).name();
      
      calculate_alt(trait_index, marker_index, *result);
      calculate_null(trait_index, marker_index, *result);      
      if(my_type == sf_LINKAGE)
      {
        calculate_posteriors(trait_index, marker_index, *result);
      }
      
      my_results.push_back(result_ptr(result)); 
    }
    
    alt_calculated = true;
    completed = true;
  }
}


// - Maximize likelihood while estimating the proportion of families w.
//   linkage (allowing for heterogeneity).
//
void
ss_smiths_faraways_test::calculate_alt(size_t trait_index, size_t marker_index, 
                                              ss_smiths_faraways_result& result)
{
  if(! alt_calculated)    // Keep a static copy of results.
  {
    ss_alt_result*  alt_result = new ss_alt_result;
    alt_result->marker = result.marker;  
  
    mle_sub_model  alt_mle(true, true);
    alt_mle.set_strict_limits();
    
    result.alt_ln_like = QNAN;
    
    if(likelihood_finite(alt_mle, trait_index, marker_index, "likelihood calculation assuming heterogeneity"))
    {
      MAXFUN::Results  data = maximize(alt_mle, trait_index, marker_index);
      result.alt_ln_like = max_checked_value(data);

      if(! SAGE::isnan(result.alt_ln_like))
      {
        result.alt_thetas.male_theta = alt_mle.male_theta();
        alt_result->alt_thetas.male_theta = result.alt_thetas.male_theta;      
        result.alt_thetas.female_theta = alt_mle.female_theta();
        alt_result->alt_thetas.female_theta = result.alt_thetas.female_theta;            
        result.alpha = alt_mle.alpha();
        alt_result->alpha = result.alpha;
        result.set_var_cov(data);            
        alt_result->var_cov = result.var_cov;                      
      }
    }
    
    alt_result->alt_ln_like = result.alt_ln_like;        
    
    result.alt_theta_ubs.male_theta = alt_mle.male_theta_ub();
    alt_result->alt_theta_ubs.male_theta = result.alt_theta_ubs.male_theta;      
    result.alt_theta_ubs.female_theta = alt_mle.female_theta_ub();  
    alt_result->alt_theta_ubs.female_theta = result.alt_theta_ubs.female_theta;
      
    alt_results.push_back(ss_alt_result_ptr(alt_result));
  }
  else    // Get results from static copy.
  {
    vector<ss_alt_result_ptr>::const_iterator  r_iter = find_if(alt_results.begin(), alt_results.end(), 
                                                               bind2nd(has_marker<ss_alt_result_ptr>(), result.marker));
    assert(r_iter != alt_results.end());
    
    result.alt_ln_like = (*r_iter)->alt_ln_like;
    result.alt_thetas.male_theta = (*r_iter)->alt_thetas.male_theta;
    result.alt_thetas.female_theta = (*r_iter)->alt_thetas.female_theta;    
    result.alpha = (*r_iter)->alpha;
    result.alt_theta_ubs.male_theta = (*r_iter)->alt_theta_ubs.male_theta;    
    result.alt_theta_ubs.female_theta = (*r_iter)->alt_theta_ubs.female_theta;
    result.var_cov = (*r_iter)->var_cov;        
  }
}

// - Maximize likelihood for the homogeneous case.
//
void
ss_smiths_faraways_test::calculate_null(size_t trait_index, size_t marker_index, 
                                                ss_smiths_faraways_result& result)
{
  mle_sub_model  null_mle(true, false);
  null_mle.set_strict_limits();
  
  if(my_type == sf_LINKAGE)
  {
    null_mle.set_male_theta(NULL_THETA);
    null_mle.set_female_theta(NULL_THETA);
    mped_calculator  null_mp_calc(my_mped, null_mle, trait_index, marker_index);
  
    result.null_ln_like = null_mp_calc.likelihood().get_log();
  }
  else
  {
    result.null_ln_like = QNAN;
    
    if(likelihood_finite(null_mle, trait_index, marker_index, "Smith's test"))
    {
      MAXFUN::Results  data = maximize(null_mle, trait_index, marker_index);
      result.null_ln_like = max_checked_value(data);
    }
      
    if(! SAGE::isnan(result.null_ln_like))
    {
      result.null_thetas.male_theta = null_mle.male_theta();
      result.null_thetas.female_theta = null_mle.female_theta();
    }
  }
  
  result.null_theta_ubs.male_theta = null_mle.male_theta_ub();
  result.null_theta_ubs.female_theta = null_mle.female_theta_ub();
}

// - Calculate posterior likelihoods for each pedigree as a whole AND for each
//   subpedigree w/i a given pedigree.  Subpedigrees may be identified
//   by a combination of pedigree id and name of an arbitrary member.
//   All calculations based on the alpha and theta(s) that was estimated for
//   the entire multipedigree.
//                                            -per discussion w. RCE 1/23/3.
//
void
ss_smiths_faraways_test::calculate_posteriors(size_t trait_index, size_t marker_index, 
                                                     ss_smiths_faraways_result& result)
{
  if(SAGE::isnan(result.alt_thetas.male_theta)   || 
     SAGE::isnan(result.alt_thetas.female_theta) ||
     SAGE::isnan(result.alpha)                     ) 
  {
    return;
  }

  mle_sub_model  mle(true, false);
  mle.set_male_theta(result.alt_thetas.male_theta);
  mle.set_female_theta(result.alt_thetas.female_theta);  

  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    string  ped_name = ped_iter->name();
    pedigree_posterior  ped_p(ped_name);
    double              ped_alt_ln_like = 0;
    double              ped_null_ln_like = 0;
  
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      peeler  p(*subped_iter, mle, trait_index, marker_index);
      subped_calculator  sp_calc(p);
      double  alt_ln_like  = sp_calc.likelihood().get_log();
      double  null_ln_like = sp_calc.unlinked_likelihood().get_log();
      
      double  term1 = result.alpha * exp(alt_ln_like);
      double  term2 = (1 - result.alpha) * exp(null_ln_like);
      double  posterior = term1 / (term1 + term2);
      
      string  member_name = subped_iter->member_begin()->name();
      
      ped_alt_ln_like += alt_ln_like;
      ped_null_ln_like += null_ln_like;
      ped_p.sub_posteriors.push_back(subpedigree_posterior(member_name, posterior)); 
    }
    
    double  ped_term1 = result.alpha * exp(ped_alt_ln_like);
    double  ped_term2 = (1 - result.alpha) * exp(ped_null_ln_like);
    ped_p.posterior = ped_term1 / (ped_term1 + ped_term2);
      
    result.posteriors.push_back(ped_p);
  }
}

void
ss_smiths_faraways_test::write_summary(ostream& out) const
{
assert(completed);

  if(my_type == sf_HOMOGENEITY)
  {
    out << "\n\n\n" << ss_smiths_result::meta
                    << ss_smiths_result::columns;
  }
  else 
  {
    out << "\n\n\n" << ss_faraways_result::meta
                    << ss_faraways_result::columns;  
  }
  
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
ss_smiths_faraways_test::write_detail(ostream& out) const
{
  assert(completed);

  if(my_type == sf_LINKAGE)
  {
    ios::fmtflags old_flags = out.flags();
  
    out << setfill(' ');
    out << left;
    
    // - Meta header.
    //
    out << setw(_LOCUS.tw()) << ""
        << POSTERIORS_FAM_ << "\n"  
        << setw(_LOCUS.tw()) << ""
        << EVIDENCE_ << "\n\n" << endl;
  
    pedigree_const_iterator  ped_iter;
    for(ped_iter = my_mped.pedigree_begin(); ped_iter != my_mped.pedigree_end(); ++ped_iter)
    {
      subpedigree_const_iterator  subped_iter;
      for(subped_iter = ped_iter->subpedigree_begin(); 
          subped_iter != ped_iter->subpedigree_end();
          ++subped_iter)
      {
        out << setw(_LOCUS.tw()) << ""
            << "Constituent Pedigree in Pedigree " << ped_iter->name() 
            << " Containing Member " << subped_iter->member_begin()->name()
            << endl;
            
        out << non_ss_faraways_result::detail_columns;
        for(size_t i = 0; i < my_results.size(); ++i)
        {
          ss_faraways_result* ptr = static_cast<ss_faraways_result*>(my_results[i].get());
          ptr->write_family_detail(out, ped_iter->name(), subped_iter->member_begin()->name());
        }
        
        out << "\n" << endl;
      }
    }    
  
    out.flags(old_flags); 
  }
  
  write_vc_matrix(out);
}

void
ss_smiths_faraways_test::write_vc_matrix(ostream& out) const
{
  assert(completed);
  
  if(my_type == sf_LINKAGE)
  {
    out << "\n\n\n" << ss_faraways_result::vc_meta << endl; 
  }
  else
  {
    out << "\n\n\n" << ss_smiths_result::vc_meta << endl;     
  }
  
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


