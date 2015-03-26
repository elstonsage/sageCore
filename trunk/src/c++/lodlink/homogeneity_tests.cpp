//============================================================================
// File:      homogeneity_tests.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   1/13/3 - created.                         djb
//                                                                          
// Notes:     Implementation of linkage homogeneity test classes.   
//               
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/homogeneity_tests.h"

namespace SAGE
{

namespace LODLINK
{

// - For Morton's tests, estimate recombination fraction(s) between 0 and 1 inclusive.
//   -RCE 2-3-3
//
//============================================================================
// IMPLEMENTATION:  non_ss_mortons_test
//============================================================================
//
void
non_ss_mortons_test::calculate()
{
  do_task_calculations<non_ss_mortons_test, non_ss_mortons_result>(*this);
}

// - Maximize likelihood over all groups w. recombination fraction for
//   each group estimated separately.
//
void
non_ss_mortons_test::calculate_alt(size_t trait_index, size_t marker_index, 
                                              non_ss_mortons_result& result)
{
  assert(! my_instructions.groups.empty());
  
  double total_ln_like = 0;
  map<string, group>::const_iterator  g_iter = my_instructions.groups.begin();
  for(; g_iter != my_instructions.groups.end(); ++g_iter)
  {
    mle_sub_model  alt_mle(false, false);
    alt_mle.set_relaxed_limits();
    result.group_theta_ub = alt_mle.average_theta_ub();  
    
    double  group_ln_like = QNAN;
    bool  success = false;
    
    if(likelihood_finite(g_iter->first, g_iter->second, alt_mle, trait_index, marker_index, "Morton's test"))
    {
      //cout << "\ngroup address " << &(g_iter->second) << endl;
      //cout << "sub_model_address " << &alt_mle << endl;
    
      MAXFUN::Results  data = maximize(g_iter->second, alt_mle, trait_index, marker_index);
      group_ln_like = max_checked_value(data); 
      
      success = result.group_results.insert(make_pair(g_iter->first, 
                                make_pair(alt_mle.average_theta(), group_ln_like))).second;
    }
    else
    {
      success = result.group_results.insert(make_pair(g_iter->first, 
                                make_pair(QNAN, group_ln_like))).second;
    }
    
    assert(success);
    total_ln_like += group_ln_like;
  }
  
  result.alt_ln_like = total_ln_like;
}

// - Maximize likelihood over entire multipedigree.
//
void
non_ss_mortons_test::calculate_null(size_t trait_index, size_t marker_index, 
                                               non_ss_mortons_result& result)
{
  mle_sub_model  null_mle(false, false);
  null_mle.set_relaxed_limits();
  
  result.null_ln_like = QNAN;
  if(likelihood_finite(null_mle, trait_index, marker_index, "Morton's test"))
  {
    MAXFUN::Results  data = maximize(null_mle, trait_index, marker_index);
    result.null_ln_like = max_checked_value(data);
  }
  
  if(! SAGE::isnan(result.null_ln_like))
  {
    result.null_theta = null_mle.average_theta(); 
  }
  
  result.null_theta_ub = null_mle.average_theta_ub();  
}

void
non_ss_mortons_test::write_summary(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << non_ss_mortons_result::meta
                  << non_ss_mortons_result::columns;
              
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
non_ss_mortons_test::write_detail(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ') << setw(_LOCUS.tw()) << ""
      << left  << "Morton's Homogeneity Test" << "\n\n" << endl;   
  
  // - For each group.
  //
  map<string, group>::const_iterator g_iter = my_instructions.groups.begin();
  for(; g_iter != my_instructions.groups.end(); ++g_iter)
  {
    out << setfill(' ') << setw(_LOCUS.tw()) << ""
        << left  << "Group " << g_iter->first << "\n\n" 
        << non_ss_mortons_result::detail_columns; 
    
    // - For ea. locus.
    //
    for(size_t r = 0; r < my_results.size(); ++r)
    {
      non_ss_mortons_result* ptr = static_cast<non_ss_mortons_result*>(my_results[r].get());
      ptr->write_group_detail(out, g_iter->first);
    }
    
    out << "\n" << endl;
  }
  
  out.flags(old_flags);  
}

//============================================================================
// IMPLEMENTATION:  ss_mortons_test
//============================================================================
//
void
ss_mortons_test::calculate()
{
  do_task_calculations<ss_mortons_test, ss_mortons_result>(*this);
}

// - Maximize likelihood over all groups w. recombination fractions for
//   each group estimated separately.
//
void
ss_mortons_test::calculate_alt(size_t trait_index, size_t marker_index, 
                                              ss_mortons_result& result)
{
  assert(! my_instructions.groups.empty());
  
  double total_ln_like = 0;
  map<string, group>::const_iterator  g_iter = my_instructions.groups.begin();
  for(; g_iter != my_instructions.groups.end(); ++g_iter)
  {
    mle_sub_model  alt_mle(true, false);
    alt_mle.set_relaxed_limits();
    result.group_theta_ubs.male_theta = alt_mle.male_theta_ub();
    result.group_theta_ubs.female_theta = alt_mle.female_theta_ub();
    
    double  group_ln_like = QNAN;
    bool  success = false;
    
    if(likelihood_finite(g_iter->first, g_iter->second, alt_mle, trait_index, marker_index, "Morton's test"))
    {
      MAXFUN::Results  data = maximize(g_iter->second, alt_mle, trait_index, marker_index);
      group_ln_like = max_checked_value(data); 
    
      success = result.group_results.insert(make_pair(g_iter->first, 
                                  make_pair(theta_pair(alt_mle.male_theta(), alt_mle.female_theta()), 
                                  group_ln_like))).second;
    }
    else
    {
      success = result.group_results.insert(make_pair(g_iter->first, 
                                  make_pair(theta_pair(QNAN, QNAN), 
                                  group_ln_like))).second;
    }
    
    assert(success);
    total_ln_like += group_ln_like;
  }
  
  result.alt_ln_like = total_ln_like;
}

// - Maximize likelihood over entire multipedigree.
//
void
ss_mortons_test::calculate_null(size_t trait_index, size_t marker_index, 
                                               ss_mortons_result& result)
{
  mle_sub_model  null_mle(true, false);
  null_mle.set_relaxed_limits();
  
  result.null_ln_like = QNAN;
  if(likelihood_finite(null_mle, trait_index, marker_index, "Morton's index"))
  {
    MAXFUN::Results  data = maximize(null_mle, trait_index, marker_index);
    result.null_ln_like = max_checked_value(data);
  }
  
  if(! SAGE::isnan(result.null_ln_like))
  {
    result.null_thetas.male_theta = null_mle.male_theta();
    result.null_thetas.female_theta = null_mle.female_theta();
  }
  
  result.null_theta_ubs.male_theta = null_mle.male_theta_ub();
  result.null_theta_ubs.female_theta = null_mle.female_theta_ub();
}

void
ss_mortons_test::write_summary(ostream& out) const
{
  assert(completed);
  
  out << "\n\n\n" << ss_mortons_result::meta
                  << ss_mortons_result::columns;
              
  for(size_t r = 0; r < my_results.size(); ++r)
  {
    my_results[r]->write_summary(out);
  }
}

void
ss_mortons_test::write_detail(ostream& out) const
{
  ios::fmtflags old_flags = out.flags();
  
  out << setfill(' ') << setw(_LOCUS.tw()) << ""
      << left  << "Morton's Homogeneity Test" << "\n\n" << endl;   
  
  // - For ea. group.
  //
  map<string, group>::const_iterator g_iter = my_instructions.groups.begin();  
  for(; g_iter != my_instructions.groups.end(); ++g_iter)
  {
    out << setfill(' ') << setw(_LOCUS.tw()) << ""
        << left  << "Group " << g_iter->first << "\n\n" 
        << ss_mortons_result::detail_columns; 
    
    // - For ea. locus.
    //
    for(size_t r = 0; r < my_results.size(); ++r)
    {
      ss_mortons_result* ptr = static_cast<ss_mortons_result*>(my_results[r].get());
      ptr->write_group_detail(out, g_iter->first);
    }
    
    out << "\n" << endl;
  }
  
  out.flags(old_flags);  
}

}
}



