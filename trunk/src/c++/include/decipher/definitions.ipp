//============================================================================
// File:      defintitions.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/28/5 created        -djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


inline string
double_to_string(double value)
{
  string  string_value = doub2str(value);
  
  if(string_value == "inf")
  {
    string_value = "+infinity";
  }
  else if(string_value == "-inf")
  {
    string_value = "-infinity";
  }
  else if(string_value == "nan")
  {
    string_value = "---";
  }
  
  return string_value;
}


//============================================================================
// IMPLEMENTATION:  output_state
//============================================================================
//
inline
output_state::output_state()
      : msg_interrupted(false), allow_progress_msg(true), first_block(true)
{}


//============================================================================
// IMPLEMENTATION:  p_values
//============================================================================
//
inline
p_values::p_values(const string& name, double asymp, double empir)
      : outer_sub_pop_name(""), asymptotic(QNAN), composite_ln_like(QNAN),
        whole_ln_like(QNAN), test_statistic(QNAN), degrees_of_freedom((size_t)(-1)),
        empirical(QNAN), permutations(0)
{}


//============================================================================
// IMPLEMENTATION:  ld_record
//============================================================================
//
inline
ld_record::ld_record(const string& m1, const string& m2, double linkage_disequilibrium)
      : marker1(m1), marker2(m2), ld(linkage_disequilibrium)
{}

