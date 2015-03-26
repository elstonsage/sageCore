//============================================================================
// File:      parser.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/5/2 created        -djb
//                                                                          
// Notes:     Inline implementation of class, parser.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
inline const instructions&
parser::user_instructions() const
{
  return my_instructions;
}

inline size_t
parser::analysis_id() const
{
  return my_analysis_id;
}


//============================================================================
// IMPLEMENTATION:  lods_ptrs
//============================================================================
//
inline
parser::lods_ptrs::lods_ptrs()
      : option(0), sex_specific(0), male_female(0), average(0)
{}


