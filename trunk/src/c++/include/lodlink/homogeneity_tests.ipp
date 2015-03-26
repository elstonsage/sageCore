//============================================================================
// File:      homogeneity_tests.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/13/3 created        -djb
//                                                                          
// Notes:     Inline implementation of homogeneity test classes.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  non_ss_mortons_test
//============================================================================
//
inline
non_ss_mortons_test::non_ss_mortons_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                              const instructions& instr)
      : task(errors, mped, instr)
{}

inline void  
non_ss_mortons_test::announce_start() const
{
  cout << "Performing Morton's test for homogeneity ..." << endl;
}


//============================================================================
// IMPLEMENTATION:  ss_mortons_test
//============================================================================
//
inline
ss_mortons_test::ss_mortons_test(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                       const instructions& instr)
      : task(errors, mped, instr)
{}

inline void  
ss_mortons_test::announce_start() const
{
  cout << "Performing Morton's test for homogeneity ..." << endl;
}

