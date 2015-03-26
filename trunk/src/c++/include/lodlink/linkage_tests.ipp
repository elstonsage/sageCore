//============================================================================
// File:      linkage_tests.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/25/2 created        -djb
//                                                                          
// Notes:     Inline implementation of linkage test classes.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  non_ss_lod_ratio_test
//============================================================================
//
inline
non_ss_lod_ratio_test::non_ss_lod_ratio_test(cerrorstream& errors, 
                                             const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : task(errors, mped, instr)
{}

inline
non_ss_lod_ratio_test::~non_ss_lod_ratio_test()
{}

inline void  
non_ss_lod_ratio_test::announce_start() const
{
  cout << "Performing lod ratio test for linkage ..." << endl;
}


//============================================================================
// IMPLEMENTATION:  ss_lod_ratio_test
//============================================================================
//
inline
ss_lod_ratio_test::ss_lod_ratio_test(cerrorstream& errors, 
                                     const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : task(errors, mped, instr)
{}

inline
ss_lod_ratio_test::~ss_lod_ratio_test()
{}

inline void  
ss_lod_ratio_test::announce_start() const
{
  cout << "Performing lod ratio test for linkage ..." << endl;
}


//============================================================================
// IMPLEMENTATION:  cleves_elston_test
//============================================================================
//
inline
cleves_elston_test::cleves_elston_test(cerrorstream& errors, 
                                       const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : task(errors, mped, instr)
{}

inline
cleves_elston_test::~cleves_elston_test()
{}

inline void  
cleves_elston_test::announce_start() const
{
  cout << "Performing Cleves-Elston test for linkage ..." << endl;
}

