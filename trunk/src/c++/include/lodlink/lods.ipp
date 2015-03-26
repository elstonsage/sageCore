//============================================================================
// File:      lods.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/13/2 created        -djb
//                                                                          
// Notes:     Inline implementation of lod score calculations classes.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  non_ss_lods
//============================================================================
//
inline
non_ss_lods::non_ss_lods(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : task(errors, mped, instr)
{}

inline void  
non_ss_lods::announce_start() const
{
  cout << "Performing lod score calculations ..." << endl;
}


//============================================================================
// IMPLEMENTATION:  ss_lods
//============================================================================
//
inline
ss_lods::ss_lods(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr)
      : task(errors, mped, instr)
{}

inline void  
ss_lods::announce_start() const
{
  cout << "Performing lod score calculations ..." << endl;
}



