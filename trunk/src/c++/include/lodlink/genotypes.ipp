//============================================================================
// File:      genotypes.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/21/3 - created.                                djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  genotype_probs
//============================================================================
//
inline
genotype_probs::genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                     const instructions& instr)
      : task(errors, mped, instr)
{}

inline void  
genotype_probs::announce_start() const
{
  cout << "Calculating genotype probabilities ..." << endl;
}


//============================================================================
// IMPLEMENTATION:  non_ss_genotype_probs
//============================================================================
//
inline
non_ss_genotype_probs::non_ss_genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                                   const instructions& instr)
      : genotype_probs(errors, mped, instr)
{}


//============================================================================
// IMPLEMENTATION:  ss_genotype_probs
//============================================================================
//
inline
ss_genotype_probs::ss_genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, 
                                           const instructions& instr)
      : genotype_probs(errors, mped, instr)
{}
