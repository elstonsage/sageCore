//============================================================================
// File:      genotype_results.ipp
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
// IMPLEMENTATION:  geno_probs
//============================================================================
//
inline
geno_probs::geno_probs()
      : aa(QNAN), ab(QNAN), bb(QNAN)
{}


//============================================================================
// IMPLEMENTATION:  genotype_result
//============================================================================
//
inline
genotype_result::genotype_result()
{}

inline
genotype_result::~genotype_result()
{}

inline void  
genotype_result::write_vc_matrix(ostream& out) const
{
  // INTENTIONALLY EMPTY
}

inline void  
genotype_result::write_summary(ostream& out) const
{
  // INTENTIONALLY EMPTY
}

inline void  
genotype_result::write_detail(ostream& out) const
{
  // INTENTIONALLY EMPTY
}

//============================================================================
// IMPLEMENTATION:  non_ss_genotype_result
//============================================================================
//
inline
non_ss_genotype_result::non_ss_genotype_result()
      : theta(QNAN)
{}      

//============================================================================
// IMPLEMENTATION:  ss_genotype_result
//============================================================================
//
inline
ss_genotype_result::ss_genotype_result()
      : thetas(QNAN, QNAN)
{}      

