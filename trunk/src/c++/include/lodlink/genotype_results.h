#ifndef LODLINK_GENOTYPE_RESULTS_H
#define LODLINK_GENOTYPE_RESULTS_H
//============================================================================
// File:      genotype_results.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/21/3 - created.                                   djb
//                                                                          
// Notes:     results classes for genotype probabilities.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <map>
#include "lodlink/results.h"

namespace SAGE
{

namespace LODLINK
{

struct geno_probs
{
  geno_probs();

  double  aa;
  double  ab;
  double  bb;
};


//----------------------------------------------------------------------------
//  Class:    genotype_result
//                                                                          
//  Purpose:  base class for results of genotype probabilities tasks.
//                                                                          
//----------------------------------------------------------------------------
//
struct genotype_result : public task_result
{
  typedef FPED::FilteredMultipedigree::member_const_pointer  ind_ptr;
  typedef FPED::FilteredMultipedigree::member_type  member;
  typedef std::map<ind_ptr, geno_probs>  geno_map;
  typedef geno_map::const_iterator  map_iter;

  genotype_result();
  virtual ~genotype_result() = 0;

  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;  
  void  write_vc_matrix(ostream& out) const;
  
  // Data members.
  string  trait;
  string  marker;
  geno_map  probs;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_genotype_result
//                                                                          
//  Purpose:  results of sex-averaged genotype probabilities task.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_genotype_result : public genotype_result
{
  non_ss_genotype_result();

  static void  build_headers();

  void  write_genotypes(ostream& out, const member& ind) const;  
   
  // Data members.
  double  theta;
   
  static size_t  member_header_offset;
  static header  detail_columns;     
};


//----------------------------------------------------------------------------
//  Class:    ss_genotype_result
//                                                                          
//  Purpose:  base class for results of genotype probabilities tasks.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_genotype_result : public genotype_result
{
  ss_genotype_result();

  static void  build_headers();

  void  write_genotypes(ostream& out, const member& ind) const;
  
  // Data members.
  theta_pair  thetas;
   
  static size_t  member_header_offset;
  static header  detail_columns;
};


#include "lodlink/genotype_results.ipp"
}
}

#endif
