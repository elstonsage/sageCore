#ifndef LODLINK_GENOTYPES_H
#define LODLINK_GENOTYPES_H
//============================================================================
// File:      genotypes.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/21/3 - created.                                   djb
//                                                                          
// Notes:     classes for calculating genotype probabilities.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "lodlink/tasks.h"
#include "lodlink/genotype_results.h"

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    genotype_probs
//                                                                          
//  Purpose:  base class for classes calculating genotype probabilities.
//                                                                          
//----------------------------------------------------------------------------
//
class genotype_probs : public task
{
  public:
    typedef FPED::FilteredMultipedigree::member_type  member_type;
    
    void  announce_start() const;
    
  protected:
    enum  types { aa, ab, bb };
  
    genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    
    void  calculate_genotypes(peeler& plr, const member_type& ind, genotype_result& result);
    
    void  write_summary(ostream& out) const;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_genotype_probs
//                                                                          
//  Purpose:  calculating genotype probabilities for sex averaged recombin-
//            ation fraction.
//                                                                          
//----------------------------------------------------------------------------
//
class non_ss_genotype_probs : public genotype_probs
{
    friend void do_task_calculations<non_ss_genotype_probs, non_ss_genotype_result>(non_ss_genotype_probs&);

  public:
    non_ss_genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);

    void  calculate();
    void  write_detail(ostream& out) const;
      void  write_detail_header(ostream& out) const;
    
  private:
  
    // - Following two function names are misnomers in this context.  They are here to make
    //   possible the use of do_task_calculations().
    //  
    void  calculate_alt(size_t trait_index, size_t marker_index, non_ss_genotype_result& result);
      double  estimate_theta(size_t trait_index, size_t marker_index);
    void  calculate_null(size_t trait_index, size_t marker_index, non_ss_genotype_result& result);  
};


//----------------------------------------------------------------------------
//  Class:    ss_genotype_probs
//                                                                          
//  Purpose:  calculating genotype probabilities for sex specific recombin-
//            ation fractions.
//                                                                          
//----------------------------------------------------------------------------
//
class ss_genotype_probs : public genotype_probs
{
    friend void do_task_calculations<ss_genotype_probs, ss_genotype_result>(ss_genotype_probs&);  

  public:
    ss_genotype_probs(cerrorstream& errors, const FPED::FilteredMultipedigree& mped, const instructions& instr);
    
    void  calculate();
    void  write_detail(ostream& out) const;
      void  write_detail_header(ostream& out) const;
    
  private:
    
    // - Following two function names are misnomers in this context.  They are here to make
    //   possible the use of do_task_calculations().
    //
    void  calculate_alt(size_t trait_index, size_t marker_index, ss_genotype_result& result);
      theta_pair  estimate_thetas(size_t trait_index, size_t marker_index);    
    void  calculate_null(size_t trait_index, size_t marker_index, ss_genotype_result& result);  
};


#include "lodlink/genotypes.ipp"
}
}

#endif
