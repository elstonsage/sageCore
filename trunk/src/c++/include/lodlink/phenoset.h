#ifndef LODLINK_PHENOSET_H
#define LODLINK_PHENOSET_H
//============================================================================
// File:      phenoset.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/28/2 - created.                                   djb
//                                                                          
// Notes:     Defines a class to represent the penetrant joint genotypes
//            at two loci for a given individual.      
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "fped/fped.h"
#include "mlocus/penmodel.h"
#include "lodlink/geno_elim.h"
#include "lodlink/definitions.h"

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    phenoset
//                                                                          
//  Purpose:  represents the set of all valid joint genotypes for a given
//            individual.
//                                                                          
//----------------------------------------------------------------------------
//
class phenoset
{
  public:
    typedef FPED::FilteredMultipedigree::member_type  member_type;
    typedef FPED::FilteredMultipedigree::subpedigree_type  subped_type;
    typedef MLOCUS::penetrance_model::phased_penetrance_iterator  phased_penetrance_iterator;  
  
    class phenoset_iterator
    {
      public:
        friend class phenoset;
      
        joint_pen_iter      operator*() const;
        phenoset_iterator&  operator++();
        bool                operator==(const phenoset_iterator& other) const;
        bool                operator!=(const phenoset_iterator& other) const;
      
      private:
        phenoset_iterator(const phased_penetrance_iterator& trait_iter,
                          const phased_penetrance_iterator& marker_iter,
                          const phenoset* ph_set);
        
        // Data members.
        phased_penetrance_iterator  my_trait_iter;
        phased_penetrance_iterator  my_marker_iter;
        const phenoset*             my_phenoset;
    };
    
    friend class phenoset_iterator;

    phenoset(const subped_type& sp, size_t trait, size_t marker, const member_type& ind);
    
    phenoset_iterator  begin() const;
    phenoset_iterator  end() const;
    
  private:
    phased_penetrance_iterator  trait_iter_begin() const;
    phased_penetrance_iterator  marker_iter_begin() const;
    phased_penetrance_iterator  trait_iter_end() const;
    phased_penetrance_iterator  marker_iter_end() const;
  
    // Data members.
    size_t  my_trait;
    size_t  my_marker;
    const member_type&  my_ind;
    
    MLOCUS::penetrance_model  my_trait_pm;
    MLOCUS::penetrance_model  my_marker_pm;
    size_t  my_trait_phenotype;
    size_t  my_marker_phenotype;
};


#include "lodlink/phenoset.ipp"
}
}

#endif


