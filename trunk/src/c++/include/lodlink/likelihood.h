#ifndef LODLINK_LIKELIHOOD_H
#define LODLINK_LIKELIHOOD_H
//============================================================================
// File:      likelihood.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/2/2 - created.                                   djb
//                                                                          
// Notes:     defines calculator classes for calculating subpedigree and  
//            multi-pedigree likelihoods.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "fped/fped.h"
#include "maxfun/maxfun.h"
#include "lodlink/peeler.h"
#include "lodlink/instructions.h"

namespace SAGE
{

namespace LODLINK
{

void  check_parameters(const MaxFunction::parameter_vector& theta, const mle_sub_model& mle);

//----------------------------------------------------------------------------
//  Class:    subped_calculator
//                                                                          
//  Purpose:  calculate likelihood of a subpedigree.
//                                                                          
//----------------------------------------------------------------------------
//
class subped_calculator : public MaxFunction
{
  public:
    subped_calculator(peeler& p);
    
    log_double  likelihood();
    log_double  unlinked_likelihood();    

  private:
    virtual double evaluate(parameter_vector& theta);
    virtual int    update_bounds(parameter_vector& theta);

    log_double  calc_unlinked_likelihood();
    
    subped_calculator(const subped_calculator& other);
    subped_calculator& operator=(const subped_calculator& other);
    
    // Data members.
    peeler&  my_peeler;
    
    log_double  my_unlinked_likelihood;
    bool        unlinked_likelihood_cached;
};


//----------------------------------------------------------------------------
//  Class:    ped_calculator
//                                                                          
//  Purpose:  calculate likelihood of a pedigree.
//                                                                          
//----------------------------------------------------------------------------
//
class ped_calculator : public MaxFunction
{
  public:
    typedef SAGE::FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_const_iterator;
  
    ped_calculator(const FPED::Pedigree& ped, const mle_sub_model& mle,
                    size_t trait, size_t marker);
    log_double  likelihood();
    log_double  unlinked_likelihood();

  private:
    virtual double evaluate(parameter_vector& theta);
    virtual int    update_bounds(parameter_vector& theta);

    log_double  calc_unlinked_likelihood();
    
    ped_calculator(const ped_calculator& other);
    ped_calculator& operator=(const ped_calculator& other);
    
    // Data members.
    const FPED::Pedigree&  my_ped;
    const mle_sub_model&   my_mle;
    size_t                 my_trait;
    size_t                 my_marker;
    
    log_double  my_unlinked_likelihood;
    bool        unlinked_likelihood_cached;    
};

//----------------------------------------------------------------------------
//  Class:    group_calculator
//                                                                          
//  Purpose:  calculate likelihood of an arbitrary set of pedigrees.
//                                                                          
//----------------------------------------------------------------------------
//
typedef FPED::FilteredMultipedigree::pedigree_const_pointer  pedigree_const_pointer;

class filtered_ped_cmp
{
  public:
    bool  operator()(const pedigree_const_pointer arg1,
                     const pedigree_const_pointer arg2 ) const
    {
      return  arg1->name() < arg2->name();
    }
};

typedef set<pedigree_const_pointer, filtered_ped_cmp>  filtered_group;

class group_calculator : public MaxFunction
{
  public:
    group_calculator(const group& g, const FPED::FilteredMultipedigree& mped, const mle_sub_model& mle,
                    size_t trait, size_t marker);
    log_double  likelihood();

  private:
    void  build_group(const group& g, const FPED::FilteredMultipedigree& mped);
    double evaluate(parameter_vector& theta);
    int    update_bounds(parameter_vector& theta);

    group_calculator(const group_calculator& other);
    group_calculator& operator=(const group_calculator& other);
    
    // Data members.
    filtered_group         my_group;
    const mle_sub_model&   my_mle;
    size_t                 my_trait;
    size_t                 my_marker;
};


//----------------------------------------------------------------------------
//  Class:    mped_calculator
//                                                                          
//  Purpose:  calculate likelihood of a multi-pedigree.
//                                                                          
//----------------------------------------------------------------------------
//
class mped_calculator : public MaxFunction
{
  public:
    typedef SAGE::FPED::FilteredMultipedigree::pedigree_const_iterator     pedigree_const_iterator;
    typedef SAGE::FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_const_iterator;
  
    mped_calculator(const FPED::FilteredMultipedigree& mped, const mle_sub_model& mle,
                    size_t trait, size_t marker);
    log_double  likelihood();
    log_double  unlinked_likelihood();

  private:
    virtual double evaluate(parameter_vector& theta);
    virtual int    update_bounds(parameter_vector& theta);

    log_double  calc_unlinked_likelihood();
    
    mped_calculator(const mped_calculator& other);
    mped_calculator& operator=(const mped_calculator& other);
    
    // Data members.
    const FPED::FilteredMultipedigree&  my_mped;
    const mle_sub_model&     my_mle;
    size_t                   my_trait;
    size_t                   my_marker;
    
    log_double  my_unlinked_likelihood;
    bool        unlinked_likelihood_cached;    
};

#include "lodlink/likelihood.ipp"

}
}

#endif



