#ifndef LODLINK_PEELER_H
#define LODLINK_PEELER_H
//============================================================================
// File:      peeler.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                   djb
//                                                                          
// Notes:     defines a peeling class for calculating pedigree 
//            anterior and posterior likelihoods.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//#define DEBUG_PEELER

#include <string>
#include "fped/fped.h"
#include "peeling/peeler3.h"
#include "lodlink/cache.h"
#include "lodlink/trans_calculator.h"
#include "lodlink/phenoset.h"

using std::string;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    peeler
//                                                                          
//  Purpose:  provide functions for calculating pedigree likelihood for two
//            point linkage by the method of Fernando, Stricker and Elston.  
//            1993.
//                                                                          
//----------------------------------------------------------------------------
//
class peeler : public peeling::peeler<joint_pen_iter, log_double>
{
  public:
    typedef FPED::FilteredMultipedigree::member_type               member_type;  
    typedef FPED::FilteredMultipedigree::subpedigree_type          subped_type;
    typedef FPED::FilteredMultipedigree::member_const_pointer      member_const_pointer;
    typedef FPED::FilteredMultipedigree::sibling_const_iterator    sibling_const_iterator;
    typedef FPED::FilteredMultipedigree::offspring_const_iterator  offspring_const_iterator;
    typedef FPED::FilteredMultipedigree::mate_const_iterator       mate_const_iterator;
    
    typedef MLOCUS::penetrance_model::phased_penetrance_iterator  phased_penetrance_iterator;
    typedef phenoset::phenoset_iterator  phenoset_iterator;
  
    peeler(const subped_type& subped, const mle_sub_model& mle,
           size_t trait, size_t marker);
    
    // Accessors.
    const subped_type&       subpedigree() const;
    size_t                   trait() const;
    size_t                   marker() const;
    const trans_calculator&  tcalc() const;
    
  
    const log_double&  
    internal_anterior(const member_type& ind, 
                      const joint_pen_iter& jpi, log_double& result);
    const log_double&  
    internal_posterior(const member_type& ind, 
                       const joint_pen_iter& jpi, log_double& result);
    const log_double&  
    internal_posterior_with_mate(const member_type& ind, const member_type& mate, 
                                 const joint_pen_iter& jpi, log_double& result);
    const log_double&  
    internal_posterior_except_mate(const member_type& ind, const member_type& mate, 
                                   const joint_pen_iter& jpi, log_double& result);
    const log_double&  
    internal_anterior_terminal(const member_type& ind, 
                               const joint_pen_iter& jpi, log_double& result);
    const log_double&  
    internal_posterior_terminal(const member_type& ind, 
                                const joint_pen_iter& jpi, log_double& result);
  
  private:
    log_double  likelihood_of_offspring(const joint_genotype& ind_jg, const joint_genotype& mate_jg,
                                        const member_type& ind, const member_type& mate);
    log_double  likelihood_of_sibs(const joint_genotype& mjg, const joint_genotype&,
                                   const member_type& ind);
    log_double  sum_ind_and_posterior(const joint_genotype& mjg, const joint_genotype& fjg,
                                      const member_type& ind);
    
    peeler(const peeler& other);
    peeler&  operator=(const peeler& other);
                                  
    // Debugging functions.
    static void  print_enter(const member_type& ind, 
                      const joint_pen_iter& jpi, const string& func_name);
    static void  print_enter_w_mate(const member_type& ind, const member_type& mate,
                             const joint_pen_iter& jpi, const string& func_name);
    static void  print_exit(double result, const string& func_name);
    
    // Data members.
    trans_calculator  my_tcalc;
    size_t  my_trait;
    size_t  my_marker;
};

#include "lodlink/peeler.ipp"

}
}

#endif



