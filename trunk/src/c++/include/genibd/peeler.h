#ifndef GENIBD_PEELER_H
#define GENIBD_PEELER_H

//============================================================================
// File:      peeler.h
//
// Author:    Dan Baechle
//
// History:   10/18/2 - created.                                 djb
//            Modified for GENIBD                                yjs Apr 2004
//
// Notes:     defines a peeling class for calculating pedigree 
//            anterior and posterior likelihoods.
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

//#define DEBUG_PEELER

#include <string>
#include "peeling/peeler3.h"
#include "genibd/cache.h"
#include "genibd/trans_calculator.h"

using std::string;

namespace SAGE
{

namespace GENIBD
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
class peeler : public peeling::peeler<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>
{
  public:

    peeler(const subped_type& subped, const inheritance_model& model);
    
    // Accessors.
    const subped_type&       subpedigree() const;
    const inheritance_model& marker() const;
  
    const log_double&  
    internal_anterior(const member_type& ind, 
                      const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
    const log_double&  
    internal_posterior(const member_type& ind, 
                       const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
    const log_double&  
    internal_posterior_with_mate(const member_type& ind, const member_type& mate, 
                                 const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
    const log_double&  
    internal_posterior_except_mate(const member_type& ind, const member_type& mate, 
                                   const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
    const log_double&  
    internal_anterior_terminal(const member_type& ind, 
                               const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
    const log_double&  
    internal_posterior_terminal(const member_type& ind, 
                                const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result);
  
    log_double  likelihood_of_sibs_except_sib(const ind_genotype& mg, const ind_genotype& fg,
                                              const member_type& ind, const member_type& ex_sib);

    log_double  likelihood_of_sibs(const ind_genotype& mg, const ind_genotype& fg,
                                   const member_type& ind);

  private:

    log_double  internal_anterior_phased(const member_type& ind, 
                                         const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, bool flip_allele = false);

    log_double  likelihood_of_offspring(const ind_genotype& ind_g, const ind_genotype& mate_g,
                                        const member_type& ind, const member_type& mate);
    log_double  sum_ind_and_posterior(const ind_genotype& mg, const ind_genotype& fg,
                                      const member_type& ind);
    
    peeler(const peeler& other);
    peeler&  operator=(const peeler& other);
                                  
    // Debugging functions.
    static void  print_enter(const member_type& ind, 
                             const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, const string& func_name);
    static void  print_enter_w_mate(const member_type& ind, const member_type& mate,
                                    const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, const string& func_name);
    static void  print_exit(double result, const string& func_name);
    
    // Data members.
    //
    const inheritance_model&   my_model;
    trans_calculator           my_tcalc;
};

#include "genibd/peeler.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
