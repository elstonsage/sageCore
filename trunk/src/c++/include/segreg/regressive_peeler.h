#ifndef REGRESSIVE_PEELER_H
#define REGRESSIVE_PEELER_H
//===================================================================
//
//  File:	regressive_peeler.h
//
//  Author:	Stephen Gross
//
//  History:    sag Initial implementation		Jul 26 01
//		sag Added documentation.		Aug 02 01
//
//  Copyright (c) 2001, R. C. Elston
//  All rights reserved
//===================================================================

/** @file
 *		The regressive peeler is a specialized, inherited version
 *		of the peeler class. It is defined as:
 *
 *		class regressive_peeler : 
 *			public peeling::peeler<genotype_index,log_double>
 *
 *		It is designed to provide peeling functionality for a
 *		regressive model; hence it uses only genotype as its 
 *		datatype, and returns a probability value (a log_double).
 */


#include "peeling/peeler3.h"
#include "segreg/peeling_caches.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/LikelihoodElements.h"
#include "segreg/types/TypeDescription.h"

namespace SAGE { namespace SEGREG {

/** @class regressive_peeler
  *
  * Calculates Equations 11 - 22 (non-FPMM)
  *
  * The regressive_peeler calculates the equations given in _Mathematical
  * Notations, Assumptions, and Formula Derivations for Computation of the
  * Likelihood for a Pedigree with Regressive Models_.  In specific, it
  * calculates equations 11-22 from that document for the non-FPMM case.
  *
  * Some of the equations have been modified from that document to allow for
  * the programming, but for the most part, each programmatic function
  * calculates a single mathematical function.
  *
  * Some of the mathematical functions are (fill in rest):
  *                          
  * \f[ ant_i(u_i) = \left\{
  *     \begin{array}{ll}
  *       \sum_{u_m} \left\{  S_{mf}(u_m) 
  *       \sum_{u_f} \left[ (S_{fm}(u_f | u_m) SL1(i) ) \right] \right\} &
  *       \mbox{if $i$ is a nonfounder, and} \\
  *       SL1(i) & \mbox{if $i$ is a founder.}
  *     \end{array} \right.
  * \f]
  *
  * \f[ ant{^*}{_i}(u_i|u_s) = \left\{
  *     \begin{array}{ll}
  *       \sum_{u_m} \left\{ S_{mf}(u_m) \sum_{u_f}
  *       \left[ S_{fm}(u_f|u_m) SL3(i) \right] \right\} &
  *       \mbox{if $i$ is a nonfounder, and} \\
  *       SL3(i) & \mbox{if $i$ is a founder.}
  *     \end{array} \right.
  * \f]
  *
  * \f[ S_{is}(u_i) = ant_i(u_i) \prod_{j \in S_i,j \not = s} pos_{ij}(u_i)
  *   \f]
  *
  * \f[ S_{is}(u_i|u_s) = ant{^*}{_i(u_i|u_s)} \prod_{j \in S_i,j \not = s}
  *     pos_{ij}(u_i) \f]
  */
class regressive_peeler : public 
	peeling::peeler<TypeDescription::State,log_double,
        peeling::individual_cache<TypeDescription::State,log_double> >
{
  public:

    //==================================================================
    // Typedefs:
    //==================================================================

    typedef RegPenetranceCalculator::penetrance_info   penetrance_info;

    //==================================================================
    // Constructor:
    //==================================================================

    regressive_peeler(const subped_type&        subped,
                      const LikelihoodElements& lelt);

    ~regressive_peeler() { } 

    //==================================================================
    // Public get functions:
    //==================================================================

    const result_type & anterior_with_mate         (const member_type& ind, const member_type& mate, 
   	                                            const data_type& g,   const data_type& h);
    const result_type & partial_parental_likelihood(const member_type& ind, const member_type& mate, 
                                                    const data_type& g);
    const result_type & partial_parental_likelihood(const member_type& ind, const member_type& mate, 
                                                    const data_type& g,   const data_type& h);
    
    //==================================================================
    // Internal equations:
    //==================================================================
  protected:
    const result_type & internal_partial_parental_likelihood
                                                      (const member_type& ind, const member_type& mate, 
                                                       const data_type& g,   result_type& ippl1);
    const result_type & internal_partial_parental_likelihood
                                                      (const member_type& ind, const member_type& mate,
                                                       const data_type& g,   const data_type&   h, 
                                                                             result_type&       ippl2);
    const result_type & internal_anterior             (const member_type& ind, const data_type&   g,    
                                                                             result_type&       ia);
    const result_type & internal_anterior_with_mate   (const member_type& ind, const member_type& mate,
                                                       const data_type& g,   const data_type&   h,
                                                                             result_type&       iawm);
    const result_type & internal_anterior_terminal    (const member_type& ind, const data_type&   g, 
                                                                             result_type&       iatg);
    const result_type & internal_posterior            (const member_type& ind, const data_type&   g, 
                                                                             result_type&       ip);
    const result_type & internal_posterior_with_mate  (const member_type& ind, const member_type& mate,
                                                       const data_type& g,   result_type&       ipwm);
    const result_type & internal_posterior_except_mate(const member_type& ind, const member_type& mate, 
                                                       const data_type& g,   result_type&       ipem);
    const result_type & internal_posterior_terminal   (const member_type& ind, const data_type&   g,
                                                                             result_type&       ipt);

    log_double calculate_founder_anterior
                                (const member_type&          ind,
                                 const data_type&            g,
                                 const PenetranceContext&    context);
    log_double calculate_nonfounder_anterior
                                (const member_type&          ind,
                                 const data_type&            g,
                                 const PenetranceContext&    context);
    double     nonfounder_likelihood
                                (const member_type&          ind,
                                 const data_type&            g,
                                 const PenetranceContext&    context);

    log_double sibship_likelihood
                                (const PenetranceContext&    context);
    log_double sibship_likelihood_excluding_sib
                                (const PenetranceContext&    context,
                                 const member_type&          excluded_sib);
                                                                             
    //==================================================================
    // Private data members:
    //==================================================================
    
    const LikelihoodElements&                   my_likelihood_elts;
};

// end namespace
}
}

#include "segreg/regressive_peeler.ipp"

#endif
