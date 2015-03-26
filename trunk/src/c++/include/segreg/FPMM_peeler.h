#ifndef FPMM_PEELER_H
#define FPMM_PEELER_H
//===================================================================
//
//  File:	FPMM_peeler.h
//
//  Author:	Stephen Gross
//
//  History:    sag Initial implementation		Jul 26 01
//		sag Added documentation.		Aug 02 01
//
//  Copyright (c) 2001, R. C. Elston
//  All rights reserved
//===================================================================



#include "peeling/peeler3.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/SL_calculator.h"
#include "segreg/peeling_caches.h"
#include "segreg/model.h"

namespace SAGE {
namespace SEGREG {

//======================================================================
//
//  Class: FPMM_peeler
//
//======================================================================
/** The FPMM peeler is a specialized, inherited version of the peeler class.
 *  It is defined as:
 *
 *  class FPMM_peeler : public peeling::peeler<genetic_info,log_double>
 *
 *  It is designed to provide peeling functionality for an FPMM model; hence
 *  it uses the struct genetic_info as its datatype, and returns a
 *  probability value (a log_double). The genetic_info struct is defined in
 *  segreg_datatypes.h, and has two data members: genotype, and
 *  polygenotype.
 */
class FPMM_peeler : public 
	peeling::peeler<genetic_info,log_double,
        peeling::individual_cache<genetic_info,log_double> >
{
  public:
    typedef polygenic_penetrance_calculator::penetrance_info   penetrance_info;

    int tabs;
    int thresh;
    void t1(const string & s)
    {
/*
      tabs++;
      if(tabs < thresh)
      {
        for(int i = 0; i < tabs; i++) cout << "--";
        cout << s << endl;
      }
*/
    }

    void t2(const string & s)
    {
/*
      if(tabs < thresh)
      {
        for(int i = 0; i < tabs; i++) cout << "--";
        cout << s << endl;
      }
      tabs--;
*/
    }


//    FPMM_peeler(const subped_type& subped);
    FPMM_peeler(const subped_type& subped, const model & modp);

    /// Public set function for using the FPMM SL calculator:
    void set_SL(FPMM_SL * SL_calculator) { SL_calc = SL_calculator; } 

    log_double partial_parental_likelihood(
				       const member_type&         indiv,
                                       genotype_index genotype_indiv,
                                       size_t         polygenotype_indiv,
                                       const member_type&         spouse);
    log_double partial_parental_likelihood(
				       const member_type&         indiv,
                                       genotype_index genotype_indiv,
                                       size_t         polygenotype_indiv,
                                       const member_type&         spouse,
                                       genotype_index genotype_spouse,
                                       size_t         polygenotype_spouse);

  private:

    //==================================================================
    // Likelihood equations:
    //==================================================================


    const result_type & internal_anterior              (const member_type& ind, const data_type & g,
                                                        result_type & ia);
    const result_type & internal_anterior_with_mate    (const member_type& ind, 
                                                        const member_type& spouse, 
                                                        const data_type & g, 
                                                        const data_type & h,
                                                        result_type & iawm);
    const result_type & internal_anterior_terminal     (const member_type& ind,
	                                                const data_type & g, 
                                                        result_type & iatg);
    const result_type & internal_posterior             (const member_type& ind, const data_type & g,
	                                                result_type & ip);
    const result_type & internal_posterior_with_mate   (const member_type& ind, const member_type& mate,
	                                                const data_type & g, 
                                                        result_type & ipwm);
    const result_type & internal_posterior_except_mate (const member_type& ind, const member_type& mate,
	                                                const data_type & g,
                                                        result_type & ipem);
    const result_type & internal_posterior_terminal    (const member_type& ind,
                                                        const data_type & g,
                                                        result_type & ipt);

    const model   * my_model;
          FPMM_SL * SL_calc;
};

}} // End namespace

#endif
