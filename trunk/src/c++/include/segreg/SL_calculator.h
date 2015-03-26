#ifndef SL_CALCULATOR_H
#define SL_CALCULATOR_H
//========================================================================
//  File: 	SL_calculator.h
//
//  Author:	Kai He, Stephen Gross
//
//  History:	sag  Initial implementation.		Jul 20 2001
//		sag  Added documentation.		Aug 02 2001
//
//  Notes:	See doxygen comments
//
// Copyright (c) 2001 R. C. Elston
// All rights reserved
//========================================================================



#include <mped/sp.h>
#include <mped/mp.h>
#include "peeling/peeler3.h"
#include "numerics/log_double.h"
#include "segreg/model.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/freq_sub_model.h"
#include "segreg/transmission_sub_model.h"
#include "segreg/polygenic_penetrance_calculator.h"
#include "segreg/peeling_caches.h"
#include "segreg/polygenic_transition_calculator.h"

namespace SAGE {
namespace SEGREG {

class FPMM_SL
{
  public:
    typedef FPED::Member               member_type;
    typedef FPED::MemberIterator       member_iterator;
    typedef FPED::MemberPointer        member_pointer;
    typedef FPED::MemberConstPointer   member_const_pointer;
    typedef FPED::MemberConstIterator  member_const_iterator;
    typedef FPED::Pedigree             pedigree_type;
    typedef FPED::PedigreePointer      pedigree_pointer;
    typedef FPED::PedigreeConstPointer pedigree_const_pointer;

    typedef polygenic_penetrance_calculator::penetrance_info   penetrance_info;

    typedef peeling::peeler<genetic_info,log_double,
                      peeling::individual_cache<genetic_info,
                                                log_double> >  peeler_type;

    struct i_info
    {
      i_info();
      i_info(const i_info&);
      log_double& operator() (int, int);
      log_double data[3][MAX_POLYGENOTYPE];
    };

    struct i_s_info
    {
      i_s_info(member_const_pointer id);
      i_s_info(const i_s_info&);
      log_double& operator() (int, int, int, int);
      log_double data[3][MAX_POLYGENOTYPE][3][MAX_POLYGENOTYPE];
      member_const_pointer     spouse_id;
    };

    struct i_m_f_info
    {
      i_m_f_info();
      i_m_f_info(const i_m_f_info&);
      log_double& operator()(int i, int j, int k, int l, int m, int n);
      log_double data[3][MAX_POLYGENOTYPE][3][MAX_POLYGENOTYPE][3][MAX_POLYGENOTYPE];
    };

    FPMM_SL(const FPED::Multipedigree & ped_data,
            const model               & modp,
            bool                        use_ascertainmentp);
//    FPMM_SL(cerrorstream & errors = sage_cerr);
//    FPMM_SL(const FPMM_SL &);
//    FPMM_SL & operator=(const FPMM_SL &);

    void set_peeler(peeler_type* peeler_pointer_param)
    {
      peeler_pointer = peeler_pointer_param;
      size_t num_of_inds = peeler_pointer->get_subpedigree().member_count();

      my_founder_SL2_4    .clear();
      my_founder_SL6      .clear();
      my_nonfounder_SL2_4 .clear();
      my_nonfounder_SL6   .clear();
      my_founder_SL2_4    .resize(num_of_inds);
      my_founder_SL6      .resize(num_of_inds);
      my_nonfounder_SL2_4 .resize(num_of_inds);
      my_nonfounder_SL6   .resize(num_of_inds);
    }
    peeler_type* get_peeler() { return peeler_pointer; }

    log_double     founder_SL2    (const penetrance_info&          ps_indiv) const;   
    log_double     founder_SL4    (const penetrance_info&          ps_indiv) const;
    log_double     founder_SL6    (const penetrance_info&          ps_indiv, 
                                   const penetrance_info&          ps_spouse) const;
    log_double     nonfounder_SL2 (const penetrance_info&          ps_indiv, 
                                   const penetrance_info&          ps_mother,
                                   const penetrance_info&          ps_father) const;
    log_double     nonfounder_SL4 (const penetrance_info&          ps_indiv, 
                                   const penetrance_info&          ps_mother, 
                                   const penetrance_info&          ps_father) const;
    log_double     nonfounder_SL6 (const penetrance_info&          ps_indiv, 
                                   const penetrance_info&          ps_spouse) const;

    int update() { return pen.update(); }

    // -------------------------------------------------
    // unconnected likelihood for unconnected individual
    // -------------------------------------------------

    log_double unconnected_likelihood(member_const_pointer memit,
				      genotype_index       genotype,
                                      int                  polygenotype)
    {
	log_double p(1.0);
        
        // Get the population frequencies of the genotype and polygenotype

	double psi = mod.freq_sub_model.prob(genotype);
        double omg = mod.fpmm_sub_model.pop_freq(polygenotype);

        // Get the penetrance.

        p = psi * omg * penetrance(memit,genotype,polygenotype);

	return p;
    }

    double penetrance(member_const_pointer memit, genotype_index genotype, int polygenotype)
    {
        penetrance_info temp;
        temp.member  = memit;
        temp.genotype=genotype;
        temp.polygenotype=polygenotype;
	return pen.get_polygenic_penetrance(temp);
    }

  protected:
    peeler_type * peeler_pointer;

  private:
    log_double int_founder_SL2_4    (const penetrance_info&          ps_indiv) const;   
    log_double int_founder_SL6      (const penetrance_info&          ps_indiv, 
                                     const penetrance_info&          ps_spouse) const;

    log_double int_nonfounder_SL2_4 (const penetrance_info&          ps_indiv, 
                                     const penetrance_info&          ps_mother,
                                     const penetrance_info&          ps_father) const;
    log_double int_nonfounder_SL6   (const penetrance_info&          ps_indiv, 
                                     const penetrance_info&          ps_spouse) const;

    mutable polygenic_transition_calculator trans_poly_calc;

    mutable vector<i_info>                my_founder_SL2_4;
    mutable vector<vector<i_s_info> >     my_founder_SL6;
    mutable vector<i_m_f_info>            my_nonfounder_SL2_4;
    mutable vector<vector<i_s_info> >     my_nonfounder_SL6;

    polygenic_penetrance_calculator    pen;
    const model &                       mod;
};


}} // End namespace

#include "segreg/SL_calculator.ipp"

#endif
