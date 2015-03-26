#ifndef MLM_PEELER_H
#define MLM_PEELER_H
//===================================================================
//  File:       mlm_peeler.h
//
//  Author:     Kai He
//
//  History:   
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===================================================================


#include <cmath>
#include "numerics/functions.h"
#include "numerics/log_double.h"
#include "mped/sp.h"
#include "mped/mp.h"
#include "mped/mp_utilities.h"
#include "segreg/model.h"
#include "peeling/peeler3.h"
#include "peeling/cache3.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/peeling_caches.h"
#include "segreg/member_calculator.h"
#include "segreg/mlm_fra_util.h"
#include "segreg/types/TypeDescription.h"

/**@class mlm_peeler
 * Implementing MLM model Equations 73 - 84
 * The mlm_peeler calculates the likelihood for the pedigree P and the genotype Ui 
 * for Multivariate logistic model which is used to deal with discrete traits.
 * The model uses current peeling algorithm - anterior and posterior to travasl each
 * family and individual in the pedigree data, calculates L(Ui), L(Ui|Us), P(Ui|Um,Uf) and
 * logistic penetrance function.
 * In order to archive the goal, there are many functions for simplifying computation.
 *
 * \f[Logistic \ form \ of \ Penetrance \ function \ for \ MLM \ model:\f]
 * \f[\frac{e^{\theta_u(i)y_i}}{1 + e^{\theta_u(i)}}\f]
 * \f[Likelihood \ of \ MLM \ model:\f]
 * \f[L(P,\ u_i)=\sum_{\Lambda_0}\{L(F_0)[\prod_{\Omega_1}\sum_{\Lambda_{1j}}\frac{L(F_{1j})}{L(1j)}
 * [...\prod_{\Omega_{kj}}\sum_{\Lambda_{kj}}\frac{L(F_{kj})}{L(kj)}]...]\}\f]
 *
 *
 */

namespace SAGE   { 
namespace SEGREG {

class MlmLikelihoodElements
{
  public:
  
    MlmLikelihoodElements(const FPED::Multipedigree& mped,
                          const model&               mdl,
                          bool                       use_ascertainment);
                       
    int update();
    
    const TypeDescription& get_type_description() const;
    
    double get_frequency    (const TypeDescription::State&         state) const;
    double get_frequency    (const TypeDescription::StateIterator& state) const;
    double get_transmission (const TypeDescription::State&         istate,
                             const TypeDescription::State&         mstate,
                             const TypeDescription::State&         fstate)    const;
    double get_transmission (const TypeDescription::StateIterator& istate,
                             const TypeDescription::StateIterator& mstate,
                             const TypeDescription::StateIterator& fstate)    const;
    
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::State&         state) const;
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::StateIterator& state) const;

/*
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::State&         state,
                          const PenetranceContext&              context)  const;
    double get_penetrance(const FPED::Member&                   ind,
                          const TypeDescription::StateIterator& state,
                          const PenetranceContext&              context)  const;
*/
  
  private:
  
    void build_type_description();
  
    boost::function<double (const genotype_index&)>  my_freq_function;
    boost::function<double (const genotype_index&,
                            const genotype_index&,
                            const genotype_index&)>  my_transm_function;
                            
    TypeDescription         my_type_description;
};


class mlm_peeler : public 
	           peeling::peeler<TypeDescription::State, log_double,
                   peeling::individual_cache<TypeDescription::State, log_double> >
{
  public:
//    typedef regressive_penetrance_calculator::penetrance_info   penetrance_info;

    typedef FPED::PedigreeConstPointer      pedigree_const_pointer;
    typedef FPED::SubpedigreeConstPointer   subpedigree_const_pointer;
    typedef FPED::FamilyConstPointer        family_const_pointer;
    typedef FPED::MemberConstPointer        member_const_pointer;
 
    typedef FPED::PedigreeConstIterator     pedigree_const_iterator;
    typedef FPED::SubpedigreeConstIterator  subpedigree_const_iterator;
    typedef FPED::FamilyConstIterator       family_const_iterator;
    typedef FPED::MemberConstIterator       member_const_iterator;
    typedef FPED::OffspringConstIterator    offspring_const_iterator;
    typedef transmission_sub_model                        tsm;

    // This structure is for Css_sum(...) function which is part of P_factor

    struct Mem_trait
    {
      member_const_pointer member;
      int                  trait;
      double               LP;
    };

    mlm_peeler(const subped_type&           subped,
               binary_member_calculator*    bmc,
	       const FPED::Multipedigree&   rmp,
               const model&                 md,
               const MlmLikelihoodElements& lelt,
               bool                         use_ascertainmentp=false);

    ~mlm_peeler() { }

    //=================================================================================================
    //
    // peeler algorithm
    //
    //=================================================================================================    

    /// These functions are for test_segreg, and are not used in the actual code.
    //@{
//    double get_trait                                 (const member_type& );
//    double get_genotype_specific_susceptibility      (genotype_index     );
//    bool   is_member_valid                           (const member_type& );
//    void   member                                    (                   );
    //@}

  private:
    // Data members

    const model&               mod;
    binary_member_calculator*  mcc;
    
    const MlmLikelihoodElements& my_lelements;
    
    PED_CALC::ApproximateFamResidAdj
        <genotype_index, FPED::Multipedigree,
         boost::counting_iterator<int> >  my_fam_resid_calc;

    double family_rho  (const member_type&, const member_type&, const data_type&, const data_type&);

    //=================================================================================================
    //
    // Internal functions for mlm_peeler algorithm
    //
    //=================================================================================================

    // These functions are derived from peeler and defined for mlm_peeler use

    const result_type& internal_anterior             (const member_type&,         const data_type&, result_type&);
    const result_type& internal_posterior            (const member_type&,         const data_type&, result_type&);
    const result_type& internal_anterior_terminal    (const member_type&,         const data_type&, result_type&);
    const result_type& internal_posterior_terminal   (const member_type&,         const data_type&, result_type&);
    const result_type& internal_posterior_with_mate  (const member_type&, const member_type&, const data_type&, result_type&);
    const result_type& internal_posterior_except_mate(const member_type&, const member_type&, const data_type&, result_type&);

};

}} // end namespaces

#include "segreg/mlm_peeler.ipp"

#endif
