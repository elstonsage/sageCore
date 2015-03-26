//===========================================================================
//  File:       mlm_peeler.ipp
//
//  Author:     Kai He                          
//
//  History:
//                                                 
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===========================================================================

#ifndef MLM_PEELER_H
#include "segreg/mlm_peeler.h"
#endif

namespace SAGE   { 
namespace SEGREG {

inline
  MlmLikelihoodElements::MlmLikelihoodElements
    (const FPED::Multipedigree& ped_data,
     const model&               mdl,
     bool                       use_ascertainment)
  : my_freq_function(boost::bind(&genotype_frequency_sub_model::prob,
                                 boost::ref(mdl.freq_sub_model), _1)),
    my_transm_function(boost::bind(&transmission_sub_model::prob,
                                 boost::ref(mdl.transm_sub_model), _1, _2, _3))
{
  build_type_description();
}

inline int
  MlmLikelihoodElements::update()
{
//  int err;

//  err = my_pen_calc.update();  if(err) return err;

  return 0;
}

inline const TypeDescription&
  MlmLikelihoodElements::get_type_description() const
{
  return my_type_description;
}

inline double
  MlmLikelihoodElements::get_frequency
    (const TypeDescription::State& state) const
{
  return my_freq_function(state.get_index());
}

inline double
  MlmLikelihoodElements::get_frequency
    (const TypeDescription::StateIterator& state) const
{
  return get_frequency(*state);
}

inline double 
  MlmLikelihoodElements::get_transmission
    (const TypeDescription::State& istate,
     const TypeDescription::State& mstate,
     const TypeDescription::State& fstate)    const
{
  return my_transm_function(istate.get_index(), mstate.get_index(), fstate.get_index());
}

inline double 
  MlmLikelihoodElements::get_transmission
    (const TypeDescription::StateIterator& istate,
     const TypeDescription::StateIterator& mstate,
     const TypeDescription::StateIterator& fstate)    const
{
  return get_transmission(*istate, *mstate, *fstate);
}

inline double
  MlmLikelihoodElements::get_penetrance
    (const FPED::Member&                ind,
     const TypeDescription::State&      state) const
{
//  return my_pen_calc.get_penetrance(
//            RegPenetranceCalculator::penetrance_info(ind, state.get_index()));
  return 0.0;
}

inline double
  MlmLikelihoodElements::get_penetrance
    (const FPED::Member&                        ind,
     const TypeDescription::StateIterator&      state) const
{
//  return get_penetrance(ind, *state);
  return 0.0;
}

inline void
  MlmLikelihoodElements::build_type_description()
{
  my_type_description.set_name("Main Type");
  
  my_type_description.add_state("AA");
  my_type_description.add_state("AB");
  my_type_description.add_state("BB");
}
  
//===========================================================================
//
//  mlm_peeler(...) CONSTRUCTOR
//
//===========================================================================
inline
mlm_peeler::mlm_peeler
    (const subped_type&           subped,
     binary_member_calculator*    bmc, 
     const FPED::Multipedigree&   ped_data,
     const model&                 md,
     const MlmLikelihoodElements& lelt,
     bool                         use_ascertainmentp)
  :  peeling::peeler<TypeDescription::State, log_double,
            peeling::individual_cache<TypeDescription::State, log_double> >(subped),
     mod(md),
     mcc(bmc),
     my_lelements(lelt),
     my_fam_resid_calc
        (ResidualGetter(md.resid_sub_model),
         boost::bind(&binary_member_calculator::get_aff_status,
                     boost::ref(mcc), _1),
         boost::bind(&binary_member_calculator::get_expected_susc,
                     boost::ref(mcc), _1, _2),
         boost::bind(&transmission_sub_model::prob, boost::ref(md.transm_sub_model),
                     _4, _2, _3),
         boost::counting_iterator<int>(index_AA),
         boost::counting_iterator<int>(index_INVALID))
{
  // Set each individual cach

  subped_type::member_const_iterator mem = subped.member_begin();

  for( ; mem != subped.member_end(); ++mem)
  {
    my_cache.get_individual_cache(*mem).set_member(*mem);
  }
}

/// Calculates \f$\rho\f$ (Equ. 78, SEGREG Formula Documentation)
inline double 
mlm_peeler::family_rho
    (const member_type& i,
     const member_type& s,
     const data_type&   p1,
     const data_type&   p2)
{
  double  P   = 0.0;

  member_const_pointer mother, father;
  genotype_index m_geno, f_geno;

  // Sort the parents into mother and father

  if(i.is_female())
  {
     mother = &i;
     father = &s;

     m_geno = p1.get_index();
     f_geno = p2.get_index();
  }
  else
  {
     mother = &s;
     father = &i;

     m_geno = p2.get_index();
     f_geno = p1.get_index();
  }

  // Get pedigree statistics

  family_const_pointer   fam = my_subpedigree.pedigree()->family_find(*mother,*father);   

  P = my_fam_resid_calc.calculate_adjustment(*fam, m_geno, f_geno);
  
  return P;
}

//===========================================================================

}}
