//============================================================================
// File:      likelihood.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/2/2 - created.                                djb
//                                                                          
// Notes:     implementation of likelihood calculator classes.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/likelihood.h"

namespace SAGE
{

namespace LODLINK
{

void
check_parameters(const MaxFunction::parameter_vector& theta, const mle_sub_model& mle)
{
  const MAXFUN::ParameterMgr* parameter_manager = mle.parameter_mgr();
  assert(parameter_manager);
  
  size_t  parameter_count = static_cast<size_t>(parameter_manager->getParamCount());
  
  assert(parameter_count <= theta.size());
  for(size_t i = 0; i < parameter_count; ++i)
  {
    assert(theta[i] == (*parameter_manager)(i));
  }
}

//============================================================================
// IMPLEMENTATION:  subped_calculator
//============================================================================
//
log_double
subped_calculator::likelihood()
{
  const FPED::FilteredMultipedigree::subpedigree_type&  subped = my_peeler.subpedigree();
  FPED::FilteredMultipedigree::member_const_iterator  m_iter = subped.member_begin();
  size_t  trait = my_peeler.trait();
  size_t  marker = my_peeler.marker();  

  log_double  like(0);
          
  phenoset  ph_set(subped, trait, marker, *m_iter);
  phenoset::phenoset_iterator  ph_iter = ph_set.begin();
  for(; ph_iter != ph_set.end(); ++ph_iter)
  {
    log_double like_term(1);
    like_term *= my_peeler.anterior(*m_iter, *ph_iter);
    like_term *= (*ph_iter).penetrance();
    like_term *= my_peeler.posterior(*m_iter, *ph_iter);  
    
    like += like_term;
  }
  
  const mle_sub_model&  my_mle = my_peeler.tcalc().mle();
  if(my_mle.uses_alpha())
  {
    double  my_alpha = my_mle.alpha();
    assert(my_alpha != QNAN);
    like = my_alpha * like + (1 - my_alpha) * unlinked_likelihood();
  }

  return like;
}

log_double
subped_calculator::unlinked_likelihood()
{
  return  unlinked_likelihood_cached ? my_unlinked_likelihood : calc_unlinked_likelihood();
}

// - Calculate subpedigree likelihood for the unlinked case.
//
log_double
subped_calculator::calc_unlinked_likelihood()
{
  const mle_sub_model&  linked_mle = my_peeler.tcalc().mle();
  bool  ss = linked_mle.is_sex_specific();
  
  mle_sub_model  unlinked_mle(ss, false);
  
  if(ss)
  {
    unlinked_mle.set_male_theta(NULL_THETA);
    unlinked_mle.set_female_theta(NULL_THETA);  
  }
  else
  {
    unlinked_mle.set_average_theta(NULL_THETA);
  }
  
  peeler  unlinked_peeler(my_peeler.subpedigree(), unlinked_mle,
                           my_peeler.trait(), my_peeler.marker());
  subped_calculator  unlinked_calculator(unlinked_peeler);
  
  my_unlinked_likelihood = unlinked_calculator.likelihood();
  unlinked_likelihood_cached = true;
  
  return  my_unlinked_likelihood;
}


//============================================================================
// IMPLEMENTATION:  ped_calculator
//============================================================================
//
log_double
ped_calculator::likelihood()
{
  log_double  like(1);
  
  subpedigree_const_iterator  subped_iter = my_ped.subpedigree_begin();
  {
    for(; subped_iter != my_ped.subpedigree_end(); ++subped_iter)
    {
      peeler  p(*subped_iter, my_mle, my_trait, my_marker);
      subped_calculator  sp_calc(p);
      like *= sp_calc.likelihood();
    }
  }
  
  return  like;  
}

log_double
ped_calculator::unlinked_likelihood()
{
  return  unlinked_likelihood_cached ? my_unlinked_likelihood : calc_unlinked_likelihood();
}

log_double
ped_calculator::calc_unlinked_likelihood()
{
  log_double  like(1);
  
  subpedigree_const_iterator  subped_iter = my_ped.subpedigree_begin();
  {
    for(; subped_iter != my_ped.subpedigree_end(); ++subped_iter)
    {
      peeler  p(*subped_iter, my_mle, my_trait, my_marker);
      subped_calculator  sp_calc(p);
      like *= sp_calc.unlinked_likelihood();
    }
  }
  
  my_unlinked_likelihood = like;
  unlinked_likelihood_cached = true;
  
  return  like;  
}


//============================================================================
// IMPLEMENTATION:  group_calculator
//============================================================================
//
void
group_calculator::build_group(const group& g, const FPED::FilteredMultipedigree& mped)
{
  group::const_iterator  g_iter     = g.begin();
  group::const_iterator  g_end_iter = g.end();
  for(; g_iter != g_end_iter; ++g_iter)
  {
    pedigree_const_pointer  pptr = mped.pedigree_find(*g_iter);
    
    assert(pptr);
    my_group.insert(pptr);
  }
}

log_double
group_calculator::likelihood()
{
  log_double  like(1);
  
  set<FPED::PedigreeConstPointer>::const_iterator  p_iter = my_group.begin();
  {
    for(; p_iter != my_group.end(); ++p_iter)
    {
      ped_calculator  ped_calc(**p_iter, my_mle, my_trait, my_marker);
      like *= ped_calc.likelihood();
    }
  }
  
  return  like;  
}


//============================================================================
// IMPLEMENTATION:  mped_calculator
//============================================================================
//
log_double
mped_calculator::likelihood()
{
  log_double  like(1);
  
  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    ped_calculator  p_calc(*ped_iter, my_mle, my_trait, my_marker);
    like *= p_calc.likelihood();
  }
  
  return  like;  
}

log_double
mped_calculator::unlinked_likelihood()
{
  return  unlinked_likelihood_cached ? my_unlinked_likelihood : calc_unlinked_likelihood();
}

log_double
mped_calculator::calc_unlinked_likelihood()
{
  log_double  like(1);
  
  pedigree_const_iterator  ped_iter = my_mped.pedigree_begin();
  for(; ped_iter != my_mped.pedigree_end(); ++ped_iter)
  {
    {
      ped_calculator  p_calc(*ped_iter, my_mle, my_trait, my_marker);
      like *= p_calc.unlinked_likelihood();
    }
  }
  
  my_unlinked_likelihood = like;
  unlinked_likelihood_cached = true;
  
  return  like;  
}

}
}






