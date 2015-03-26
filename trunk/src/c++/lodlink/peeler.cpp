//============================================================================
// File:      peeler.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/17/2 - created.                                djb
//                                                                          
// Notes:     Lodlink peeler implementation.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/peeler.h"

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  peeler
//============================================================================
//
peeler::peeler(const subped_type& subped, const mle_sub_model& mle,
               size_t trait, size_t marker)
      : peeling::peeler<joint_pen_iter, log_double>(subped), my_tcalc(mle),
        my_trait(trait), my_marker(marker)
{
  typedef FPED::FilteredMultipedigree::member_const_iterator  member_const_iterator;

  const MLOCUS::penetrance_model&  tpm = ge_models::get_model(subped, trait);
  const MLOCUS::penetrance_model&  mpm = ge_models::get_model(subped, marker);
  
  member_const_iterator  ind_iter;
  for(ind_iter = my_subpedigree.member_begin(); 
      ind_iter != my_subpedigree.member_end(); ++ind_iter)
  {
    size_t  tph_id = ind_iter->subindex() + 1;
    assert(tph_id != MLOCUS::NPOS);
    
    size_t  mph_id = ind_iter->subindex() + 1;
    assert(mph_id != MLOCUS::NPOS);
    
    my_cache.get_individual_cache(*ind_iter).build(*ind_iter, tph_id, mph_id, tpm, mpm);
  }
}                       

// - Equation (3).  Fernando, Stricker, Elston.  1993.
//
// - Base peeler insures that individual is not an anterior terminal.
//
const log_double&  
peeler::internal_anterior(const member_type& ind, 
                          const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, jpi, "internal_anterior()");
#endif

  member_const_pointer  mother = ind.get_mother();
  member_const_pointer  father = ind.get_father();
  
  assert(mother && father);
  
  joint_genotype  ind_genotype(jpi);
 
  phenoset  mothers_phenoset(my_subpedigree, my_trait, my_marker, *mother);
  phenoset  fathers_phenoset(my_subpedigree, my_trait, my_marker, *father);
  
  log_double  mother_sum(0);
  phenoset_iterator  m_iter = mothers_phenoset.begin();
  for(; m_iter != mothers_phenoset.end(); ++m_iter)
  {
    log_double  mother_term(1);
    mother_term *= anterior(*mother, *m_iter);
    mother_term *= (*m_iter).penetrance();
    mother_term *= posterior_except_mate(*mother, *father, *m_iter);
  
    log_double  father_sum(0);
    phenoset_iterator  f_iter = fathers_phenoset.begin();
    for(; f_iter != fathers_phenoset.end(); ++f_iter)
    {
      // - Check value of transition first.  This improves performance as 
      //   there will be many instances where the transition probability is 0.
      //
      double  trans = my_tcalc.transition(*m_iter, *f_iter, ind_genotype);
      if(trans == 0)
      {
        continue;
      }
      
      log_double  father_term(1);
      father_term *= anterior(*father, *f_iter);
      father_term *= (*f_iter).penetrance();
      father_term *= posterior_except_mate(*father, *mother, *f_iter);
      father_term *= log_double(trans);
      father_term *= likelihood_of_sibs(*m_iter, *f_iter, ind); 
      
      father_sum += father_term;
    }
    
    mother_term *= father_sum;
    mother_sum += mother_term;
  }
  
#ifdef DEBUG_PEELER
  print_exit(mother_sum.get_double(), "internal_anterior()");
#endif  
  
  return result = mother_sum;
}

const log_double&
peeler::internal_posterior(const member_type& ind,
                           const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, jpi, "internal_posterior()");
#endif

  log_double  r(1);
  mate_const_iterator  iter;
  for(iter = ind.mate_begin(); iter != ind.mate_end(); ++iter)
  {
    r *= posterior_with_mate(ind, iter->mate(), jpi);
  }
  
#ifdef DEBUG_PEELER
  print_exit(r.get_double(), "internal_posterior()");
#endif

  return result = r;
}

// - Equation (4).  Fernando, Stricker, Elston.  1993.
//
const log_double&  
peeler::internal_posterior_with_mate(const member_type& ind, const member_type& mate, 
                                     const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter_w_mate(ind, mate, jpi, "internal_posterior_with_mate()");
#endif

  joint_genotype  ind_genotype(jpi);
  phenoset  mates_phenoset(my_subpedigree, my_trait, my_marker, mate);
  
  log_double  mate_sum(0);
  phenoset_iterator  m_iter = mates_phenoset.begin();
  for(; m_iter != mates_phenoset.end(); ++m_iter)
  {
    log_double  mate_term(1);
    mate_term *= anterior(mate, *m_iter);
    mate_term *= (*m_iter).penetrance();
    mate_term *= posterior_except_mate(mate, ind, *m_iter);  
    mate_term *= likelihood_of_offspring(ind_genotype, *m_iter, ind, mate);
    
    mate_sum += mate_term; 
  }
  
#ifdef DEBUG_PEELER
  print_exit(mate_sum.get_double(), "internal_posterior_with_mate()");
#endif
  
  return result = mate_sum;
}

const log_double&
peeler::internal_posterior_except_mate(const member_type& ind, const member_type& mate,
                            const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, jpi, "internal_posterior_except_mate()");
#endif

  log_double  r(1);
  mate_const_iterator  iter;
  for(iter = ind.mate_begin(); iter != ind.mate_end(); ++iter)
  {
    if(iter->mate().index() != mate.index())
    {
      r *= posterior_with_mate(ind, iter->mate(), jpi);
    }
  }
  
#ifdef DEBUG_PEELER
  print_exit(r.get_double(), "internal_posterior_except_mate()");
#endif

  return result = r;
}

// - Posterior likelihood of a pedigree 'leaf'.
//
const log_double&
peeler::internal_posterior_terminal(const member_type& ind,
                                    const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, jpi, "internal_posterior_terminal()");
#endif

  log_double  r(1);

#ifdef DEBUG_PEELER
  print_exit(r.get_double(), "internal_posterior_terminal()");
#endif

  return result = r;
}

// - Anterior likelihood of a founder.
//
const log_double&  
peeler::internal_anterior_terminal(const member_type& ind, 
                               const joint_pen_iter& jpi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, jpi, "internal_anterior_terminal()");
#endif

  log_double  r = log_double(joint_genotype(jpi).frequency());
  
#ifdef DEBUG_PEELER
  print_exit(r.get_double(), "internal_anterior_terminal()");
#endif

  return result = r;
}


//
// ----------------------------  ancillary functions  -------------------------
//

// - Likelihood of offspring of a given individual w. a particular mate.
//   Line 2 of equation (4) of Fernando, Stricker, Elston.  1993.
//
log_double
peeler::likelihood_of_offspring(const joint_genotype& ind_jg, const joint_genotype& mate_jg,
                                const member_type& ind, const member_type& mate)
{
  log_double  r(1);
  
  offspring_const_iterator  iter;
  for(iter = ind.offspring_begin(mate); iter != ind.offspring_end(); ++iter)
  {
    assert(! ind.is_sex_unknown());
    
    if(ind.is_male())
    {
      r *= sum_ind_and_posterior(mate_jg, ind_jg, *iter);
    }
    else
    {
      r *= sum_ind_and_posterior(ind_jg, mate_jg, *iter);
    }
  }
  
  return r;
}

// - Likelihood of sibs of given individual.  Line 4 of equation (3) of Fernando, Stricker,
//   Elston.  1993.
//
log_double
peeler::likelihood_of_sibs(const joint_genotype& mjg, const joint_genotype& fjg,
                           const member_type& ind)
{
  log_double  r(1);
  
  sibling_const_iterator  iter;
  for(iter = ind.sibling_begin(); iter != ind.sibling_end(); ++iter)
  {
    r *= sum_ind_and_posterior(mjg, fjg, *iter);
  }
  
  return r;
}

//
// -----------------------  debugging functions  -------------------------
//
void
peeler::print_enter(const member_type& ind, 
                    const joint_pen_iter& jpi, const string& func_name)
{
  cout << "ENTERING " << func_name << " ...\n"
       << "individual " << ind.name() 
       << ", ";
       
  joint_genotype  jg(jpi);
  jg.print();
       
  cout << "\n" << endl;  
}
    
void
peeler::print_enter_w_mate(const member_type& ind, const member_type& mate,
                   const joint_pen_iter& jpi, const string& func_name) 
{
  cout << "ENTERING " << func_name << " ...\n"
       << "individual " << ind.name()
       << ", mate " << mate.name() 
       << ", ";
       
  joint_genotype  jg(jpi);
  jg.print();
       
  cout << "\n" << endl;
}

void
peeler::print_exit(double result, const string& func_name)
{
  cout << "EXITING " << func_name << " ...\n"
       << "return value " << result
       << "\n" << endl;
}

}
}


