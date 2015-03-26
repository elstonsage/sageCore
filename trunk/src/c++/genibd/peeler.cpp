//============================================================================
// File:      peeler.cpp
//
// Author:    Dan Baechle
//
// History:   10/17/2 - created.                                 djb
//            Modified for GENIBD                                yjs Apr 2004
//
// Notes:     Lodlink peeler implementation.
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

#include "genibd/peeler.h"

namespace SAGE
{

namespace GENIBD
{

//============================================================================
// IMPLEMENTATION:  peeler
//============================================================================
//
peeler::peeler(const subped_type& subped, const inheritance_model& model)
      : peeling::peeler<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>(subped), my_model(model)
{
  fmember_const_iterator  ind_iter = my_subpedigree.member_begin();

  for( ; ind_iter != my_subpedigree.member_end(); ++ind_iter )
  {
    size_t  mph_id = ind_iter->subindex() + 1;
    assert(mph_id != MLOCUS::NPOS);

    //cout << "building '" << ind_iter->name()
    //     << "', mph_id = '" << mph_id
    //     << "', model = " << my_model.name() << endl;

    my_cache.get_individual_cache(*ind_iter).build(*ind_iter, mph_id, my_model);
  }
}                       

// - Equation (3).  Fernando, Stricker, Elston.  1993.
//
// - Base peeler insures that individual is not an anterior terminal.
//
const log_double&  
peeler::internal_anterior(const member_type& ind, 
                          const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_anterior()");
#endif

  log_double p = internal_anterior_phased(ind, ipi, false);

#ifdef DEBUG_PEELER
  print_exit(mother_sum.get_double(), "internal_anterior()");
#endif  
  
  return result = p;
}

const log_double&
peeler::internal_posterior(const member_type& ind,
                           const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_posterior()");
#endif

  log_double  r(1);
  fmate_const_iterator  iter;
  for( iter = ind.mate_begin(); iter != ind.mate_end(); ++iter )
  {
    r *= posterior_with_mate(ind, iter->mate(), ipi);
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
                                     const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter_w_mate(ind, mate, ipi, "internal_posterior_with_mate()");
#endif

  ind_genotype ind_geno(ipi);

  size_t mate_index = mate.subindex() + 1;
  
  log_double  mate_sum(0);

  for(MLOCUS::penetrance_model::phased_penetrance_iterator mate_iter =
        my_model.phased_penetrance_begin(mate_index);
      mate_iter != my_model.phased_penetrance_end(mate_index);
      ++mate_iter )
  {
    // If x-linked and mate is the male, skip non-Y genotypes
    if(    my_model.is_x_linked()
        && mate.is_male()
        && !is_Y_genotype(mate_iter.phased_geno()) )
      continue;

    ind_genotype mate_geno(mate_iter);

    log_double total_offspring = likelihood_of_offspring(ind_geno, mate_geno, ind, mate);

    if( total_offspring.get_double() != 0. )
    {
      log_double  mate_term(1);

      mate_term *= anterior(mate, mate_iter);
      mate_term *= *mate_iter;
      mate_term *= posterior_except_mate(mate, ind, mate_iter);  
      mate_term *= total_offspring;

      mate_sum += mate_term; 
    }
  }
  
#ifdef DEBUG_PEELER
  print_exit(mate_sum.get_double(), "internal_posterior_with_mate()");
#endif
  
  return result = mate_sum;
}

const log_double&
peeler::internal_posterior_except_mate(const member_type& ind, const member_type& mate,
                                       const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_posterior_except_mate()");
#endif

  log_double  r(1);
  fmate_const_iterator  iter;
  for( iter = ind.mate_begin(); iter != ind.mate_end(); ++iter )
  {
    if( iter->mate().index() != mate.index() )
    {
      r *= posterior_with_mate(ind, iter->mate(), ipi);
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
                                    const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_posterior_terminal()");
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
                                   const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, log_double& result)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_anterior_terminal()");
#endif

  log_double  r = log_double(ind_genotype(ipi).frequency());

#ifdef DEBUG_PEELER
  print_exit(r.get_double(), "internal_anterior_terminal()");
#endif

  return result = r;
}


//
// ----------------------------  ancillary functions  -------------------------
//

log_double
peeler::internal_anterior_phased(const member_type& ind, 
                                 const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, bool flip_allele)
{
#ifdef DEBUG_PEELER
  print_enter(ind, ipi, "internal_anterior_phased()");
#endif

  FPED::MemberConstPointer  mother = ind.parent1();
  FPED::MemberConstPointer  father = ind.parent2();
  
  if(mother->is_male() || father->is_female())
  {
    std::swap(mother, father);
  }
 
  ind_genotype  ind_geno(ipi);
 
  size_t mother_index = mother->subindex() + 1;
  size_t father_index = father->subindex() + 1;
  
  log_double  mother_sum(0);

  for(MLOCUS::penetrance_model::phased_penetrance_iterator moth_pen_iter
        = my_model.phased_penetrance_begin(mother_index);
      moth_pen_iter != my_model.phased_penetrance_end(mother_index);
      ++moth_pen_iter )
  {
    ind_genotype mi_geno(moth_pen_iter);

    log_double  father_sum(0);

    for(MLOCUS::penetrance_model::phased_penetrance_iterator fath_pen_iter
          = my_model.phased_penetrance_begin(father_index);
        fath_pen_iter != my_model.phased_penetrance_end(father_index);
        ++fath_pen_iter )
    {
      // If x-linked, skip father's non Y-chromosomal genotypes
      if(    my_model.is_x_linked()
          && !is_Y_genotype(fath_pen_iter.phased_geno()) )
        continue;

      ind_genotype fi_geno(fath_pen_iter);

      // - Check value of transition first.  This improves performance as 
      //   there will be many instances where the transition probability is 0.
      //
      double  trans = my_tcalc.transition(mi_geno, fi_geno, ind_geno);

      if( trans == 0. )
      {
        continue;
      }

      log_double  total_sib = likelihood_of_sibs(mi_geno, fi_geno, ind);

      if( total_sib.get_double() != 0. )
      {
        log_double  father_term(1);

        father_term *= anterior(*father, fath_pen_iter);
        father_term *= *fath_pen_iter;
        father_term *= log_double(trans);
        father_term *= posterior_except_mate(*father, *mother, fath_pen_iter);
        father_term *= total_sib;

        father_sum += father_term;
      }
    }

    if( father_sum.get_double() != 0. )
    {
      log_double  mother_term(1);

      mother_term *= anterior(*mother, moth_pen_iter);
      mother_term *= *moth_pen_iter;
      mother_term *= posterior_except_mate(*mother, *father, moth_pen_iter);    
      mother_term *= father_sum;

      mother_sum += mother_term;
    }
  }
  
#ifdef DEBUG_PEELER
  print_exit(mother_sum.get_double(), "internal_anterior_phased()");
#endif  
  
  return mother_sum;
}


// - Likelihood of offspring of a given individual w. a particular mate.
//   Line 2 of equation (4) of Fernando, Stricker, Elston.  1993.
//
log_double
peeler::likelihood_of_offspring(const ind_genotype& ind_g, const ind_genotype& mate_g,
                                const member_type&  ind,   const member_type&  mate)
{
  log_double  r(1);
  
  foffspring_const_iterator iter = ind.offspring_begin(mate);
  for( ; iter != ind.offspring_end(); ++iter )
  {
    if( ind.is_male() )
    {
      r *= sum_ind_and_posterior(mate_g, ind_g, *iter);
    }
    else
    {
      r *= sum_ind_and_posterior(ind_g, mate_g, *iter);
    }
  }
  
  return r;
}

// - Likelihood of sibs of given individual.  Line 4 of equation (3) of Fernando, Stricker,
//   Elston.  1993.
//
// Calcuates the posterior probability of the sibs of an individual,
//  given the genotypes of the parents.
//
log_double
peeler::likelihood_of_sibs(const ind_genotype& mg, const ind_genotype& fg,
                           const member_type& ind)
{
  log_double r(1);
  
  fsibling_const_iterator iter = ind.sibling_begin();
  for( ; iter != ind.sibling_end(); ++iter )
  {
    r *= sum_ind_and_posterior(mg, fg, *iter);
  }
  
  return r;
}

// Calculates the same as sib_prob, except it ignores the sibling ex_sib
//  (this is used for calculating all sibs except those we've conditioned
//  for, in sibling, avuncular, grandparental, and cousin pairs
//
log_double
peeler::likelihood_of_sibs_except_sib(const ind_genotype& mg, const ind_genotype& fg,
                                      const member_type& ind, const member_type& ex_sib)
{
  log_double r(1);
  
  fsibling_const_iterator iter = ind.sibling_begin();
  for( ; iter != ind.sibling_end(); ++iter )
  {
    if( iter->name() == ex_sib.name() )
      continue;

    r *= sum_ind_and_posterior(mg, fg, *iter);
  }
  
  return r;
}

//
// -----------------------  debugging functions  -------------------------
//
void
peeler::print_enter(const member_type& ind, 
                    const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, const string& func_name)
{
  cout << "ENTERING " << func_name << " ...\n"
       << "individual " << ind.name() 
       << ", ";
       
  ind_genotype  jg(ipi);
  jg.print();
       
  cout << "\n" << endl;  
}
    
void
peeler::print_enter_w_mate(const member_type& ind, const member_type& mate,
                   const MLOCUS::penetrance_model::phased_penetrance_iterator& ipi, const string& func_name) 
{
  cout << "ENTERING " << func_name << " ...\n"
       << "individual " << ind.name()
       << ", mate " << mate.name() 
       << ", ";
       
  ind_genotype  jg(ipi);
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

} // end of namespace GENIBD

} // end of namespace SAGE


