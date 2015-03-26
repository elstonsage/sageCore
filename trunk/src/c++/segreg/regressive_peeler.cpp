//===================================================================
//
//  File:	regressive_peeler.cpp
//
//  Author:	Stephen Gross
//
//  History:    sag Initial implementation		Jul 26 01
//
//  Copyright (c) 2001, R. C. Elston
//  All rights reserved
//===================================================================

#include "mped/mp_utilities.h"
#include "segreg/regressive_peeler.h"

using namespace std;

namespace SAGE
{
namespace SEGREG
{

void bombout(string& func, regressive_peeler::result_type & x)
{
  if(x.get_double() == 0)
  {
    cout << func << " bombs at zero " << endl;
    exit(0);
  }
}
//===================================================================
//
//  partial_parental_likelihood  Sis(Ui)
//
//===================================================================

/** This is Eq. 11, modified for recursion.
  *
  * \f[ S_{is}(u_i) = ant_i(u_i) \prod_{j \in S_i,j \not = s} pos_{ij}(u_i)
  *   \f]
  */
const regressive_peeler::result_type&
regressive_peeler::internal_partial_parental_likelihood(
	const member_type& ind, const member_type& mate, const data_type& g, result_type& ippl1)
{
  // Eq. 15b

  result_type anterior1  = anterior(ind,g);
  result_type posterior1 = posterior_except_mate(ind,mate,g);
  ippl1                  = anterior1 * posterior1;
  return ippl1;
}

//===================================================================
//
//  partial_parental_likelihood  Sis(Ui|Us)
//
//===================================================================

/**
  * This is Eq. 15, modified for recursion.
  *
  * \f[
  *   S_{is}(u_i|u_s) = ant{^*}{_i(u_i|u_s)} \prod_{j \in S_i,j \not = s} pos_{ij}(u_i)
  * \f]
  */
const regressive_peeler::result_type&
regressive_peeler::internal_partial_parental_likelihood(
	const member_type& ind, const member_type& mate, const data_type& g, const data_type& h, result_type& ippl2)
{
  // Eq. 15a

  log_double anterior1  = anterior_with_mate(ind,mate,g,h);
  log_double posterior1 = posterior_except_mate(ind,mate,g);
  log_double ppl        = anterior1 * posterior1;
  ippl2                 = ppl;

  return ippl2;
}

//===================================================================
//
//  internal_anterior(...)
//
//===================================================================

/**
  * Calculates Eq. 13 for non-founders:
  *
  * \f[ ant_i(u_i) = 
  *       \sum_{u_m} \left\{  S_{mf}(u_m) 
  *       \sum_{u_f} \left[ (S_{fm}(u_f | u_m) SL1(i) ) \right] \right\} 
  * \f]
  */
const regressive_peeler::result_type & 
regressive_peeler::internal_anterior(const member_type& ind, const data_type & g, result_type & ia)
{
  PenetranceContext context = my_likelihood_elts.get_penetrance_context();
  
  context.set_nuclear_family(*ind.family());
  
  log_double anterior = calculate_nonfounder_anterior(ind, g, context);

  if(SAGE::isnan(anterior.get_double())) ia = 0.0;
  else                             ia = anterior;

  return ia;
}


//===================================================================
//
//  internal anterior with mate (ant * equation)
//
//===================================================================

/**
  * Calculates Eq. 17 for non-founders:
  *
  * \f[ ant{^*}{_i}(u_i|u_s) = 
  *       \sum_{u_m} \left\{ S_{mf}(u_m) \sum_{u_f}
  *       \left[ S_{fm}(u_f|u_m) SL3(i) \right] \right\}
  * \f]
  */
const regressive_peeler::result_type &
regressive_peeler::internal_anterior_with_mate(const member_type& ind, const member_type& spouse, 
  const data_type & g, const data_type & h, result_type & iawm)
{
  PenetranceContext context = my_likelihood_elts.get_penetrance_context();
  
  context.set_link_couple(ind, spouse);
  context.set_link_spouse_state(h);

  log_double    anterior(0.0);

  if(ind.is_founder())
  {
    anterior = calculate_founder_anterior(ind, g, context);
  }
  else // individual is NOT a founder
  {
    context.set_nuclear_family(*ind.family());
  
    anterior = calculate_nonfounder_anterior(ind, g, context);
  }

  if(SAGE::isnan(anterior.get_double())) iawm = 0.0;
  else                                   iawm = anterior;

  return iawm;
}

//===================================================================
//
//  internal_anterior_terminal(...)
//
//===================================================================
/**
  * Calculates Eq. 13 for founders:
  *
  * \f[ ant_i(u_i) =  SL1(i)  \f]
  */
const regressive_peeler::result_type & 
regressive_peeler::internal_anterior_terminal
    (const member_type& ind,
     const data_type &  g, 
     result_type &      iatg)
{
  PenetranceContext context = my_likelihood_elts.get_penetrance_context();
  
  iatg = calculate_founder_anterior(ind,g,context);

  if(SAGE::isnan(iatg.get_double())) iatg = 0.0;

  return iatg;
}

//===================================================================
//
//  internal_posterior(...)
//
//===================================================================
const regressive_peeler::result_type & 
regressive_peeler::internal_posterior(const member_type& ind, const data_type & g, result_type & ip)
{
  log_double posterior(1.0);

  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
  {
    log_double pwm = posterior_with_mate(ind,mate_loop->mate(),g);

    posterior *= pwm;
  }

  if(SAGE::isnan(posterior.get_double())) ip = 0.0;
  else                              ip = posterior;

  return ip;
}

//===================================================================
//
//  internal_posterior_with_mate(...)
//
//===================================================================
const regressive_peeler::result_type & 
regressive_peeler::internal_posterior_with_mate(
        const member_type& ind, const member_type& mate, const data_type & g, result_type & ipwm)
{
  // Create our context
  PenetranceContext context = my_likelihood_elts.get_penetrance_context();

  FPED::FamilyConstPointer fam =
      ind.pedigree()->family_find(ind, mate);

  context.set_nuclear_family(*fam);
  
  if(ind.is_female())  context.set_mother_state(g);
  else                 context.set_father_state(g);

  // Calculate the liklihood by summing over each possible state of the mate
  
  log_double    posterior(0.0);
  
  const TypeDescription& tdesc = my_likelihood_elts.get_type_description();

  for(TypeDescription::StateIterator mate_geno = tdesc.begin();
      mate_geno != tdesc.end();
      ++mate_geno)
  {
    log_double ppl1      = partial_parental_likelihood(mate,ind,*mate_geno,g);

    if(ind.is_female())  context.set_father_state(*mate_geno);
    else                 context.set_mother_state(*mate_geno);

    log_double sib_like  = sibship_likelihood(context); 

    posterior += ppl1 * sib_like;
  }

  if(SAGE::isnan(posterior.get_double())) ipwm = 0.0;
  else                              ipwm = posterior;

  return ipwm;
}

//===================================================================
//
//  internal_posterior_except_mate(...)
//
//===================================================================
const regressive_peeler::result_type & 
regressive_peeler::internal_posterior_except_mate(
  const member_type& ind, const member_type& mate, const data_type & g, result_type & ipem)
{
  log_double posterior(1.0);

  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
    if(&mate_loop->mate() != &mate) posterior *= posterior_with_mate(ind,mate_loop->mate(),g);

  if(SAGE::isnan(posterior.get_double())) ipem = 0.0;
  else                              ipem = posterior;

  return ipem;
}

//===================================================================
//
//  internal_posterior_terminal();
//
//===================================================================
const regressive_peeler::result_type & 
regressive_peeler::internal_posterior_terminal(const member_type& ind,
                                               const data_type & g,
                                               result_type & ipt)
{
  ipt = 1.0;
  return ipt;
}

log_double
  regressive_peeler::calculate_nonfounder_anterior
    (const member_type&          ind,
     const data_type&            g,
     const PenetranceContext&    context)
{
  PenetranceContext local_context = context;
  
  log_double    mother_total(0.0);

  const TypeDescription& tdesc = my_likelihood_elts.get_type_description();

  for(TypeDescription::StateIterator mother_geno = tdesc.begin();
      mother_geno != tdesc.end(); 
      ++mother_geno)
  {
    log_double mother_partial = partial_parental_likelihood(
                                   context.get_mother(),
                                   context.get_father(),
                                   *mother_geno);
                                   
                             
    log_double father_total(0.0);

    for(TypeDescription::StateIterator father_geno = tdesc.begin();
        father_geno != tdesc.end(); 
        ++father_geno)
    {
      log_double father_partial = partial_parental_likelihood(
                                     context.get_father(),
                                     context.get_mother(),
                                     *father_geno,
                                     *mother_geno);
                                
      local_context.set_mother_state(*mother_geno);
      local_context.set_father_state(*father_geno);
      
      double     ind_likelihood     = nonfounder_likelihood(ind, g, local_context);
      log_double sibling_likelihood = sibship_likelihood_excluding_sib(local_context, ind);
      
      father_total += father_partial * ind_likelihood * sibling_likelihood;
    }
    mother_total += mother_partial * father_total;
  }

  log_double anterior = mother_total;

  return anterior;
}

double
  regressive_peeler::nonfounder_likelihood
    (const member_type&          ind,
     const data_type&            geno,
     const PenetranceContext&    context)
{
    const TypeDescription::State& mstate = context.get_mother_state();
    const TypeDescription::State& fstate = context.get_father_state();

    double trans_geno_indiv = my_likelihood_elts.get_transmission
                                    (geno,mstate,fstate);
              
    double penetrance_indiv = my_likelihood_elts.get_penetrance
                                    (ind, geno, context);

    return trans_geno_indiv * penetrance_indiv;
}

log_double
  regressive_peeler::sibship_likelihood
    (const PenetranceContext& context)
{
  // Create the accumulator for the likelihood to be returned.
  log_double result(1.0);
  
  const FPED::Family& fam = context.get_family();

  // Loop through children, accumulating likelihood
  for(FPED::OffspringConstIterator child = fam.offspring_begin();
      child != fam.offspring_end();
      ++child)
  {
    log_double child_prob(0.0);

    const TypeDescription& tdesc = my_likelihood_elts.get_type_description();

    for(TypeDescription::StateIterator child_geno = tdesc.begin();
        child_geno != tdesc.end(); 
        ++child_geno)
    {
      double child_likelihood = nonfounder_likelihood(*child, *child_geno, context);

      log_double posterior_child = posterior(*child,*child_geno);

      child_prob += child_likelihood * posterior_child;
    }

    result *= child_prob;
  }

  return result;

}

log_double
regressive_peeler::sibship_likelihood_excluding_sib
    (const PenetranceContext& context,
     const member_type&       excluded_sib)
{
  // Sibling loop begins here

  log_double sibling_loop_total(1.0);

  member_type::sibling_const_iterator sib_loop = excluded_sib.sibling_begin();

  for( ; sib_loop != excluded_sib.sibling_end(); ++sib_loop)
  {
    log_double state_loop_total(0);

    const TypeDescription& tdesc = my_likelihood_elts.get_type_description();

    for(TypeDescription::StateIterator state_loop = tdesc.begin();
        state_loop != tdesc.end(); 
        ++state_loop)
    {
      double sib_likelihood = nonfounder_likelihood(*sib_loop, *state_loop, context);

      log_double posterior_sib  = posterior(*sib_loop,*state_loop);

      state_loop_total += sib_likelihood * posterior_sib;
    }
    sibling_loop_total *= state_loop_total;
  }

  log_double SL1 =  sibling_loop_total;

  return SL1;
}

}
}
