//===================================================================
//
//  File:	FPMM_peeler.cpp
//
//  Author:	Stephen Gross
//
//  History:    sag Initial implementation		Jul 26 01
//
//  Copyright (c) 2001, R. C. Elston
//  All rights reserved
//===================================================================

#include "mped/mp_utilities.h"
#include "segreg/FPMM_peeler.h"

using namespace std;

namespace SAGE {
namespace SEGREG {

//===================================================================
//
//  FPMM_peeler(...) constructor #2
//
//===================================================================
FPMM_peeler::FPMM_peeler(const subped_type& subped, const model & modp)
  : peeling::peeler<genetic_info,log_double,
    peeling::individual_cache<genetic_info,log_double> > (subped) 
{
  my_model = &modp;
  tabs     = 0;
  thresh   = 3;
}

//===================================================================
//
//  partial_parental_likelihood  Sis(Ui,Vi)
//
//===================================================================
log_double
FPMM_peeler::partial_parental_likelihood(
	const member_type& indiv, genotype_index genotype_indiv, 
        size_t polygenotype_indiv, const member_type& spouse)
{
//  t1("ppl1");

  genetic_info indiv_data(genotype_indiv,polygenotype_indiv);

  log_double anterior1  = anterior(indiv,indiv_data);
  log_double posterior1 = posterior_except_mate(indiv,spouse,indiv_data);
  log_double ppl        = anterior1 * posterior1;

//  t2("ppl1");
  return ppl;
}

//===================================================================
//
//  partial_parental_likelihood  Sis(Ui,Vi|Us,Vs)
//
//===================================================================
log_double
FPMM_peeler::partial_parental_likelihood(
  const member_type& indiv,  genotype_index genotype_indiv,  size_t polygenotype_indiv,
  const member_type& spouse, genotype_index genotype_spouse, size_t polygenotype_spouse)
{
//  t1("ppl2");

  genetic_info indiv_data (genotype_indiv, polygenotype_indiv );
  genetic_info spouse_data(genotype_spouse,polygenotype_spouse);

  log_double anterior1  = internal_anterior_with_mate(indiv,spouse,indiv_data,spouse_data,anterior1);
  log_double posterior1 = posterior_except_mate      (indiv,spouse,indiv_data                      );
  log_double ppl        = anterior1 * posterior1;

//  t2("ppl2");
  return ppl;
}

//===================================================================
//
//  internal_anterior(...)
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_anterior(const member_type& ind, const data_type & g, result_type & ia)
{
//  t1("int ant");
  // Establish return value:
  log_double anterior(0.0);

  // Create penetrance_info's to pass on to SL equations:

  penetrance_info indiv (ind,g.genotype,g.polygenotype);

  penetrance_info mother(*ind.get_mother());
  penetrance_info father(*ind.get_father());

  // Start calculating the anterior:

  log_double mother_genotype_loop_total     (0.0),
             mother_polygenotype_loop_total (0.0),
             father_genotype_loop_total     (0.0),          
             father_polygenotype_loop_total (0.0),
             pplikelihood1                  (0.0),
             pplikelihood2                  (0.0),
             sibling_likelihood             (0.0);

  for(mother.genotype = index_AA; mother.genotype != index_INVALID; ++mother.genotype)
  {
    mother_polygenotype_loop_total = 0.0;
    for(mother.polygenotype = 0; mother.polygenotype < my_model->fpmm_sub_model.max_pgt(); ++mother.polygenotype)
    {
      father_genotype_loop_total = 0.0;
      for(father.genotype = index_AA; father.genotype != index_INVALID; ++father.genotype)
      {
        father_polygenotype_loop_total = 0.0;
        for(father.polygenotype = 0; father.polygenotype < my_model->fpmm_sub_model.max_pgt(); ++father.polygenotype)
        {
          sibling_likelihood  = SL_calc->nonfounder_SL2(indiv, mother, father);
          pplikelihood2       = partial_parental_likelihood(
                                *father.member,father.genotype,father.polygenotype,
                                *mother.member,mother.genotype,mother.polygenotype);

          father_polygenotype_loop_total += pplikelihood2 * sibling_likelihood;
        }
        father_genotype_loop_total += father_polygenotype_loop_total;
      }
      pplikelihood1 = partial_parental_likelihood(*mother.member,mother.genotype,mother.polygenotype,*father.member);

      mother_polygenotype_loop_total += pplikelihood1 * father_genotype_loop_total;
    }
    mother_genotype_loop_total += mother_polygenotype_loop_total;
  }

  anterior = mother_genotype_loop_total;

  if(SAGE::isnan(anterior.get_double())) ia = 0.0;
  else                             ia = anterior;

  // Return the calculated anterior

//  t2("int ant");
  return ia;
}

//===================================================================
//
//  internal anterior with mate (ant*)
//
//===================================================================
const FPMM_peeler::result_type &
FPMM_peeler::internal_anterior_with_mate(const member_type& ind, const member_type& spouse, 
  const data_type & g, const data_type & h, result_type & iawm)
{
  return iawm = anterior(ind, g);
}

//===================================================================
//
//  internal_anterior_terminal(...)
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_anterior_terminal(const member_type& ind,const data_type & g, 
                                              result_type & iatg)
{
  t1("int ant term");
  penetrance_info indiv(ind,g.genotype,g.polygenotype);

  iatg = SL_calc->founder_SL2(indiv);

  if(SAGE::isnan(iatg.get_double())) iatg = 0.0;

  t2("int ant term");
  return iatg;
}

//===================================================================
//
//  internal_posterior(...)
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_posterior(const member_type& ind, const data_type & g, 
                                      result_type & ip)
{
  t1("int post");
  log_double posterior_loop_total(1.0);

  // Mate loop begins here

  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
    posterior_loop_total *= posterior_with_mate(ind,mate_loop->mate(),g);

  if(SAGE::isnan(posterior_loop_total.get_double())) ip = 0.0;
  else                                         ip = posterior_loop_total;

  t2("int post");
  return ip;
}

//===================================================================
//
//  internal_posterior_except_mate(...)
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_posterior_except_mate(const member_type& ind, const member_type& mate,
                                                  const data_type & g,
                                                  result_type & ipem)
{
  t1("int post exc mate");
  log_double posterior_loop_total(1.0);

  // Mate loop begins here

  member_type::mate_const_iterator mate_loop = ind.mate_begin();

  for( ; mate_loop != ind.mate_end(); ++mate_loop)
    if(&mate_loop->mate() != &mate)
      posterior_loop_total *= posterior_with_mate(ind,mate_loop->mate(),g);

  if(SAGE::isnan(posterior_loop_total.get_double())) ipem = 0.0;
  else                                         ipem = posterior_loop_total;

  t2("int post exc mate");
  return ipem;
}

//===================================================================
//
//  internal_posterior_with_mate(...)
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_posterior_with_mate(
        const member_type&            ind, const member_type&        mate,
        const data_type & g,   result_type & ipwm)
{
  t1("int post with mate");
  log_double posterior(0.0);
  log_double ppl1(0.0);
  log_double sibling_likelihood(0.0);
  log_double genotype_spouse_total(0.0);
  log_double polygenotype_spouse_total(0.0);

  penetrance_info indiv_data(ind,g.genotype,g.polygenotype);
  penetrance_info spouse_data(mate);

  // Genotype spouse loop begins here

  for(int genotype_spouse_loop = 0; genotype_spouse_loop < 3; ++genotype_spouse_loop)
  {
    polygenotype_spouse_total = 0;
    spouse_data.genotype      = (genotype_index)genotype_spouse_loop;

    // Polygenotype spouse loop begins here
    for(spouse_data.polygenotype = 0; 
        spouse_data.polygenotype < my_model->fpmm_sub_model.max_pgt(); 
      ++spouse_data.polygenotype)
    {
      ppl1               = partial_parental_likelihood(
                           mate,             spouse_data.genotype,spouse_data.polygenotype,
                           *indiv_data.member,indiv_data .genotype,indiv_data .polygenotype);

      sibling_likelihood = SL_calc->nonfounder_SL6(indiv_data,spouse_data);

      polygenotype_spouse_total += ppl1 * sibling_likelihood;

    } // Polygenotype spouse loop ends here

    genotype_spouse_total += polygenotype_spouse_total;

  } // Genotype spouse loop ends here

  posterior = genotype_spouse_total;

  if(SAGE::isnan(posterior.get_double())) ipwm = 0.0;
  else                              ipwm = posterior;

  t2("int post with mate");
  return ipwm;
}

//===================================================================
//
//  internal_posterior_terminal();
//
//===================================================================
const FPMM_peeler::result_type & 
FPMM_peeler::internal_posterior_terminal(const member_type& ind,
                                               const data_type & g,
                                               result_type & ipt)
{
  t1("int post term");
  ipt = 1.0;
  t2("int post term");
  return ipt;
}

}
}
