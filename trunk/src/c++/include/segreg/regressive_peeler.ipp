//===================================================================
//
//  File:	regressive_peeler.ipp
//
//  Author:	Stephen Gross
//
//===================================================================

#ifndef REGRESSIVE_PEELER_H
#include "segreg/regressive_peeler.h"
#endif

namespace SAGE { namespace SEGREG {

//===================================================================
//
//  regressive_peeler(...)
//
//===================================================================
inline
regressive_peeler::regressive_peeler
  (const subped_type&             subped,
   const LikelihoodElements&      lelt)
  : peeling::peeler<TypeDescription::State,log_double,
    peeling::individual_cache<TypeDescription::State,log_double> > (subped),
    my_likelihood_elts(lelt)
{
  // Set each individual cache

  subped_type::member_const_iterator mem = subped.member_begin();

  for( ; mem != subped.member_end(); ++mem)
  {
    my_cache.get_individual_cache(*mem).set_member(*mem);
  }
}

inline const regressive_peeler::result_type& 
regressive_peeler::anterior_with_mate(
	const member_type& ind, const member_type& mate, const data_type & g, const data_type & h)
{
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  if(c.anterior_with_mate_cached(mate, g, h)) return c.anterior_with_mate(mate, g, h);
  
  result_type& r = c.anterior_with_mate(mate, g, h);

  internal_anterior_with_mate(ind, mate, g, h, r);
  
  return r;
}

inline const regressive_peeler::result_type&
regressive_peeler::partial_parental_likelihood(
	const member_type& ind, const member_type& mate, const data_type& g)
{
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  if(c.ppl_cached(mate, g)) return c.ppl(mate, g);
  
  result_type& r = c.ppl(mate, g);

  internal_partial_parental_likelihood(ind, mate, g, r);
  
  return r;
}

inline const regressive_peeler::result_type& 
regressive_peeler::partial_parental_likelihood(
	const member_type& ind, const member_type& mate, const data_type& g, const data_type& h)
{
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  if(c.ppl_with_mate_cached(mate, g, h)) return c.ppl_with_mate(mate, g, h);
  
  result_type& r = c.ppl_with_mate(mate, g, h);

  internal_partial_parental_likelihood(ind, mate, g, h, r);
  
  return r;
}

inline log_double 
regressive_peeler::calculate_founder_anterior
    (const member_type&          ind,
     const data_type&            g,
     const PenetranceContext&    context)
{
    // Calculate psi and penetrance values:
    double penetrance = my_likelihood_elts.get_penetrance(ind, g, context);
    double psi        = my_likelihood_elts.get_frequency(g);

    return log_double(psi * penetrance);
}

}}
