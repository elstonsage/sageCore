//===================================================================
//
//  File:	peeling_caches.ipp
//
//  Author:	Stephen Gross
//
//  History:	sag Initial implementation		Jul 30 2001
//
//  Copyright (c) 2001 R. C. Elston292
//  All rights reserved
//
//===================================================================

#ifndef PEELING_CACHES_H
#include "segreg/peeling_caches.h"
#endif

namespace SAGE
{
namespace peeling
{
  
//===================================================================
//  individual_cache()  [ Constructor ]
//===================================================================
inline
individual_cache<SEGREG::TypeDescription::State,log_double>::individual_cache()
{
  double QNAN = std::numeric_limits<double>::quiet_NaN();
  my_anterior[0]  = my_anterior[1]  = my_anterior[2]  = QNAN;
  my_posterior[0] = my_posterior[1] = my_posterior[2] = QNAN;
}

inline void
individual_cache<SEGREG::TypeDescription::State,log_double>::set_member(const member_type& m)
{
  // Determine the mate count

  size_t mate_count = m.mate_count();

  my_mate_data.reserve(mate_count);

  member_type::mate_const_iterator mt = m.mate_begin();
 
  for( ; mt != m.mate_end(); ++mt)
  {
    my_mate_data.push_back(mate_data(&mt->mate()));
  }
}

//===================================================================
//  anterior_cached(...)
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior_cached(
  const data_type & gen_info) const
{
  if(SAGE::isnan(my_anterior[gen_info.get_index()].get_double()))  
    return false;
  else
    return true;
}

//===================================================================
//  anterior_with_mate_cached(...)
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior_with_mate_cached(
  const member_type& mate, const data_type& g, const data_type& h) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return !SAGE::isnan(i->anterior_with_mate[g.get_index()][h.get_index()].get_double());
}

//=================================================================== 
//  posterior_cached(...) 
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_cached(
  const data_type & gen_info) const
{
  if(SAGE::isnan(my_posterior[gen_info.get_index()].get_double()))
    return false;
  else
    return true;
}

//===================================================================
//  posterior_with_mate_cached(...) 
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_with_mate_cached(
  const member_type& mate, const data_type & gen_info) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return !SAGE::isnan(i->posterior_with_mate[gen_info.get_index()].get_double());
}

//===================================================================
//  posterior_except_mate_cached(...) 
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_except_mate_cached(
  const member_type& mate, const data_type & gen_info) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return !SAGE::isnan(i->posterior_except_mate[gen_info.get_index()].get_double());
}

//===================================================================
//  ppl_cached(...) 
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl_cached(const member_type& mate, const data_type& g) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return !SAGE::isnan(i->partial_parental_likelihood[g.get_index()].get_double());
}

//===================================================================
//  ppl_with_mate_cached(...) 
//===================================================================
inline bool
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl_with_mate_cached(
	const member_type& mate, const data_type& g, const data_type& h) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return !SAGE::isnan(i->ppl_with_mate[g.get_index()][h.get_index()].get_double());
}

//===================================================================
//  SET anterior(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior(const data_type & p)
{ return my_anterior[p.get_index()]; }

//===================================================================
//  SET anterior_with_mate(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior_with_mate(
	const member_type& mate, const data_type & g, const data_type & h)
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->anterior_with_mate[g.get_index()][h.get_index()];
}

//===================================================================
//  SET posterior(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior(const data_type & p)
{ return my_posterior[p.get_index()]; }

//===================================================================
//  SET posterior_with_mate(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_with_mate(
  const member_type& mate, const data_type & p)
{ 
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->posterior_with_mate[p.get_index()];
}

//===================================================================
//  SET posterior_except_mate(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_except_mate(
  const member_type& mate, const data_type & p)
{ 
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->posterior_except_mate[p.get_index()];
}

//===================================================================
//  SET ppl(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl(
	const member_type& mate, const data_type& g)
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->partial_parental_likelihood[g.get_index()];
}

//===================================================================
//  SET ppl_with_mate(...) 
//===================================================================
inline individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl_with_mate(
	const member_type& mate, const data_type& g, const data_type& h)
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->ppl_with_mate[g.get_index()][h.get_index()];
}

//===================================================================
//  GET anterior(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior(
  const data_type & p) const
{ return my_anterior[p.get_index()]; }

//===================================================================
//  GET anterior_with_mate(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::anterior_with_mate(
	const member_type& mate, const data_type & g, const data_type & h) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->anterior_with_mate[g.get_index()][h.get_index()];
}

//===================================================================
//  GET posterior(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior(
  const data_type & p) const

{ return my_posterior[p.get_index()]; }

//===================================================================
//  GET posterior_with_mate(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_with_mate(
  const member_type& mate, const data_type & p) const

{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->posterior_with_mate[p.get_index()];
}

//===================================================================
//  GET posterior_except_mate(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::posterior_except_mate(
  const member_type& mate, const data_type & p) const

{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->posterior_except_mate[p.get_index()];
}

//===================================================================
//  GET ppl(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl(
	const member_type& mate, const data_type& g) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->partial_parental_likelihood[g.get_index()];
}

//===================================================================
//  GET ppl_with_mate(...) 
//===================================================================
inline const individual_cache<SEGREG::TypeDescription::State,log_double>::result_type &
individual_cache<SEGREG::TypeDescription::State,log_double>::ppl_with_mate(
	const member_type& mate, const data_type& g, const data_type& h) const
{
  mate_data_const_iterator i = find(mate);

  if(i == my_mate_data.end()) SAGE_internal_error();  // Should never happen

  return i->ppl_with_mate[g.get_index()][h.get_index()];
}

inline individual_cache<SEGREG::TypeDescription::State,log_double>::mate_data_const_iterator
  individual_cache<SEGREG::TypeDescription::State,log_double>::find(const member_type& mate) const
{
  for(mate_data_const_iterator i  = my_mate_data.begin();
                               i != my_mate_data.end();
                             ++i)
  {
    if(i->mate == &mate)
     return i;
  }

  return my_mate_data.end();  
}

inline individual_cache<SEGREG::TypeDescription::State,log_double>::mate_data_iterator
  individual_cache<SEGREG::TypeDescription::State,log_double>::find(const member_type& mate)
{
  for(mate_data_iterator i  = my_mate_data.begin();
                               i != my_mate_data.end();
                             ++i)
  {
    if(i->mate == &mate)
     return i;
  }

  return my_mate_data.end();  
}

//===================================================================
//  polygenotype_mate_info()  [ Constructor ]
//===================================================================
inline
individual_cache<genetic_info,
  log_double>::polygenotype_mate_info::polygenotype_mate_info(
  const member_type* mate_param)
{
  mate = mate_param;
  log_double QNAN = log_double(numeric_limits<double>::quiet_NaN());

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      with_mate[i][j] = except_mate[i][j] = QNAN;
}

//===================================================================
//  polygenotype_mate_info(...) COPY constructor
//===================================================================
inline
individual_cache<genetic_info,log_double>::polygenotype_mate_info::polygenotype_mate_info
  (const polygenotype_mate_info & other)
{
  mate = other.mate;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
    {
      with_mate  [i][j] = other.with_mate  [i][j];  
      except_mate[i][j] = other.except_mate[i][j];  
    }
}

//===================================================================
//  polygenotype_mate_info::operator=(...)
//===================================================================
inline
individual_cache<genetic_info,log_double>::polygenotype_mate_info &
individual_cache<genetic_info,log_double>::polygenotype_mate_info::operator=
  (const polygenotype_mate_info & other)
{
  mate = other.mate;

  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
    {
      with_mate  [i][j] = other.with_mate  [i][j];  
      except_mate[i][j] = other.except_mate[i][j];  
    }

  return *this;
}

//===================================================================
//  individual_cache()  [ Constructor ]
//===================================================================
inline
individual_cache<genetic_info,log_double>::individual_cache()
{
  log_double QNAN = log_double(numeric_limits<double>::quiet_NaN());
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
      my_anterior[i][j] = my_posterior[i][j] = QNAN;
}

//===================================================================
//  individual_cache(...) [ COPY constructor ]
//===================================================================
inline
individual_cache<genetic_info,log_double>::individual_cache(const individual_cache& x)
{
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < MAX_POLYGENOTYPE; ++j)
    {
      my_anterior [i][j] = x.my_anterior [i][j];
      my_posterior[i][j] = x.my_posterior[i][j];
    }

  posterior_on_mate = x.posterior_on_mate;
}

//===================================================================
//  anterior_cached(...)
//===================================================================
inline
bool
individual_cache<genetic_info,log_double>::anterior_cached(
  const data_type & gen_info) const
{
  if(SAGE::isnan(my_anterior[gen_info.genotype][gen_info.polygenotype].get_double()))
    return false;
  else
    return true;
}

//===================================================================
//  posterior_cached(...)
//===================================================================
inline
bool
individual_cache<genetic_info,log_double>::posterior_cached(
  const data_type & gen_info) const
{
  if(SAGE::isnan(my_posterior[gen_info.genotype][gen_info.polygenotype].get_double())) 
    return false;
  else 
    return true;
}

//===================================================================
//  posterior_with_mate_cached(...) 
//===================================================================
inline
bool
individual_cache<genetic_info,log_double>::posterior_with_mate_cached(
  const member_type& mate, const data_type & gen_info) const
{
  for(gen_const_iterator 
      i  = posterior_on_mate.begin ();
      i != posterior_on_mate.end   (); ++i)

    if((i->mate == &mate) && 
       (!SAGE::isnan(i->with_mate[gen_info.genotype][gen_info.polygenotype].get_double())))
      return true;

  return false;
}

//===================================================================
//  posterior_except_mate_cached(...) 
//===================================================================
inline
bool
individual_cache<genetic_info,log_double>::posterior_except_mate_cached(
  const member_type& mate, const data_type & gen_info) const
{
  for(gen_const_iterator 
      i  = posterior_on_mate.begin ();
      i != posterior_on_mate.end   (); ++i)

    if((i->mate == &mate) && 
       (!SAGE::isnan(i->except_mate[gen_info.genotype][gen_info.polygenotype].get_double())))
      return true;

  return false;
}

//===================================================================
//  SET anterior(...) 
//===================================================================
inline
individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::anterior(
  const data_type & p)
{
  return my_anterior[p.genotype][p.polygenotype];
}

//===================================================================
//  SET posterior(...) 
//===================================================================
inline
individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior(const data_type & p)
{
  return my_posterior[p.genotype][p.polygenotype];
}

//===================================================================
//  SET posterior_with_mate(...) 
//===================================================================
inline
individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior_with_mate(
  const member_type& mate, const data_type & p)
{ 
  for(gen_iterator 
      i  = posterior_on_mate.begin (); 
      i != posterior_on_mate.end   (); ++i)

    if(i->mate == &mate) 
      return i->with_mate[p.genotype][p.polygenotype];

  posterior_on_mate.push_back(polygenotype_mate_info(&mate));
  return posterior_with_mate(mate,p);
}

//===================================================================
//  SET posterior_except_mate(...) 
//===================================================================
inline
individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior_except_mate(
  const member_type& mate, const data_type & p)
{ 
  for(gen_iterator 
      i  = posterior_on_mate.begin (); 
      i != posterior_on_mate.end   (); ++i)

    if(i->mate == &mate) 
      return i->except_mate[p.genotype][p.polygenotype];

  posterior_on_mate.push_back(polygenotype_mate_info(&mate));
  return posterior_except_mate(mate,p);
}

//===================================================================
//  GET anterior(...) 
//===================================================================
inline
const individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::anterior(
  const data_type & p) const
{
  return my_anterior[p.genotype][p.polygenotype]; 
}

//===================================================================
//  GET posterior(...) 
//===================================================================
inline
const individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior(
  const data_type & p) const
{
  return my_posterior[p.genotype][p.polygenotype]; 
}

//===================================================================
//  GET posterior_with_mate(...) 
//===================================================================
inline
const individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior_with_mate(
  const member_type& mate, const data_type & p) const

{
  genotype_index Ui = p.genotype;
  size_t         Vi = p.polygenotype;

  for(gen_const_iterator 
      i  = posterior_on_mate.begin (); 
      i != posterior_on_mate.end   (); ++i)

    if(i->mate == &mate) return i->with_mate[Ui][Vi];

  return posterior_on_mate.begin()->with_mate[Ui][Vi];
}

//===================================================================
//  GET posterior_except_mate(...) 
//===================================================================
inline
const individual_cache<genetic_info,log_double>::result_type &
individual_cache<genetic_info,log_double>::posterior_except_mate(
  const member_type& mate, const data_type & p) const
{
  genotype_index Ui = p.genotype;
  size_t         Vi = p.polygenotype;

  for(gen_const_iterator 
      i  = posterior_on_mate.begin (); 
      i != posterior_on_mate.end   (); ++i)

    if(i->mate == &mate) return i->except_mate[Ui][Vi];

  return posterior_on_mate.begin()->except_mate[Ui][Vi];
}



}}
