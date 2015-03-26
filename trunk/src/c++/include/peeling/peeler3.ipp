//========================
// Function Instantiations
//========================

template <class Data, class Result, class IndCache>
peeler<Data, Result, IndCache>::peeler(const subped_type& subped)
  : my_subpedigree(subped), my_cache(subped)
{ }

template <class Data, class Result, class IndCache>
peeler<Data, Result, IndCache>::~peeler() { }

template <class Data, class Result, class IndCache>
const typename peeler<Data, Result, IndCache>::result_type& peeler<Data, Result, IndCache>::anterior
    (const member_type& ind, const data_type& g)
{
  // Get the individual cache
  
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  // Test if it's already available
  
  if(c.anterior_cached(g)) return c.anterior(g);
  
  // If not, is the person a founder?
  
  result_type& r = c.anterior(g);

  if(ind.parent_begin() == ind.parent_end()) internal_anterior_terminal(ind, g, r);
  else                                       internal_anterior         (ind, g, r);
  
  // And return the results.

  return r;
}

template <class Data, class Result, class IndCache>
const typename peeler<Data, Result, IndCache>::result_type& peeler<Data, Result, IndCache>::posterior
   (const member_type& ind, const data_type& g)
{
  // Get the individual cache
  
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  // Test if it's already available
  
  if(c.posterior_cached(g)) return c.posterior(g);
  
  // If not, is the person a leaf?
  
  result_type& r = c.posterior(g);

  if(ind.mate_count() == 0) internal_posterior_terminal(ind, g, r);
  else                      internal_posterior         (ind, g, r);
  
  // And return the results.

  return r;
}


template <class Data, class Result, class IndCache>
const typename peeler<Data, Result, IndCache>::result_type& peeler<Data, Result, IndCache>::posterior_with_mate
   (const member_type& ind, const member_type& mate, const data_type& g)
{
  // Get the individual cache
  
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  // Test if it's already available

  if(c.posterior_with_mate_cached(mate, g))
    return c.posterior_with_mate(mate, g);
  
  result_type& r = c.posterior_with_mate(mate, g);

  // Calculate
  
  internal_posterior_with_mate(ind, mate, g, r);
  
  // And return the results.

  return r;
}

template <class Data, class Result, class IndCache>
const typename peeler<Data, Result, IndCache>::result_type&
  peeler<Data, Result, IndCache>::posterior_except_mate
    (const member_type& ind, const member_type& mate, const data_type& g)
{
  // Get the individual cache
  
  individual_cache_type& c = my_cache.get_individual_cache(ind);

  // Test if it's already available
  
  if(c.posterior_except_mate_cached(mate, g)) 
    return c.posterior_except_mate(mate, g);
  
  result_type& r = c.posterior_except_mate(mate, g);

  // Calculate
  
  internal_posterior_except_mate(ind, mate, g, r);
  
  // And return the results.

  return r;
}


