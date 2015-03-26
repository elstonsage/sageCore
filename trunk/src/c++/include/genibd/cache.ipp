//============================================================================
// IMPLEMENTATION:  pv_equal
//============================================================================
//
template<class T1, class T2> inline bool 
pv_equal<T1, T2>::operator()(T1 first, T2 second) const
{
  return first.first == second;
}

//============================================================================
// IMPLEMENTATION:  individual_cache
//============================================================================
//
inline
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::individual_cache()
{}

// - Elements of m_indices correspond to possible marker genotypes 
//   and index into anteriors, posteriors, etc.
//   Client code should only use penetrance iterators whose
//   genotypes are penetrant to access cached values.
//
inline void
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::build
(const member_type& ind, size_t mph_id, const MLOCUS::inheritance_model& mm)
{
  //cout << "build cache for " << ind.name() << endl;

  size_t  pmg_count = 0;           // Penetrant marker genotype count.

  data_type  mp_iter = mm.phased_penetrance_begin(mph_id);
  
  for( ; mp_iter != mm.phased_penetrance_end(mph_id); ++mp_iter )
  {
    //cout << "("
    //     << mp_iter.geno_id() << ":" << mp_iter.phased_geno().name() << ", "
    //     << mp_iter.phenotype_id() << ":" << mp_iter.pheno().name() << ") ";

    m_indices.insert(make_pair(mp_iter.geno_id(), pmg_count++));
  }

  //cout << endl << "pmg_count = " << pmg_count << endl;

  my_anteriors.resize(pmg_count, log_double(QNAN));
  my_posteriors.resize(pmg_count, log_double(QNAN));
  
  size_t  p_size = ind.mate_count();

  my_posteriors_with_mate.resize(  p_size, make_pair((size_t)(-1), 
                                   vector<log_double>(pmg_count, log_double(QNAN))));
  my_posteriors_except_mate.resize(p_size, make_pair((size_t)(-1), 
                                   vector<log_double>(pmg_count, log_double(QNAN))));
  
  FPED::MateConstIterator  iter;
  size_t i = 0;
  for( iter = ind.mate_begin(); iter != ind.mate_end(); ++i, ++iter )
  {
    assert(i < my_posteriors_with_mate.size());
    my_posteriors_with_mate[i].first = iter->mate().subindex();
    
    assert(i < my_posteriors_except_mate.size()); 
    my_posteriors_except_mate[i].first = iter->mate().subindex(); 
  }
}

//
// ------------------------- cached? ------------------------------------
//
inline bool
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::anterior_cached(const data_type& upi) const
{
  size_t  index = genotype_index(upi);
  assert(index < my_anteriors.size());
  bool  cached = ! SAGE::isnan(my_anteriors[index].get_double());
  
  //cout << (cached ? " ANTERIOR CACHED " : "ANTERIOR NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_cached(const data_type& upi) const
{
  size_t  index = genotype_index(upi);
  assert(index < my_posteriors.size());
  bool  cached = ! SAGE::isnan(my_posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR CACHED " : " POSTERIOR NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_with_mate_cached
(const member_type& mate_index, const data_type& upi) const
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
  bool  cached = ! SAGE::isnan(posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR WITH MATE CACHED " : " POSTERIOR WITH MATE NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_except_mate_cached
(const member_type& mate_index, const data_type& upi) const
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
  bool  cached = ! SAGE::isnan(posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR EXCEPT MATE CACHED " : " POSTERIOR EXCEPT MATE NOT CACHED ") << endl;
  
  return cached;
}

//
// ----------------------- get non-const reference ------------------------
//
inline log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::anterior(const data_type& upi) 
{
  size_t  index = genotype_index(upi);
  assert(index < my_anteriors.size());
  
  return my_anteriors[index];
}

inline log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior(const data_type& upi) 
{
  size_t  index = genotype_index(upi);
  assert(index < my_posteriors.size());
  
  return my_posteriors[index];
}

inline log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_with_mate
(const member_type& mate_index, const data_type& upi)
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

inline log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_except_mate
(const member_type& mate_index, const data_type& upi)
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

//
// ----------------------- get const reference ------------------------
//
inline const log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::anterior(const data_type& upi) const 
{
  size_t  index = genotype_index(upi);
  assert(index < my_anteriors.size());
  
  return my_anteriors[index];
}

inline const log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior(const data_type& upi) const
{
  size_t  index = genotype_index(upi);
  assert(index < my_posteriors.size());
  
  return my_posteriors[index];
}

inline const log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_with_mate
(const member_type& mate_index, const data_type& upi) const
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

inline const log_double&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posterior_except_mate
(const member_type& mate_index, const data_type& upi) const
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.subindex());
  
  size_t  index = genotype_index(upi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

//
// ----------------------- ancillary functions ----------------------------
//
inline size_t
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::genotype_index(const data_type& upi) const
{
  map<size_t, size_t>::const_iterator  col = m_indices.find(upi.geno_id());
  assert(col != m_indices.end());
  
  return col->second;
}

// - Find posteriors corresponding to a specific mate.
//
inline vector<log_double>&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posteriors_with_mate(size_t mate) const
{
  vector<posterior_vector>::iterator  iter;
  iter = find_if(my_posteriors_with_mate.begin(), my_posteriors_with_mate.end(),
                 bind2nd(pv_equal<posterior_vector, size_t>(), mate));
  assert(iter != my_posteriors_with_mate.end());
  
  return iter->second;
}

// - Find posteriors corresponding to a specific mate.
//
inline vector<log_double>&
individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator, log_double>::posteriors_except_mate(size_t mate) const
{
  vector<posterior_vector>::iterator  iter;
  iter = find_if(my_posteriors_except_mate.begin(), my_posteriors_except_mate.end(),
                 bind2nd(pv_equal<posterior_vector, size_t>(), mate));
  assert(iter != my_posteriors_except_mate.end());
  
  return iter->second;
}


