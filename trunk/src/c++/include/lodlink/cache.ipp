//============================================================================
// File:      cache.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/17/2 - created.                                djb
//                                                                          
// Notes:     inlines for individual_cache.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



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
individual_cache<joint_pen_iter, log_double>::individual_cache()
{}

// - Elements of t_indices and m_indices correspond to possible trait and 
//   marker genotypes respectively.  Their sum is an index into anteriors, 
//   posteriors, etc.  Client code should only use penetrance iterators whose
//   genotypes are penetrant to access cached values.
//
inline void
individual_cache<joint_pen_iter, log_double>::build
(const member_type& ind, size_t tph_id, size_t mph_id, 
 const MLOCUS::penetrance_model& tm, const MLOCUS::penetrance_model& mm)
{
  size_t  pmg_count = 0;           // Penetrant marker genotype count.
  MLOCUS::penetrance_model::phased_penetrance_iterator  mp_iter;
  
  for(mp_iter = mm.phased_penetrance_begin(mph_id); 
      mp_iter != mm.phased_penetrance_end(mph_id); ++mp_iter)
  {
    m_indices.insert(make_pair(mp_iter.geno_id(), pmg_count++));
  }
  
  size_t  ptg_count = 0;           // Penetrant trait genotype count.
  MLOCUS::penetrance_model::phased_penetrance_iterator  tp_iter;
  
  for(tp_iter = tm.phased_penetrance_begin(tph_id);
      tp_iter != tm.phased_penetrance_end(tph_id); ++ tp_iter)
  {
    t_indices.insert(make_pair(tp_iter.geno_id(), ptg_count++ * pmg_count));
  }
  
  size_t  phenoset_size = pmg_count * ptg_count;
  
  my_anteriors.resize(phenoset_size, log_double(QNAN));
  my_posteriors.resize(phenoset_size, log_double(QNAN));
  
  size_t  p_size = ind.mate_count();
  my_posteriors_with_mate.resize(p_size, make_pair((size_t)(-1), 
                                   vector<log_double>(phenoset_size, log_double(QNAN))));
  my_posteriors_except_mate.resize(p_size, make_pair((size_t)(-1), 
                                   vector<log_double>(phenoset_size, log_double(QNAN))));
  
  FPED::Pedigree::mate_const_iterator  iter;
  size_t  i = 0;
  for(iter = ind.mate_begin(); iter != ind.mate_end(); ++i, ++iter)
  {
    assert(i < my_posteriors_with_mate.size());
    my_posteriors_with_mate[i].first = iter->mate().index();
    
    assert(i < my_posteriors_except_mate.size()); 
    my_posteriors_except_mate[i].first = iter->mate().index(); 
  }
}

//
// ------------------------- cached? ------------------------------------
//
inline bool
individual_cache<joint_pen_iter, log_double>::anterior_cached(const joint_pen_iter& jpi) const
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_anteriors.size());
  bool  cached = ! SAGE::isnan(my_anteriors[index].get_double());
  
  //cout << (cached ? " ANTERIOR CACHED " : "ANTERIOR NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<joint_pen_iter, log_double>::posterior_cached(const joint_pen_iter& jpi) const
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_posteriors.size());
  bool  cached = ! SAGE::isnan(my_posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR CACHED " : " POSTERIOR NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<joint_pen_iter, log_double>::posterior_with_mate_cached
(const member_type& mate_index, const joint_pen_iter& jpi) const
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
  bool  cached = ! SAGE::isnan(posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR WITH MATE CACHED " : " POSTERIOR WITH MATE NOT CACHED ") << endl;
  
  return cached;
}

inline bool
individual_cache<joint_pen_iter, log_double>::posterior_except_mate_cached
(const member_type& mate_index, const joint_pen_iter& jpi) const
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
  bool  cached = ! SAGE::isnan(posteriors[index].get_double());
  
  //cout << (cached ? " POSTERIOR EXCEPT MATE CACHED " : " POSTERIOR EXCEPT MATE NOT CACHED ") << endl;
  
  return cached;
}

//
// ----------------------- get non-const reference ------------------------
//
inline log_double&
individual_cache<joint_pen_iter, log_double>::anterior(const joint_pen_iter& jpi) 
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_anteriors.size());
  
  return my_anteriors[index];
}

inline log_double&
individual_cache<joint_pen_iter, log_double>::posterior(const joint_pen_iter& jpi) 
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_posteriors.size());
  
  return my_posteriors[index];
}

inline log_double&
individual_cache<joint_pen_iter, log_double>::posterior_with_mate
(const member_type& mate_index, const joint_pen_iter& jpi)
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

inline log_double&
individual_cache<joint_pen_iter, log_double>::posterior_except_mate
(const member_type& mate_index, const joint_pen_iter& jpi)
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

//
// ----------------------- get const reference ------------------------
//
inline const log_double&
individual_cache<joint_pen_iter, log_double>::anterior(const joint_pen_iter& jpi) const 
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_anteriors.size());
  
  return my_anteriors[index];
}

inline const log_double&
individual_cache<joint_pen_iter, log_double>::posterior(const joint_pen_iter& jpi) const
{
  size_t  index = phenoset_index(jpi);
  assert(index < my_posteriors.size());
  
  return my_posteriors[index];
}

inline const log_double&
individual_cache<joint_pen_iter, log_double>::posterior_with_mate
(const member_type& mate_index, const joint_pen_iter& jpi) const
{
  vector<log_double>&  posteriors = posteriors_with_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

inline const log_double&
individual_cache<joint_pen_iter, log_double>::posterior_except_mate
(const member_type& mate_index, const joint_pen_iter& jpi) const
{
  vector<log_double>&  posteriors = posteriors_except_mate(mate_index.index());
  
  size_t  index = phenoset_index(jpi);
  assert(index < posteriors.size());
    
  return  posteriors[index];
}

//
// ----------------------- ancillary functions ----------------------------
//
inline size_t
individual_cache<joint_pen_iter, log_double>::phenoset_index(const joint_pen_iter& jpi) const
{
  std::map<size_t, size_t>::const_iterator  row = t_indices.find(jpi.trait_iter.geno_id());
  std::map<size_t, size_t>::const_iterator  col = m_indices.find(jpi.marker_iter.geno_id());
  assert(row != t_indices.end());
  assert(col != m_indices.end());
  
  return row->second + col->second;
}

// - Find posteriors corresponding to a specific mate.
//
inline vector<log_double>&
individual_cache<joint_pen_iter, log_double>::posteriors_with_mate(size_t mate) const
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
individual_cache<joint_pen_iter, log_double>::posteriors_except_mate(size_t mate) const
{
  vector<posterior_vector>::iterator  iter;
  iter = find_if(my_posteriors_except_mate.begin(), my_posteriors_except_mate.end(),
                 bind2nd(pv_equal<posterior_vector, size_t>(), mate));
  assert(iter != my_posteriors_except_mate.end());
  
  return iter->second;
}


