//============================================================================
// File:      em.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/26/4 - created.                                   djb
//                                                                          
// Notes:      
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION: member_order
//============================================================================
//
template <typename T> inline bool
member_order<T>::operator ()(T m1, T m2) const
{
  const string&  pedigree1 = m1->pedigree()->name();
  const string&  pedigree2 = m2->pedigree()->name();

  if(pedigree1 == pedigree2)
  {
    return  m1->name() < m2->name();
  }
  else
  {
    return  pedigree1 < pedigree2;
  }
}

//============================================================================
// IMPLEMENTATION:  em_haplotype
//============================================================================
//
inline
em_haplotype::em_haplotype()
      : my_sequence(0), my_old_count(0), my_new_count(0), my_static_count(0) 
{}

inline
em_haplotype::em_haplotype(const hap_seq* seq)
      : my_sequence(seq), my_old_count(0), my_new_count(0), my_static_count(0) 
{}

inline double  
em_haplotype::old_freq(size_t total_hap_count) const
{
  assert(total_hap_count != 0);
  return  my_old_count / total_hap_count;
}
    
inline double  
em_haplotype::new_freq(size_t total_hap_count) const
{
  assert(total_hap_count != 0);
  return  my_new_count / total_hap_count;
}
    
inline const hap_seq&
em_haplotype::sequence() const
{
  return  *my_sequence;
}
    
inline void  
em_haplotype::incr(double amount)
{
  my_new_count += amount;
}

inline void  
em_haplotype::reset()
{
  my_old_count = my_new_count;
  my_new_count = my_static_count;
}
    
inline void  
em_haplotype::zero()
{
  my_old_count = 0.0;
  my_new_count = 0.0;
}

// - Intra-haplotype frequency comparison.
//
inline bool  
em_haplotype::converged(size_t total_hap_count, double epsilon) const
{
  return  abs(new_freq(total_hap_count) - old_freq(total_hap_count)) < epsilon;
}

// - Inter-haplotype frequency comparison.
//
inline bool  
em_haplotype::similar(const em_haplotype& other, size_t total_hap_count, double epsilon) const
{
  return  abs(new_freq(total_hap_count) - other.new_freq(total_hap_count)) < epsilon;
}

//============================================================================
// IMPLEMENTATION:  hap_freq
//============================================================================
//
inline
hap_freq::hap_freq(const string& haplotype, double frequency)
      : hap(haplotype), freq(frequency)
{}

inline bool
hap_freq::operator >(const hap_freq& other) const
{
  if(freq != other.freq)
  {
    return  freq > other.freq;
  }
  else   
  {
    return  hap > other.hap;
  }
}

//============================================================================
// IMPLEMENTATION:  em_haplotype_map
//============================================================================
//
inline
em_haplotype_map::em_haplotype_map()
      : counts_initialized(false)
{}
    
inline bool  
em_haplotype_map::contains(const hap_seq& key) const
{
  return  my_hap_seqs.find(key) != my_hap_seqs.end();
}

// - Intra-haplotype frequency comparison.
//
inline bool  
em_haplotype_map::converged(size_t total_hap_count, double epsilon) const
{
  bool  conv = true;
  
  std::map<hap_seq, size_t>::const_iterator  iter = my_hap_seqs.begin();
  std::map<hap_seq, size_t>::const_iterator  end_iter = my_hap_seqs.end();
  for(; iter != end_iter; iter++)
  {
    size_t  index = iter->second;
    assert(my_haplotypes.find(index) != my_haplotypes.end());
    if(! my_haplotypes[index].converged(total_hap_count, epsilon))
    {
      conv = false;
      break;
    }
  }
  
  return  conv;
}

// - Inter-haplotype frequency comparison.
//
inline bool  
em_haplotype_map::similar(em_haplotype_map& other, size_t total_hap_count, double epsilon) const
{
  bool  sim = true;
  
  std::map<size_t, em_haplotype>::const_iterator  iter = my_haplotypes.begin();
  std::map<size_t, em_haplotype>::const_iterator  end_iter = my_haplotypes.end();
  for(; iter != end_iter; iter++)
  {
    const hap_seq&  seq = iter->second.sequence();
    assert(other.contains(seq));
    if(! iter->second.similar(other[other.hap_seq_to_index(seq)], total_hap_count, epsilon))
    {
      sim = false;
      break;
    }
  }
  
  return  sim;
}

inline double  
em_haplotype_map::total_freq(size_t total_hap_count) const
{
  double  total = 0.0;
  
  std::map<size_t, em_haplotype>::const_iterator  iter = my_haplotypes.begin();
  std::map<size_t, em_haplotype>::const_iterator  end_iter = my_haplotypes.end();
  for(; iter != end_iter; iter++)
  {
    total += iter->second.new_freq(total_hap_count);
  }
  
  return  total;
}

// - If haplotype doesn't exist, create it.  Return haplotype index.
//
inline size_t 
em_haplotype_map::add_haplotype(const hap_seq& key)
{
  typedef std::map<hap_seq, size_t>::iterator  seq_iterator;
  typedef std::map<size_t, em_haplotype>::iterator  hap_iterator;
  
  size_t  index;
  
  seq_iterator  s_iter = my_hap_seqs.find(key);
  seq_iterator  s_end_iter = my_hap_seqs.end();
  if(s_iter == s_end_iter)
  {
    index = my_hap_seqs.size();
    pair<seq_iterator, bool>  seq_insert_result = my_hap_seqs.insert(make_pair(key, index));
    assert(seq_insert_result.second);    
  
    pair<hap_iterator, bool>  hap_insert_result = 
              my_haplotypes.insert(make_pair(index, em_haplotype(&(seq_insert_result.first->first))));
    assert(hap_insert_result.second);
  }
  else
  {
    index = s_iter->second;
  }
  
  return  index;
}

inline em_haplotype&
em_haplotype_map::operator [](size_t index) 
{
  assert(my_haplotypes.find(index) != my_haplotypes.end());
  
  return  my_haplotypes[index];
}

inline const em_haplotype&
em_haplotype_map::get_haplotype(size_t index) const 
{
  assert(my_haplotypes.find(index) != my_haplotypes.end());
  
  return  my_haplotypes[index];
}

inline const hap_seq&
em_haplotype_map::index_to_hap_seq(size_t index) const
{
  return  my_haplotypes[index].sequence();
}

inline size_t
em_haplotype_map::hap_seq_to_index(const hap_seq& key) const
{
  if(contains(key))
  {
    return  my_hap_seqs[key];
  }
  else
  {
    return  (size_t)(-1);
  }
}

inline const std::map<size_t, em_haplotype>&
em_haplotype_map::haplotypes() const
{
  return  my_haplotypes;
}
    
inline void  
em_haplotype_map::reset()
{
  std::map<size_t, em_haplotype>::iterator  iter = my_haplotypes.begin();
  std::map<size_t, em_haplotype>::iterator  end_iter = my_haplotypes.end();
  for(; iter != end_iter; iter++)
  {
    iter->second.reset();
  }
}

inline void  
em_haplotype_map::zero()
{
  std::map<size_t, em_haplotype>::iterator  iter = my_haplotypes.begin();
  std::map<size_t, em_haplotype>::iterator  end_iter = my_haplotypes.end();
  for(; iter != end_iter; iter++)
  {
    iter->second.zero();
  }
}


//============================================================================
// IMPLEMENTATION:  comb_prob
//============================================================================
//
inline
comb_prob::comb_prob(string c, double p)
      : comb(c), prob(p)
{}

inline bool
comb_prob::operator >(const comb_prob& other) const
{
  if(prob != other.prob)
  {
    return  prob > other.prob;
  }
  else
  {
    // - Sort in ascending order if probabilities are the same!
    //
    return  comb < other.comb;
  }
}

//============================================================================
// IMPLEMENTATION:  em_phenotype
//============================================================================
//
inline
em_phenotype::em_phenotype()
{}

inline size_t
em_phenotype::count() const
{
  return  my_count;
}

inline bool
em_phenotype::is_ambiguous() const
{
  return  ambiguous;
}

inline const vector<pair<em_phenotype::combination, double> >&
em_phenotype::combinations() const
{
  return  my_combinations;
}

inline size_t
em_phenotype::comb_size() const
{
  assert(! my_combinations.empty());
  
  return  my_combinations[0].first.size();
}

inline void
em_phenotype::incr()
{
  my_count++;
}


//============================================================================
// IMPLEMENTATION:  em_phenotype::combination_iterator
//============================================================================
//
inline pair<size_t, size_t>
em_phenotype::combination_iterator::operator *()
{
  return  *my_internal_iterator;
}

inline std::map<size_t, size_t>::iterator
em_phenotype::combination_iterator::operator ++(int)
{
  return  my_internal_iterator++;
}

inline bool
em_phenotype::combination_iterator::at_end() const
{
  return  my_internal_iterator == my_distinct_haplotypes.end();
}


//============================================================================
// IMPLEMENTATION:  base_em_phenotype_map
//============================================================================
//
inline
base_em_phenotype_map::base_em_phenotype_map(output_state& ostate,
                                             APP::Output_Streams& streams, const locus_group& loci,
                                             const string& inner_sub_pop_name, 
                                             const string& outer_sub_pop_name)
      : my_output_state(ostate), my_streams(streams), my_errors(streams.errors()), my_messages(streams.messages()), 
        my_loci(loci), my_inner_sub_pop_name(inner_sub_pop_name), my_outer_sub_pop_name(outer_sub_pop_name), 
        my_total_hap_count(0), my_max_ln_likelihood(QNAN), maximized(false)
{
  build_sub_pop_name();
}

inline
base_em_phenotype_map::~base_em_phenotype_map()
{}

inline void
base_em_phenotype_map::build_sub_pop_name()
{
  if(! (my_outer_sub_pop_name.empty() || my_inner_sub_pop_name.empty()))
  {
    my_sub_pop_name = my_outer_sub_pop_name + " / " + my_inner_sub_pop_name;
  }
  else if(my_outer_sub_pop_name.empty())
  {
    my_sub_pop_name = my_inner_sub_pop_name;
  }
  else
  {
    my_sub_pop_name = "";
  }
}

inline base_em_phenotype_map::const_iterator
base_em_phenotype_map::begin() const
{
  return  my_phenotypes.begin();
}

inline base_em_phenotype_map::const_iterator
base_em_phenotype_map::end() const
{
  return  my_phenotypes.end();
}

inline base_em_phenotype_map::const_iterator
base_em_phenotype_map::final_begin() const
{
  return  my_best_phenotypes.begin();
}

inline base_em_phenotype_map::const_iterator
base_em_phenotype_map::final_end() const
{
  return  my_best_phenotypes.end();
}

inline bool
base_em_phenotype_map::empty() const
{
  return  my_phenotypes.empty();
}

inline const string&
base_em_phenotype_map::inner_sub_pop_name() const
{
  return  my_inner_sub_pop_name;
}

inline const string&
base_em_phenotype_map::outer_sub_pop_name() const
{
  return  my_outer_sub_pop_name;
}

inline const string&
base_em_phenotype_map::sub_pop_name() const
{
  return  my_sub_pop_name;
}

inline const locus_group&
base_em_phenotype_map::loci() const
{
  return  my_loci;
}

inline size_t
base_em_phenotype_map::total_hap_count() const
{
  return  my_total_hap_count;
}

inline const em_haplotype_map&
base_em_phenotype_map::haplotypes() const
{
  return  my_haplotypes;
}

inline const em_haplotype_map&
base_em_phenotype_map::final_haplotypes() const
{
  assert(maximized);
  return  my_best_haplotypes;
}

inline double
base_em_phenotype_map::max_ln_likelihood() const
{
  assert(maximized);
  return  my_max_ln_likelihood;
}

// - Return number of haplotype frequencies - 1.
//
inline size_t
base_em_phenotype_map::independent_param_count() const
{
  return  my_haplotypes.haplotypes().size() - 1;
}

inline void
base_em_phenotype_map::init_weights()
{
  vector<em_phenotype>::iterator  p_iter     = my_phenotypes.begin();
  vector<em_phenotype>::iterator  p_end_iter = my_phenotypes.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    p_iter->init_weights();
  }
}

inline void
base_em_phenotype_map::update_weights()
{
  vector<em_phenotype>::iterator  p_iter     = my_phenotypes.begin();
  vector<em_phenotype>::iterator  p_end_iter = my_phenotypes.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    p_iter->update_weights(my_total_hap_count, my_haplotypes);
  }
}

inline bool
base_em_phenotype_map::is_maximized() const
{
  return  maximized;
}


//============================================================================
// IMPLEMENTATION:  member_em_phenotype_map
//============================================================================
//
inline
member_em_phenotype_map::member_em_phenotype_map(output_state& ostate,
                                                 APP::Output_Streams& streams,
                                                 const vector<member>& members,
                                                 const locus_group& loci, 
                                                 const string& inner_sub_pop_name,
                                                 const string& outer_sub_pop_name)
      : base_em_phenotype_map(ostate, streams, loci, inner_sub_pop_name, outer_sub_pop_name), 
        my_member_pool(members)
{}      

inline
member_em_phenotype_map::~member_em_phenotype_map()
{}

//============================================================================
// IMPLEMENTATION:  pop_freq_writer
//============================================================================
//
inline
pop_freq_writer::pop_freq_writer(ostream& out, double cutoff)
      : my_out(out), my_frequencies(0), my_cutoff(cutoff), 
        my_ln_likelihood(QNAN)
{}

inline void
pop_freq_writer::cutoff_note(ostream& out)
{
    out << "Note:  Haplotypes listed have estimated frequencies greater than or equal to\n"
        << "       the cutoff or have the greatest frequency estimate.\n" << endl;
}

inline void
pop_freq_writer::set_frequencies(const set<hap_freq, greater<hap_freq> >* frequencies)
{
  my_frequencies = frequencies;
}                                 

inline void
pop_freq_writer::set_ln_likelihood(double ln_likelihood)
{
  my_ln_likelihood = ln_likelihood;
}
