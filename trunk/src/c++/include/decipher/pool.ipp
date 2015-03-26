//============================================================================
// File:      pool.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/9/6 - created.                                   djb
//                                                                          
// Notes:      
//                                                                          
// Copyright (c) 2006 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  pool_pheno_seq_generator
//============================================================================
//
inline
pool_pheno_seq_generator::pool_pheno_seq_generator(const vector<member>& members,
                                         const locus_group& loci, cerrorstream& errors,
                                         ostream& messages,
                                         const instructions& instr,
                                         output_state& ostate)
      : my_members(members), my_loci(loci), my_errors(errors), my_messages(messages),
        my_instructions(instr), my_output_state(ostate)
{}

inline cerrorstream&
pool_pheno_seq_generator::errors()
{
  return  my_errors;
}

inline ostream&
pool_pheno_seq_generator::messages()
{
  return  my_messages;
}

inline const instructions&
pool_pheno_seq_generator::instr() const
{
  return  my_instructions;
}

inline bool&
pool_pheno_seq_generator::messages_interrupted()
{
  return  my_output_state.msg_interrupted;
}

inline const vector<member>&
pool_pheno_seq_generator::members() const
{
  return  my_members;
}

inline const locus_group&
pool_pheno_seq_generator::loci() const
{
  return  my_loci;
}

inline pool_pheno_seq_generator::iterator
pool_pheno_seq_generator::begin()
{
  return  iterator(this);
}

inline pool_pheno_seq_generator::iterator
pool_pheno_seq_generator::end()
{
  return  iterator(this, true);
}


//============================================================================
// IMPLEMENTATION:  pool_pheno_seq_generator::iterator
//============================================================================
//
inline
pool_pheno_seq_generator::iterator::iterator(pool_pheno_seq_generator* generator, bool end)
     : my_generator(generator)
{
  my_member_end_iter = my_generator->members().end();

  if(end)
  {
    my_member_iter = my_generator->members().end();
  }
  else
  {
    my_member_iter = my_generator->members().begin();
    if(my_member_iter != my_member_end_iter)
    {
      const locus_group&  loci = my_generator->loci();
    
      locus_group::const_iterator  l_iter     = loci.begin();
      locus_group::const_iterator  l_end_iter = loci.end();
      for(; l_iter != l_end_iter; ++l_iter)
      {
        my_pool_pheno_seq.push_back(make_pair((*l_iter).second, map<size_t, size_t>()));
      }        
    
      if(! build_seq())
      {
        ++(*this);
      }   
    }
  }
}

inline pool_pheno_seq_generator::iterator
pool_pheno_seq_generator::iterator::operator ++()
{
  if(my_member_iter != my_member_end_iter)
  {
    ++my_member_iter;
  }
  
  if(my_member_iter != my_member_end_iter)
  {
    if(! build_seq())
    {
      ++(*this);
    }
  }
  
  return  *this;
}

inline member
pool_pheno_seq_generator::iterator::get_member() const
{
  return  *my_member_iter;
}

inline size_t
pool_pheno_seq_generator::iterator::phenotype_count() const
{
  return  my_expanded_pheno_seq_count;
}

inline size_t  
pool_pheno_seq_generator::iterator::pool_size() const
{
  return  my_generator->instr().pool_size;
}

// - How many possible genotypes are there at a locus with n possible
//   alleles and k chromosomes?  This equation courtesy of Rob Igo.
//
inline size_t
pool_pheno_seq_generator::iterator::genotype_count(size_t n, size_t k)
{
  double  count = choose(n + k -1, k);
  
  if(! SAGE::isnan(count))
  {
    return  static_cast<size_t>(count);
  }
  else
  {
    return  (size_t)(-1);
  }
}

inline bool
pool_pheno_seq_generator::iterator::operator ==(const iterator& other) const
{
  return  my_generator == other.my_generator &&
          my_member_iter == other.my_member_iter;
}

inline bool
pool_pheno_seq_generator::iterator::operator !=(const iterator& other) const
{
  return  ! operator ==(other);
}

inline pool_pheno_seq&
pool_pheno_seq_generator::iterator::operator *()
{
  return  my_pool_pheno_seq;
}

inline pool_pheno_seq*
pool_pheno_seq_generator::iterator::operator ->()
{
  return  &my_pool_pheno_seq;
}


//============================================================================
// IMPLEMENTATION:  pool_em_phenotype_map
//============================================================================
//
inline
pool_em_phenotype_map::pool_em_phenotype_map(output_state& ostate,
                                   APP::Output_Streams& streams,
                                   const FilteredMultipedigree& mped,
                                   const vector<member>& members, 
                                   const instructions& instr, const locus_group& loci,
                                   const string& inner_sub_pop_name, const string& outer_sub_pop_name)
      : member_em_phenotype_map(ostate, streams, members, loci, inner_sub_pop_name, outer_sub_pop_name), my_mped(mped), 
        my_instructions(instr)  
{
  build();
}

inline size_t  
pool_em_phenotype_map::pheno_index(member m) const
{
  std::map<member, pool_pheno_seq>::const_iterator  m_iter     = my_member_directory.find(m);
  std::map<member, pool_pheno_seq>::const_iterator  m_end_iter = my_member_directory.end();
  assert(m_iter != m_end_iter);
  
  std::map<pool_pheno_seq, size_t>::const_iterator  p_iter     = my_phenotype_directory.find(m_iter->second);
  std::map<pool_pheno_seq, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
  assert(p_iter != p_end_iter);
  
  size_t  phenotype_index = p_iter->second;
  assert(phenotype_index < my_phenotypes.size());
  
  return  phenotype_index;
}

inline const em_phenotype&  
pool_em_phenotype_map::operator [](member m) const
{
  return  my_phenotypes[pheno_index(m)];
}

inline const em_phenotype&  
pool_em_phenotype_map::final_phenotype(member m) const
{
  assert(maximized);
  return  my_best_phenotypes[pheno_index(m)];
}

inline const FilteredMultipedigree&
pool_em_phenotype_map::multipedigree() const
{
  return  my_mped;
}



