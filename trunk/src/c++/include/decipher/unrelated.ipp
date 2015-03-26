//============================================================================
// File:      unrelated.ipp
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
// IMPLEMENTATION:  pheno_seq_generator
//============================================================================
//
inline
pheno_seq_generator::pheno_seq_generator(const vector<member>& members,
                                         const locus_group& loci)
      : my_members(members), my_loci(loci)
{}

inline const vector<member>&
pheno_seq_generator::members() const
{
  return  my_members;
}

inline const locus_group&
pheno_seq_generator::loci() const
{
  return  my_loci;
}

inline pheno_seq_generator::iterator
pheno_seq_generator::begin() const
{
  return  iterator(this);
}

inline pheno_seq_generator::iterator
pheno_seq_generator::end() const
{
  return  iterator(this, true);
}


//============================================================================
// IMPLEMENTATION:  pheno_seq_generator::iterator
//============================================================================
//
inline
pheno_seq_generator::iterator::iterator(const pheno_seq_generator* generator, bool end)
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
        my_pheno_seq.push_back(make_pair((*l_iter).second, (size_t)(-1)));
      }    
    
      build_seq();
    }
  }
}

inline pheno_seq_generator::iterator
pheno_seq_generator::iterator::operator ++()
{
  ++my_member_iter;
  build_seq();
  
  return  *this;
}

inline member
pheno_seq_generator::iterator::get_member() const
{
  return  *my_member_iter;
}

inline bool
pheno_seq_generator::iterator::operator ==(const iterator& other) const
{
  return  my_generator == other.my_generator &&
          my_member_iter == other.my_member_iter;
}

inline bool
pheno_seq_generator::iterator::operator !=(const iterator& other) const
{
  return  ! operator ==(other);
}

inline pheno_seq&
pheno_seq_generator::iterator::operator *()
{
  return  my_pheno_seq;
}

inline pheno_seq*
pheno_seq_generator::iterator::operator ->()
{
  return  &my_pheno_seq;
}


//============================================================================
// IMPLEMENTATION:  unrelated_em_phenotype_map
//============================================================================
//
inline
unrelated_em_phenotype_map::unrelated_em_phenotype_map(output_state& ostate,
                                   APP::Output_Streams& streams,
                                   const FilteredMultipedigree& mped,
                                   const vector<member>& members, 
                                   const instructions& instr, const locus_group& loci,
                                   const string& inner_sub_pop_name, 
                                   const string& outer_sub_pop_name)
      : member_em_phenotype_map(ostate, streams, members, loci, inner_sub_pop_name, outer_sub_pop_name), 
        my_mped(mped), my_instructions(instr)  
{
  build();
}

inline size_t  
unrelated_em_phenotype_map::pheno_index(member m) const
{
  std::map<member, pheno_seq>::const_iterator  m_iter     = my_member_directory.find(m);
  std::map<member, pheno_seq>::const_iterator  m_end_iter = my_member_directory.end();
  assert(m_iter != m_end_iter);
  
  std::map<pheno_seq, size_t>::const_iterator  p_iter     = my_phenotype_directory.find(m_iter->second);
  std::map<pheno_seq, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
  assert(p_iter != p_end_iter);
  
  size_t  phenotype_index = p_iter->second;
  assert(phenotype_index < my_phenotypes.size());
  
  return  phenotype_index;
}

inline const em_phenotype&  
unrelated_em_phenotype_map::operator [](member m) const
{
  return  my_phenotypes[pheno_index(m)];
}

inline const em_phenotype&  
unrelated_em_phenotype_map::final_phenotype(member m) const
{
  assert(maximized);
  return  my_best_phenotypes[pheno_index(m)];
}

inline const FilteredMultipedigree&
unrelated_em_phenotype_map::multipedigree() const
{
  return  my_mped;
}


