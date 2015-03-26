//============================================================================
// File:      family.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/10/4 - created.                                   djb
//                                                                          
// Notes:      
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  family_generator
//============================================================================
//
inline
family_generator::family_generator(const vector<member>& members, 
                                     const locus_group& loci,
                                     cerrorstream& errors,
                                     ostream& messages,
                                     output_state& ostate)
      : my_members(members), my_loci(loci), my_errors(errors), my_messages(messages),
        my_output_state(ostate)
{}

inline cerrorstream&
family_generator::errors()
{
  return  my_errors;
}

inline ostream&
family_generator::messages()
{
  return  my_messages;
}

inline bool&
family_generator::messages_interrupted()
{
  return  my_output_state.msg_interrupted;
}

inline const locus_group&  
family_generator::loci() const
{
  return  my_loci;
}

inline const vector<member>&  
family_generator::members() const
{
  return  my_members;
}

inline family_generator::iterator
family_generator::begin()
{
  return  iterator(this);
}

inline family_generator::iterator
family_generator::end()
{
  return  iterator(this, true);
}


//============================================================================
// IMPLEMENTATION:  family_generator::iterator
//============================================================================
//
inline
family_generator::iterator::iterator(family_generator* generator, bool end)
     : my_generator(generator)
{
  const vector<member>&  members = my_generator->members();
  my_member_end_iter = members.end();

  if(end)
  {
    my_member_iter = my_member_end_iter;
  }
  else
  {
    my_member_iter = members.begin();
    if(my_member_iter != my_member_end_iter)
    {
      // - Shouldn't happen!
      //
      assert(MPED::mp_utilities::nuclear_family_count(*my_member_iter) != 0);

      my_subped_ptr = (*my_member_iter)->subpedigree();

      bool  build_combs_succeeded = false; 
      try
      {
        build_combs_succeeded = build_combs();
      }
      catch(const bad_alloc&)
      {
        memory_error();
        ++(*this);
        return;
      }
  
      if(! build_combs_succeeded)
      {
        recombination_error();
        ++(*this);
      }
    }
  }
}

inline family_generator::iterator
family_generator::iterator::operator ++()
{
  ++my_member_iter;
  if(my_member_iter != my_member_end_iter)
  {
  
    // - Shouldn't happen!
    //
    assert(MPED::mp_utilities::nuclear_family_count(*my_member_iter) != 0);
    
    my_subped_ptr = (*my_member_iter)->subpedigree();

    bool  build_combs_succeeded = false; 
    try
    {
      build_combs_succeeded = build_combs();
    }
    catch(const bad_alloc&)
    {
      memory_error();
      ++(*this);
      return  *this;
    }

    if(! build_combs_succeeded)
    {
      recombination_error();
      ++(*this);
    }    
  }
  
  return  *this;
}

inline void
family_generator::iterator::memory_error()
{
  bool&  interrupted = my_generator->messages_interrupted();
  if(! interrupted)
  {
    my_generator->messages() << endl;
    interrupted = true;
  }
  
  string  member_name   = (*my_member_iter)->name();
  string  pedigree_name = (*my_member_iter)->pedigree()->name();
  my_generator->errors() << priority(error) << "Not enough memory to process " 
                         << "constituent pedigree containing member "
                         << member_name << " in pedigree " << pedigree_name  
                         << ".  Skipping this constituent pedigree ..." << endl;
}

inline void
family_generator::iterator::recombination_error()
{
  bool&  interrupted = my_generator->messages_interrupted();
  if(! interrupted)
  {
    my_generator->messages() << endl;
    interrupted = true;
  }
  
  string  member_name   = (*my_member_iter)->name();
  string  pedigree_name = (*my_member_iter)->pedigree()->name();
  my_generator->errors() << priority(information) << "Constituent pedigree containing member "
                         << member_name << " in pedigree " << pedigree_name << " is not " 
                         << "consistent with an assumption of no recombination.  Skipping "
                         << "this constituent pedigree ..." << endl;
}

inline member
family_generator::iterator::get_member() const
{
  return  *my_member_iter;
}

inline family_generator::pedigree
family_generator::iterator::get_pedigree() const
{
  return  (*my_member_iter)->pedigree();
}

inline const vector<size_t>&  
family_generator::iterator::get_inconsistent_loci() const
{
  return  my_inconsistent_loci;
}

inline bool
family_generator::iterator::operator ==(const iterator& other) const
{
  return  my_generator == other.my_generator &&
          my_member_iter == other.my_member_iter;
}

inline bool
family_generator::iterator::operator !=(const iterator& other) const
{
  return  ! operator ==(other);
}

inline set<hap_seq_comb>&
family_generator::iterator::operator *()
{
  return  my_hap_seq_combs;
}

inline set<hap_seq_comb>*
family_generator::iterator::operator ->()
{
  return  &my_hap_seq_combs;
}

//============================================================================
// IMPLEMENTATION:  not_uninformative_leaf
//============================================================================
//
inline
not_uninformative_leaf::not_uninformative_leaf(
          FPED::has_informative_loci<FilteredMultipedigree::member_type>& base_filter)
      : genotyped(base_filter)
{}

// - Return true unless member is a leaf, ungenotyped and has informative siblings.
//   Will be retained.
//
inline bool
not_uninformative_leaf::operator ()(const FilteredMultipedigree::member_type& member) const
{
  bool  return_value;

  if(member.offspring_count() > 0)        // Not a 'leaf'.
  {
    return_value = true;
  }
  else if(genotyped(member))              // Informative.
  {
    return_value = true;
  }
  else if(member.sibling_count() == 0)    // Only child or singleton.
  {
    return_value = true;
  }
  else                                    // Look for informative siblings.
  {
    return_value = true;
  
    FilteredMultipedigree::sibling_const_iterator  sib_iter     = member.sibling_begin();
    FilteredMultipedigree::sibling_const_iterator  sib_end_iter = member.sibling_end();    
    for(; sib_iter != sib_end_iter; ++sib_iter)
    {
      if(genotyped(*sib_iter))
      {
        return_value = false;             // Member has informative sib.  Will not be retained.
        break;    
      }
    }
    
    if(return_value == true)              // There are multiple uninformative sibs.  Keep only one.
    {
      return_value = keeper(member);
    }    
  }

  return  return_value;
}

// - Arbitrarily pick one person from member's sibship to keep.  It must be the 
//   same person each time this is called.
//
inline bool
not_uninformative_leaf::keeper(const FilteredMultipedigree::member_type& member)
{
  set<const FilteredMultipedigree::member_type*>  sibship;
  
  FilteredMultipedigree::family_const_pointer  family = member.family(); 
  FilteredMultipedigree::offspring_const_iterator  o_iter     = family->offspring_begin();
  FilteredMultipedigree::offspring_const_iterator  o_end_iter = family->offspring_end();
  for(; o_iter != o_end_iter; ++o_iter)
  {
    sibship.insert(&(*o_iter));
  }

  return  &member == *(sibship.begin());
}

//============================================================================
// IMPLEMENTATION:  family_em_phenotype_map
//============================================================================
//
inline size_t
family_em_phenotype_map::pheno_index(member m) const
{
  std::map<member, size_t>::const_iterator  p_iter     = my_phenotype_directory.find(m);
  std::map<member, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
  assert(p_iter != p_end_iter);
  
  size_t  phenotype_index = p_iter->second;
  assert(phenotype_index < my_phenotypes.size());
  
  return  phenotype_index;
}

inline const em_phenotype&
family_em_phenotype_map::operator [](member m) const
{
  return  my_phenotypes[pheno_index(m)];
}

inline const em_phenotype&
family_em_phenotype_map::final_phenotype(member m) const
{
  assert(maximized);
  return  my_best_phenotypes[pheno_index(m)];
}

inline const FilteredMultipedigree&
family_em_phenotype_map::multipedigree() const
{
  return  my_filtered_mped;
}

template <typename T> inline void  
family_em_phenotype_map::related_sub_build(T iter, T end_iter)
{
  for(; iter != end_iter; ++iter)
  {
    try
    {
      set<hap_seq_comb>&  combs = *iter;

      assert(! combs.empty());
      my_total_hap_count += combs.begin()->size();
      my_phenotypes.push_back(em_phenotype(&my_haplotypes, combs));
      my_phenotype_directory.insert(make_pair(iter.get_member(), my_phenotypes.size() - 1));
    }
    catch(const bad_alloc&)
    {
      if(! my_output_state.msg_interrupted)
      {
        my_messages << endl;
        my_output_state.msg_interrupted = true;
      }
      
      my_errors << priority(fatal) << "Not enough memory available to process "
                << "the data.  Currently processing member " << iter.get_member()->name() 
                << " in pedigree " << iter.get_member()->pedigree()->name() 
                << ".  Exiting program ..." <<  endl;
      exit(1);
    }                
  }
}

