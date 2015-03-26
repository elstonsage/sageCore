//============================================================================
// File:      founders.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/11/5 - created.                                   djb
//                                                                          
// Notes:      
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  founders_generator
//============================================================================
//
inline
founders_generator::founders_generator(const vector<member>& members, 
                                     const locus_group& loci,
                                     cerrorstream& errors,
                                     ostream& messages,
                                     output_state& ostate)
      : my_members(members), my_loci(loci), my_errors(errors), my_messages(messages),
        my_output_state(ostate)
{}

inline cerrorstream&
founders_generator::errors()
{
  return  my_errors;
}

inline ostream&
founders_generator::messages()
{
  return  my_messages;
}

inline bool&
founders_generator::messages_interrupted()
{
  return  my_output_state.msg_interrupted;
}

inline const locus_group&  
founders_generator::loci() const
{
  return  my_loci;
}

inline const vector<member>&  
founders_generator::members() const
{
  return  my_members;
}

inline founders_generator::iterator
founders_generator::begin()
{
  return  iterator(this);
}

inline founders_generator::iterator
founders_generator::end()
{
  return  iterator(this, true);
}


//============================================================================
// IMPLEMENTATION:  founders_generator::iterator
//============================================================================
//
inline
founders_generator::iterator::iterator(founders_generator* generator, bool end)
     : my_hap_count(0), my_generator(generator)
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
      // - Singleton; shouldn't happen!
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

inline founders_generator::iterator
founders_generator::iterator::operator ++()
{
  ++my_member_iter;
  
  if(my_member_iter != my_member_end_iter)
  {
    // - Singleton; shouldn't happen!
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
founders_generator::iterator::memory_error()
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
founders_generator::iterator::recombination_error()
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
founders_generator::iterator::get_member() const
{
  return  *my_member_iter;
}

inline founders_generator::pedigree
founders_generator::iterator::get_pedigree() const
{
  return  (*my_member_iter)->pedigree();
}

inline const vector<size_t>&  
founders_generator::iterator::get_inconsistent_loci() const
{
  return  my_inconsistent_loci;
}

inline bool
founders_generator::iterator::operator ==(const iterator& other) const
{
  return  my_generator == other.my_generator &&
          my_member_iter == other.my_member_iter;
}

inline bool
founders_generator::iterator::operator !=(const iterator& other) const
{
  return  ! operator ==(other);
}

inline set<hap_seq_comb>&
founders_generator::iterator::operator *()
{
  return  my_hap_seq_combs;
}

inline set<hap_seq_comb>*
founders_generator::iterator::operator ->()
{
  return  &my_hap_seq_combs;
}
 
 
