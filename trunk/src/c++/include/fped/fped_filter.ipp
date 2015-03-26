#ifndef FPED_FILTER_H
#include "fped/fped_filter.h"
#endif

namespace SAGE
{
namespace FPED
{

// ===========================================
// FilterResults inlines
// ===========================================

inline
FilterResults::FilterResults()
  : my_included_members(new MemberPtrList()),
    my_excluded_members(new MemberPtrList())
{ }

inline
FilterResults::FilterResults(const FilterResults& r)
  : my_included_members(r.my_included_members),
    my_excluded_members(r.my_excluded_members)
{ }

inline
FilterResults::~FilterResults()
{ }

inline
FilterResults& FilterResults::operator=(const FilterResults& r)
{
  if(&r != this)
  {
    my_included_members = r.my_included_members;
    my_excluded_members = r.my_excluded_members;
  }

  return *this;
}

inline
size_t FilterResults::get_included_member_count() const
{
  return my_included_members->size();
}

inline
const FilterResults::MemberPtrList& FilterResults::get_included_members() const
{
  return *my_included_members;
}

inline
FilterResults::MemberPtrIterator FilterResults::get_included_member_begin() const
{
  return my_included_members->begin();
}

inline
FilterResults::MemberPtrIterator FilterResults::get_included_member_end  () const
{
  return my_included_members->end();
}

inline
size_t FilterResults::get_excluded_member_count() const
{
  return my_excluded_members->size();
}

inline
const FilterResults::MemberPtrList& FilterResults::get_excluded_members() const
{
  return *my_excluded_members;
}

inline
FilterResults::MemberPtrIterator FilterResults::get_excluded_member_begin() const
{
  return my_excluded_members->begin();
}

inline
FilterResults::MemberPtrIterator FilterResults::get_excluded_member_end  () const
{
  return my_excluded_members->end();
}

inline void FilterResults::uniquify()
{
  if(!my_included_members.unique())
  {
    MemberPtrListShPtr tmp(new MemberPtrList());
    
    *tmp = *my_included_members;
    
    my_included_members = tmp;
  }

  if(!my_excluded_members.unique())
  {
    MemberPtrListShPtr tmp(new MemberPtrList());
    
    *tmp = *my_excluded_members;
    
    my_excluded_members = tmp;
  }
}

inline
void FilterResults::add_member_to_included(const MPED::member_base& m)
{
  uniquify();
  my_included_members->push_back(&m);
}

inline
void FilterResults::add_member_to_excluded(const MPED::member_base& m)
{
  uniquify();
  my_excluded_members->push_back(&m);
}

inline
void FilterResults::splice(FilterResults& r)
{
  if(&r == this) return;
  
  uniquify();

  if(my_included_members.get() != r.my_included_members.get())
    my_included_members->splice(  my_included_members->begin(), *r.my_included_members);
  else
    SAGE_internal_error();

  if(my_included_members.get() != r.my_included_members.get())
    my_excluded_members->splice(  my_excluded_members->begin(), *r.my_excluded_members);
  else
    SAGE_internal_error();
}

inline OUTPUT::Table 
FilterResults::get_excluded_member_table() const
{
  OUTPUT::Table t("Excluded member list");

  t << OUTPUT::TableColumn("Pedigree id") << OUTPUT::TableColumn("Member id");

  for(MemberPtrIterator i = get_excluded_member_begin(); i != get_excluded_member_end(); ++i)
    t << (OUTPUT::TableRow() << (*i)->pedigree()->name() << (*i)->name());

  return t;
}


// ===========================================
// MPFilterer Templates and inlines
// ===========================================

template<class MPTYPE>
FilterResults
MPFilterer::add_multipedigree(FilteredMultipedigree& fped, const MPTYPE& mped)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add each pedigree
  typename MPTYPE::pedigree_const_iterator ped = mped.pedigree_begin();
  for( ; ped != mped.pedigree_end(); ++ped )
  {
    FilterResults fr = add_pedigree(fped, *ped);

    filter_results.splice(fr);
  }
  
  return filter_results;
}

template <class MPTYPE, class FILTER>
inline FilterResults
MPFilterer::add_multipedigree_filtered
    (FilteredMultipedigree& fped,
     const MPTYPE&          mped,
     FILTER                 f     )
{
  if(f(mped))
    return add_multipedigree(fped, mped);
  else
  {
    FilterResults fr;
    
    for(size_t i = 0; i != mped.member_count(); ++i)
    {
      fr.add_member_to_excluded(mped.member_index(i));
    }
    
    return fr;
  }
}

template<class MPTYPE, class FILTER>
FilterResults 
MPFilterer::add_multipedigree_filtered_by_pedigrees
        (FilteredMultipedigree& fped,
         const MPTYPE&          mped,
         FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add each pedigree
  typename MPTYPE::pedigree_const_iterator ped = mped.pedigree_begin();
  for( ; ped != mped.pedigree_end(); ++ped )
  {
    FilterResults fr = add_pedigree_filtered(fped, *ped, f);

    filter_results.splice(fr);
  }
  
  return filter_results;
}

template<class MPTYPE, class FILTER>
FilterResults 
MPFilterer::add_multipedigree_filtered_by_subpedigrees
        (FilteredMultipedigree& fped,
         const MPTYPE&          mped,
         FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add each pedigree
  typename MPTYPE::pedigree_const_iterator ped = mped.pedigree_begin();
  for( ; ped != mped.pedigree_end(); ++ped )
  {
    FilterResults fr = add_pedigree_filtered_by_subpedigrees(fped, *ped, f);

    filter_results.splice(fr);
  }
  
  return filter_results;
}

template<class MPTYPE, class FILTER>
FilterResults 
MPFilterer::add_multipedigree_filtered_by_unconnecteds
        (FilteredMultipedigree& fped, 
         const MPTYPE&          mped,
         FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add each pedigree
  typename MPTYPE::pedigree_const_iterator ped = mped.pedigree_begin();
  for( ; ped != mped.pedigree_end(); ++ped )
  {
    FilterResults fr = add_pedigree_filtered_by_unconnecteds(fped, *ped, f);

    filter_results.splice(fr);
  }
  
  return filter_results;
}

template<class MPTYPE, class FILTER>
FilterResults 
MPFilterer::add_multipedigree_filtered_by_members
        (FilteredMultipedigree& fped,
         const MPTYPE&          mped,
         FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add each pedigree
  typename MPTYPE::pedigree_const_iterator ped = mped.pedigree_begin();
  for( ; ped != mped.pedigree_end(); ++ped )
  {
    FilterResults fr = add_pedigree_filtered_by_members(fped, *ped, f);

    filter_results.splice(fr);
  }
  
  return filter_results;
}
    

template<class PTYPE>
FilterResults 
MPFilterer::add_pedigree(FilteredMultipedigree& fped, const PTYPE& ped)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add subpedigrees
  typename PTYPE::subpedigree_const_iterator sped = ped.subpedigree_begin();
  for( ; sped != ped.subpedigree_end(); ++sped )
  {
    FilterResults fr = add_subpedigree(fped, *sped);
    
    filter_results.splice(fr);
  }
  
  // Add unconnecteds
  FilterResults uncon_results = add_members(fped, ped.unconnected_begin(), ped.unconnected_end());

  filter_results.splice(uncon_results);
  
  return filter_results;
}

template<class PTYPE, class FILTER>
FilterResults
MPFilterer::add_pedigree_filtered
    (FilteredMultipedigree& fped,
     const PTYPE&           ped,
     FILTER                 f)
{
  if(f(ped))
    return add_pedigree(fped, ped);
  else
  {
    FilterResults fr;
    
    for(size_t i = 0; i != ped.member_count(); ++i)
    {
      fr.add_member_to_excluded(ped.member_index(i));
    }
    
    return fr;
  }
  
}

template<class PTYPE, class FILTER>
FilterResults
MPFilterer::add_pedigree_filtered_by_subpedigrees
    (FilteredMultipedigree& fped,
     const PTYPE&           ped,
     FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add subpedigrees
  typename PTYPE::subpedigree_const_iterator sped = ped.subpedigree_begin();
  for( ; sped != ped.subpedigree_end(); ++sped )
  {
    FilterResults fr = add_subpedigree_filtered(fped, *sped, f);
    
    filter_results.splice(fr);
  }
  
  return filter_results;
}

template<class PTYPE, class FILTER>
FilterResults
MPFilterer::add_pedigree_filtered_by_unconnecteds
    (FilteredMultipedigree& fped,
     const PTYPE&           ped,
     FILTER                 f)
{
  FilterResults filter_results = add_members_filtered(fped, ped.unconnected_begin(), ped.unconnected_end(), f);

  return filter_results;
}

template<class PTYPE, class FILTER>
FilterResults
MPFilterer::add_pedigree_filtered_by_members
    (FilteredMultipedigree& fped,
     const PTYPE&           ped,
     FILTER                 f)
{
  // Keep track of who we add and remove
  FilterResults filter_results;

  // Add subpedigrees
  typename PTYPE::subpedigree_const_iterator sped = ped.subpedigree_begin();
  for( ; sped != ped.subpedigree_end(); ++sped )
  {
    FilterResults fr = add_subpedigree_filtered_by_members(fped, *sped, f);
    
    filter_results.splice(fr);
  }
  
  // Add unconnecteds
  FilterResults uncon_results = add_members_filtered(fped, ped.unconnected_begin(), ped.unconnected_end(), f);

  filter_results.splice(uncon_results);
  
  return filter_results;
}

template<class SPTYPE>
FilterResults 
MPFilterer::add_subpedigree
    (FilteredMultipedigree& fped,
     const SPTYPE&          sped)
{
  // Add our members
  FilterResults filter_results = add_members(fped, sped.member_begin(), sped.member_end());
  
  // Add our lineages.
  for(typename SPTYPE::member_const_iterator mem  = sped.member_begin(); 
                                             mem != sped.member_end(); 
                                           ++mem)
  {
    // If the individual is a founder, there's nothing to do
    if(mem->is_founder()) continue;

    // Get the parents
    const typename SPTYPE::member_type& p1 = *mem->parent1();
    const typename SPTYPE::member_type& p2 = *mem->parent2();
    
    fped.add_lineage(mem->pedigree()->name(), mem->name(), p1.name());
    fped.add_lineage(mem->pedigree()->name(), mem->name(), p2.name());
  }

  return filter_results;
}

template<class SPTYPE, class FILTER>
FilterResults 
MPFilterer::add_subpedigree_filtered
    (FilteredMultipedigree& fped, 
     const SPTYPE&          sped,
     FILTER                 f)
{
  if(f(sped))
    return add_subpedigree(fped, sped);
  else
  {
    FilterResults fr;
    
    for(size_t i = 0; i != sped.member_count(); ++i)
    {
      fr.add_member_to_excluded(sped.member_index(i));
    }
    
    return fr;
  }

}         

template<class SPTYPE, class FILTER>
FilterResults 
MPFilterer::add_subpedigree_filtered_by_members
    (FilteredMultipedigree& fped,
     const SPTYPE&          sped,
     FILTER                 f)
{
  // Add our members
  FilterResults filter_results = add_members_filtered(fped, sped.member_begin(), sped.member_end(), f);
  
  // Create a vector of bools where the value is true if the member with that
  // subindex was included, false otherwise.
  vector<bool> included_members(sped.member_count(), false);

  for(FilterResults::MemberPtrIterator i  = filter_results.get_included_member_begin();
                                       i != filter_results.get_included_member_end();
                                     ++i)
  {
    included_members[(*i)->subindex()] = true;
  }
  
  // Check our lineages.  We add a particular lineage if and only if the parents 
  // and the child of that lineage were all included.
  for(typename SPTYPE::member_const_iterator mem  = sped.member_begin(); 
                                             mem != sped.member_end(); 
                                           ++mem)
  {
    // If the member wasn't included, there's nothing to do
    if(!included_members[mem->subindex()]) continue;

    // If the individual is a founder, there's nothing to do either
    if(!mem->parent1()) continue;

    // Get the parents
    const typename SPTYPE::member_type& p1 = *mem->parent1();
    const typename SPTYPE::member_type& p2 = *mem->parent2();
    
    // If the parents are both true, add the lineages
    if(included_members[p1.subindex()] && included_members[p2.subindex()])
    {
      fped.add_lineage(mem->pedigree()->name(), mem->name(), p1.name());
      fped.add_lineage(mem->pedigree()->name(), mem->name(), p2.name());
    }
  }

  return filter_results;
}

template<class FTYPE>
FilterResults 
MPFilterer::add_family
    (FilteredMultipedigree& fped,
     const FTYPE&           fam)
{
  // Add parents
  FilterResults filter_results = add_members(fped, fam.parent_begin(), fam.parent_end());
  
  // Add offspring
  FilterResults offspring_filter_results =
      add_members(fped, fam.offspring_begin(), fam.offspring_end());
      
  // Get the parents
  const typename FTYPE::member_type& p1 = *fam.parent1();
  const typename FTYPE::member_type& p2 = *fam.parent2();
  
  // Iterate through included children
  for(typename FTYPE::offspring_const_iterator i  = fam.offspring_begin();
      i != fam.offspring_end();
      ++i)
  {
    fped.add_lineage(i->pedigree()->name(), i->name(), p1.name());
    fped.add_lineage(i->pedigree()->name(), i->name(), p2.name());
  }

  filter_results.splice(offspring_filter_results);

  return filter_results;
}             
  
template<class FTYPE, class FILTER>
FilterResults add_family_filtered
        (FilteredMultipedigree& fped,
         const FTYPE&           fam,
         FILTER                 f)
{
  if(f(fam))
  {
    return add_family(fped,fam);
  }
  else
  {
    FilterResults fr;
    
    fr.add_member_to_excluded(*fam.parent1());
    fr.add_member_to_excluded(*fam.parent2());
    
    for(typename FTYPE::offspring_const_iterator i = fam.offspring_begin(); i != fam.offspring_end(); ++i)
    {
      fr.add_member_to_excluded(*i);
    }
    
    return fr;
    
  }
}

template<class FTYPE, class FILTER>
FilterResults 
MPFilterer::add_family_filtered_by_members
    (FilteredMultipedigree& fped,
     const FTYPE&           fam,
     FILTER                 f)
{
  // Add parents
  FilterResults filter_results
      = add_members_filtered(fped, fam.parent_begin(), fam.parent_end(), f);
  
  // Determine if both parents were included.  This determines if we include
  // lineage later
  bool parents_both_included = (filter_results.get_excluded_member_count() == 0);

  // Add offspring
  FilterResults offspring_filter_results =
      add_members_filtered(fped, fam.offspring_begin(), fam.offspring_end(), f);
      
  // If both parents included, add lineages for any child included.
  if(parents_both_included)
  {
    // Get the parents
    const typename FTYPE::member_type& p1 = *fam.parent1();
    const typename FTYPE::member_type& p2 = *fam.parent2();
    
    // Iterate through included children
    for(FilterResults::MemberPtrList::const_iterator i  = offspring_filter_results.get_included_member_begin();
                                                     i != offspring_filter_results.get_included_member_end();
                                                   ++i)
    {
      fped.add_lineage((*i)->pedigree()->name(), (*i)->name(), p1.name());
      fped.add_lineage((*i)->pedigree()->name(), (*i)->name(), p2.name());
    }
  }

  filter_results.splice(offspring_filter_results);

  return filter_results;
}

template <class MTYPE>
FilterResults
MPFilterer::add_member
    (FilteredMultipedigree& fped,
     const MTYPE&           mem)
{
  FilterResults filter_results;

  // Add the member, referencing the original source
  fped.add_member(mem.pedigree().name(), mem.name(), mem.get_detailed_sex(), FilteredMemberInfo(mem));

  // Update the included and excluded lists
  filter_results.add_member_to_included(mem);
  
  return filter_results;
}

template <class MTYPE, class FILTER>
FilterResults
MPFilterer::add_member_filtered
    (FilteredMultipedigree& fped,
     const MTYPE&           mem,
     FILTER                 f     )
{
  FilterResults filter_results;

  // If the member is to be included, the function f will return true.
  if(f(mem))
  {
    // Add the member, referencing the original source
    fped.add_member(mem.pedigree().name(), mem.name(), mem.get_detailed_sex(), FilteredMemberInfo(mem));

    // Update the included and excluded lists
    filter_results.add_member_to_included(mem);
  }
  else
  {
    // Don't add the member
    
    // Update the included and excluded lists
    filter_results.add_member_to_excluded(mem);
  }
  
  return filter_results;
}

template<class Iterator>
FilterResults
MPFilterer::add_members
  (FilteredMultipedigree& fped, Iterator b, Iterator e)
{
  FilterResults filter_results;
  
  // Iterate from begin to end, adding each member and keeping track of which ones
  // were included and which excluded
  for( ; b != e; ++b )
  {
    fped.add_member(b->pedigree()->name(), b->name(), b->get_detailed_sex(), FilteredMemberInfo(*b));
    filter_results.add_member_to_included(*b);
  }
  
  return filter_results;
}

template<class Iterator, class FILTER>
FilterResults
MPFilterer::add_members_filtered
  (FilteredMultipedigree& fped, Iterator b, Iterator e, FILTER& f)
{
  FilterResults filter_results;
  
  // Iterate from begin to end, adding each member and keeping track of which ones
  // were included and which excluded
  for( ; b != e; ++b )
  {
    if(f(*b))
    {
      fped.add_member(b->pedigree()->name(), b->name(), b->get_detailed_sex(), FilteredMemberInfo(*b));
      filter_results.add_member_to_included(*b);
    }
    else
    {
      filter_results.add_member_to_excluded(*b);
    }
  }
  
  return filter_results;
}

} // End namespace RPED
} // End namespace SAGE

