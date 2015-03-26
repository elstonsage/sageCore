#ifndef FPED_FUNC_H
#include "fped/fped_func.h"
#endif

namespace SAGE
{
namespace FPED
{

// =============================================
// always_keep_member inlines
// =============================================
  
template <class TYPE>
inline 
bool always_keep::operator()(const TYPE&) const
{
  return true;
}

// =============================================
// has_informative_loci inlines
// =============================================

template <class MTYPE>
inline 
has_informative_loci<MTYPE>::has_informative_loci
   (const typename MTYPE::multipedigree_type& mped, 
    bool use_all)
  : my_loci_check_states(mped.info().markers().size(), use_all)
{ }

template <class MTYPE>
inline 
has_informative_loci<MTYPE>::has_informative_loci(const has_informative_loci& h)
  : my_loci_check_states(h.my_loci_check_states)
{ }

template <class MTYPE>
inline 
has_informative_loci<MTYPE>::~has_informative_loci()
{ }

template <class MTYPE>
inline has_informative_loci<MTYPE>& 
  has_informative_loci<MTYPE>::operator=(const has_informative_loci& h)
{
  if(&h == *this) return *this;
  
  my_loci_check_states = h.my_loci_check_states;

  return *this;
}

template <class MTYPE>
inline void 
  has_informative_loci<MTYPE>::set_check_status_for_locus(size_t locus_index, bool check)
{
  my_loci_check_states[locus_index] = check;
}

template <class MTYPE>
inline bool 
  has_informative_loci<MTYPE>::get_check_status_for_locus(size_t locus_index) const
{
  return my_loci_check_states[locus_index];
}

template <class MTYPE>
bool 
has_informative_loci<MTYPE>::operator()(const MTYPE& member) const
{
  // Get our multipedigree and pedigree info from the member
  
  const typename MTYPE::mpinfo_type&  mped_info = member.multipedigree()->info();
  const typename MTYPE::pedinfo_type& ped_info  = member.pedigree()->info();

  // Get our member index
  
  size_t mindex = member.index();

  // Check each loci for which the check status flag is true for existing phenotype.
  // If we find one, stop.
  
  bool pheno_exist = false;

  for( size_t m = 0; !pheno_exist && m < mped_info.marker_count(); ++m )
  {
    if(my_loci_check_states[m])
    {
      const RPED::RefMarkerInfo& mi = mped_info.marker_info(m);

      if( !ped_info.phenotype_missing(mindex, m, mi) )
        pheno_exist = true;
    }
  }

  // Return whether we've found one
  
  return pheno_exist;
}

// =============================================
// has_informative_traits inlines
// =============================================

template <class MTYPE>
inline 
has_informative_traits<MTYPE>::has_informative_traits
   (const typename MTYPE::multipedigree_type& mped, 
    bool use_all)
  : my_trait_check_states(mped.info().markers().size(), use_all)
{ }

template <class MTYPE>
inline 
has_informative_traits<MTYPE>::has_informative_traits(const has_informative_traits& h)
  : my_trait_check_states(h.my_trait_check_states)
{ }

template <class MTYPE>
inline 
has_informative_traits<MTYPE>::~has_informative_traits()
{ }

template <class MTYPE>
inline has_informative_traits<MTYPE>& 
  has_informative_traits<MTYPE>::operator=(const has_informative_traits& h)
{
  if(&h == *this) return *this;
  
  my_trait_check_states = h.my_trait_check_states;

  return *this;
}

template <class MTYPE>
inline void 
  has_informative_traits<MTYPE>::set_check_status_for_trait(size_t trait_index, bool check)
{
  my_trait_check_states[trait_index] = check;
}

template <class MTYPE>
inline bool 
  has_informative_traits<MTYPE>::get_check_status_for_trait(size_t trait_index) const
{
  return my_trait_check_states[trait_index];
}

template <class MTYPE>
bool 
has_informative_traits<MTYPE>::operator()(const MTYPE& member) const
{
  // Get our multipedigree and pedigree info from the member
  
  const typename MTYPE::mpinfo_type&  mped_info = member.multipedigree()->info();
  const typename MTYPE::pedinfo_type& ped_info  = member.pedigree()->info();

  // Get our member index
  
  size_t mindex = member.index();

  // Check each trait for which the check status flag is true for existing trait value.
  // If we find one, stop.
  
  bool trait_exist = false;

  for( size_t m = 0; !trait_exist && m < mped_info.trait_count(); ++m )
  {
    if(my_trait_check_states[m])
    {
      if( !ped_info.trait_missing(mindex, m) )
        trait_exist = true;
    }
  }

  // Return whether we've found one
  
  return trait_exist;
}

// =============================================
// is_inf_within_sped_t inlines
// =============================================

template<class MTYPE, class FILTER> 
inline
  is_inf_within_sped_t<MTYPE,FILTER>::is_inf_within_sped_t
    (const FILTER& f)
  : my_filter(f),
    my_current_subpedigree(NULL),
    my_informative_states()
{ }

template<class MTYPE, class FILTER> 
inline
  is_inf_within_sped_t<MTYPE,FILTER>::is_inf_within_sped_t
    (const is_inf_within_sped_t& src)
  : my_filter(src.my_filter),
    my_current_subpedigree(src.my_current_subpedigree),
    my_informative_states(src.my_informative_states)
{ }

template<class MTYPE, class FILTER> 
inline
  is_inf_within_sped_t<MTYPE,FILTER>::~is_inf_within_sped_t()
{ }

template<class MTYPE, class FILTER> 
  is_inf_within_sped_t<MTYPE,FILTER>& is_inf_within_sped_t<MTYPE,FILTER>::operator=
    (const is_inf_within_sped_t& src)
{
  if(&src == this) return *this;
  
  my_filter = src.filter;
  
  my_current_subpedigree = src.my_current_subpedigree;
  
  my_informative_states = src.my_informative_states;
}

template<class MTYPE, class FILTER> 
bool is_inf_within_sped_t<MTYPE,FILTER>::operator()
  (const MTYPE& member) const
{
  // IF the member has no lineage information, it's informativity reduces to
  // the my_filter function
  if(member.is_unconnected())
    return my_filter(member);
  
  if(member.subpedigree() != my_current_subpedigree)
  {
    my_current_subpedigree = member.subpedigree();
    calc_subped_inf_vector();
  }

  return my_informative_states[member.subindex()];
}


template<class MTYPE, class FILTER> 
pair<typename is_inf_within_sped_t<MTYPE,FILTER>::InfChildStateEnum,
     const MPED::member_base*>
  is_inf_within_sped_t<MTYPE,FILTER>::classify_inf_children
    (const MPED::member_base&  mem) const
{
  // Keep track of any informative child we see
  
  const MPED::member_base* informative_child = NULL;

  // Iterate through progeny, keeping track of the number of informative
  // children we've found.  If we find more than one, we can stop early.

  InfChildStateEnum inf_child_state = NONE;
  
  MPED::progeny_base_const_iterator off = mem.progeny_begin();

  for( ; inf_child_state < TWO_OR_MORE && off != mem.progeny_end(); ++off )
  {
    if( my_informative_states[off->subindex()] )
    {
      inf_child_state = (InfChildStateEnum) ((size_t) inf_child_state + 1);

      informative_child = &*off;
    }
  }

  return std::make_pair(inf_child_state, informative_child);
}

template<class MTYPE, class FILTER> 
void is_inf_within_sped_t<MTYPE,FILTER>::calc_subped_inf_vector() const
{
  // Initial setup:
  
  //   Assume everyone is informative

  my_informative_states.clear();
  my_informative_states.resize(my_current_subpedigree->member_count(), true);

  //  Iterate through the members to classify them based upon initial informativity
  //  through my_filter.  Uninformative individuals are not yet known to be
  //  informative or not (members informative by my_filter always are), so are
  //  kept track of as unclassified
  
  set<size_t>    unclassified_inds;

  vector<bool>   initial_informativity(my_informative_states);
  
  typename MTYPE::member_const_iterator ind = my_current_subpedigree->member_begin();
  for( ; ind != my_current_subpedigree->member_end(); ++ind )
  {
    initial_informativity[ind->subindex()] = my_filter(*ind);
    
    if( ! initial_informativity[ind->subindex()] )
      unclassified_inds.insert(ind->subindex());
  }

  // First Pass:
  
  //   Label unclassified individuals with no children as uninformative.
  //   Any such individuals cannot be structurally informative
  
  for( set<size_t>::iterator i  = unclassified_inds.begin();
                             i != unclassified_inds.end();
                           ++i )
  {
    const MPED::member_base& mem = my_current_subpedigree->member_index(*i);

    if( mem.progeny_begin() == mem.progeny_end() )
      my_informative_states[mem.subindex()] = false;
  }

  // Second Pass:
  
  //   Loop through unclassified individuals, looking for individuals which
  //   have insufficient individuals close to them which are informative to
  //   be structurally informative.  Repeat loop until nothin changes over
  //   a pass or all the unclassified individuals have been classified.
  
  //   Create a flag that indicates whether the last iteration of the while
  //   loop changed anyone's status.
  
  bool changed = true;

  //   Evaluate our state.  As long as there are individuals in the unclassified list,
  //   and we're continuing to change status, we keep iterating through our individuals.
  
  while( unclassified_inds.size() && changed )
  {
    // Initially, assume we won't change anything.
    
    changed = false;

    // Check each individual in the unclassified list for being uninformative
    
    set<size_t>::iterator ui = unclassified_inds.begin();
    for( ; ui != unclassified_inds.end(); ++ui )
    {
      // Unclassified individuals are labeled informative unless we can prove otherwise
      // If we've already shown a member to be uninformative, we don't have to
      // prove it again. But, if the unclassified individual is currently labeled 
      // informative, we have to check it to see whether it's local state clarrifies
      // its state
      if( my_informative_states[*ui] )
      {
        const MPED::member_base& mem = my_current_subpedigree->member_index(*ui);

        // Get informative child status
        
        pair<InfChildStateEnum, const MPED::member_base*>
          inf_child_state = classify_inf_children(mem);

        // If the member has two or more informative children, it is structurally
        // informative (at least for now), so do nothing more about it until
        // the next pass.
        
        if( inf_child_state.first == TWO_OR_MORE )
          continue;

        // This member may be uninformative.  Determine if it has 0 or 1 informative
        // children
        
        if( inf_child_state.first == NONE )
        {
          // If we have no informative children, we can label this as uninformative.
        
          my_informative_states[*ui] = false;
        }
        else
        {
          // The member must have exactly one informative child.  In this case, we
          // check the parents informativity to see if it's needed for ancestral
          // connection
  
          bool has_informative_parents = mem.parent1() && 
              ( my_informative_states[mem.parent1()->subindex()] || 
                my_informative_states[mem.parent2()->subindex()] );
  
          if( has_informative_parents )
            continue;
          
          // At this point, the member has no informative parents, and only one
          // informative child.  The only other place it can get validated as
          // informative from is the mate.
          
          const MPED::member_base* mate = inf_child_state.second->parent1();
  
          if( mate->index() == mem.index() )
            mate = inf_child_state.second->parent2();
          
          // If the mate is initially informative, the member is structurally informative
          
          if( initial_informativity[mate->subindex()] )
            continue;
          
          // Otherwise, if the mate's parents are informative, the member is
          // structurally informative
          
          bool mate_has_informative_parents = mate->parent1() && 
              ( my_informative_states[mate->parent1()->subindex()] || 
                my_informative_states[mate->parent2()->subindex()] );
          
          if( mate_has_informative_parents )
            continue;
          
          // Finally, if the mate has more than one informative child (through
          // a second spouse), the member is structurally informative 
  
          inf_child_state = classify_inf_children(*mate);
  
          if( inf_child_state.first == TWO_OR_MORE )
            continue;
          
          // We've exhausted the possibilities.  This individual, and its mate,
          // have no informative parents and only one informative child.
          // Therefore, they connect nothing informative, and can be classified
          // uninformative
          
          my_informative_states[mem.subindex()]   = false;
          my_informative_states[mate->subindex()] = false;
        }
      }
      
      // If the individual is now classified as uninformative, remove them
      // from the unclassified list and indicate that something has changed status
      if( !my_informative_states[*ui] )
      {
        changed = true;

        // Remove the unclassified individual.
        
        set<size_t>::iterator temp = ui;

        --temp;

        unclassified_inds.erase(ui);

        ui = temp;
      }
    }
  }
}
template<class FILTER>
inline is_inf_within_sped_t<typename FILTER::argument_type, FILTER>
    is_inf_within_sped(const FILTER& f)
{ return is_inf_within_sped_t<typename FILTER::argument_type, FILTER>(f); }

// =============================================
// is_family_sib_pair_inf_t inlines
// =============================================

template<class MTYPE, class FILTER> 
inline
  is_family_sib_pair_inf_t<MTYPE,FILTER>::is_family_sib_pair_inf_t
    (const typename MTYPE::family_type& fam, const FILTER& f)
  : my_filter(f),
    my_current_family(NULL)
{
  set_family(fam);
}

template<class MTYPE, class FILTER> 
inline
  is_family_sib_pair_inf_t<MTYPE,FILTER>::is_family_sib_pair_inf_t
    (const is_family_sib_pair_inf_t& src)
  : my_filter(src.my_filter),
    my_current_family(src.my_current_family),
    my_family_has_informative_spairs(src.my_family_has_informative_spairs)
{ }

template<class MTYPE, class FILTER> 
inline
  is_family_sib_pair_inf_t<MTYPE,FILTER>::~is_family_sib_pair_inf_t()
{
  my_current_family = NULL;
}

template<class MTYPE, class FILTER> 
  is_family_sib_pair_inf_t<MTYPE,FILTER>& is_family_sib_pair_inf_t<MTYPE,FILTER>::operator=
    (const is_family_sib_pair_inf_t& src)
{
  if(&src == this) return *this;
  
  my_filter = src.filter;
  
  my_current_family = src.my_current_family;
  
  my_family_has_informative_spairs = src.my_family_has_informative_spairs;
}

template<class MTYPE, class FILTER> 
void is_family_sib_pair_inf_t<MTYPE,FILTER>::set_family(const typename MTYPE::family_type& fam)
{
  if(my_current_family == &fam) return;
  
  my_current_family = &fam;
  
  // Set the my_family_has_inf_spairs flag based on 
  
  my_family_has_informative_spairs = does_fam_have_enough_inf_sibs();
}

template<class MTYPE, class FILTER> 
inline
//lint -e{10}  "Expecting an identifier" on the MTYPE::family_type is due
//             to the templating confusing lint.  Spurious
const typename MTYPE::family_type& is_family_sib_pair_inf_t<MTYPE,FILTER>::get_family() const
{
  return *my_current_family;
}

template<class MTYPE, class FILTER> 
bool is_family_sib_pair_inf_t<MTYPE,FILTER>::operator()
  (const MTYPE& member) const
{
  // IF the family wasn't informative, it doesn't matter who the person was
  
  if(!my_family_has_informative_spairs) return false;

  // The family has at least one valid sibling pair.  We must classify our member
  // into either a parent, a child, or a non-member.
  
  // IF the member is a parent, it's always going to be informative
  
  if(is_member_a_parent(member)) return true;
  
  // Otherwise, if the member is a child, we have to call the informativity
  // function on it to check
  
  if(is_member_a_child(member)) return my_filter(member);
  
  // Otherwise, the member doesn't belong to the family, so it's not informative
  // at all
  
  return false;
}

template<class MTYPE, class FILTER> 
inline
bool is_family_sib_pair_inf_t<MTYPE,FILTER>::is_member_a_parent(const MTYPE& mem) const
{
  return (&mem == my_current_family->parent1() ||
          &mem == my_current_family->parent2());
}

template<class MTYPE, class FILTER> 
inline
bool is_family_sib_pair_inf_t<MTYPE,FILTER>::is_member_a_child(const MTYPE& mem) const
{
  return (mem.family() == my_current_family); 
}

template<class MTYPE, class FILTER> 
bool is_family_sib_pair_inf_t<MTYPE,FILTER>::does_fam_have_enough_inf_sibs() const
{
  // Count the number of informative sibs.  Stop at 2, since we don't care once
  // we have at least two.
  
  size_t inf_sib_count = 0;

  typename MTYPE::offspring_const_iterator off = my_current_family->offspring_begin();
  
  for( ; inf_sib_count < 2 && off != my_current_family->offspring_end(); ++off)
  {
    if(my_filter(*off)) ++ inf_sib_count;
  }
  
  // if we have at least two, return true.
  
  return 2 <= inf_sib_count; 
}

template<class FTYPE, class FILTER>
inline is_family_sib_pair_inf_t<typename FTYPE::member_type, FILTER> 
    is_family_sib_pair_inf(const FTYPE& fam, const FILTER& f)
{
  return is_family_sib_pair_inf_t<typename FTYPE::member_type, FILTER>(fam, f);
}
// =====================
// filter_to_sib_pair()
// =====================

template <class MTYPE>
void filter_to_sib_pair
    (FilteredMultipedigree& fped,
     const MTYPE&           ind1,
     const MTYPE&           ind2,
     string                 ped_name)
{
  // If the ped_name isn't provided, we must create it
  if(ped_name.size() == 0)
  {
    ped_name = ind1.pedigree()->name() + ":" + ind1.name() + " x " +
               ind2.pedigree()->name() + ":" + ind2.name();
  }
  
  // Create two dummy parents
  fped.add_member(ped_name, "~dummy2", SAGE::MPED::SEX_MALE);
  fped.add_member(ped_name, "~dummy1", SAGE::MPED::SEX_FEMALE);
  
  // add ind1 & ind2 as offsprings
  fped.add_member(ped_name, ind1.name(), ind1.get_detailed_sex(), FilteredMemberInfo(ind1));
  fped.add_lineage(ped_name, ind1.name(), "~dummy1", "~dummy2");
  
  fped.add_member(ped_name, ind2.name(), ind2.get_detailed_sex(), FilteredMemberInfo(ind2));
  fped.add_lineage(ped_name, ind2.name(), "~dummy1", "~dummy2");
}


} // End namespace RPED
} // End namespace SAGE

