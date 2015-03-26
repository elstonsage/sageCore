#ifndef COMMON_ALLELE_SET_H
#include "mcmc/common_allele_set.h"
#endif

namespace SAGE
{
namespace MCMC
{

inline
CommonAlleleSet::CommonAlleleSet()
  : my_model(),
    my_state(EMPTY),
    my_allele_genotype_counts(0, 0),
    my_allele_pair_count(0),
    my_is_dirty(false)
{
  set_to_empty();
}

inline
CommonAlleleSet::CommonAlleleSet(const MLOCUS::genotype_model& gmodel)
  : my_model(gmodel),
    my_state(EMPTY),
    my_allele_genotype_counts(my_model.allele_count(), 0),
    my_allele_pair_count(0),
    my_is_dirty(false)
{
  set_to_empty();
}

inline
CommonAlleleSet::CommonAlleleSet(const CommonAlleleSet& cas)
  : my_model(cas.my_model),
    my_state(cas.my_state),
    my_allele_genotype_counts(cas.my_allele_genotype_counts),
    my_allele_pair_count(cas.my_allele_pair_count),
    my_is_dirty(cas.my_is_dirty)
{
  my_valid_alleles[0] = cas.my_valid_alleles[0];
  my_valid_alleles[1] = cas.my_valid_alleles[1];
}

inline
CommonAlleleSet& 
CommonAlleleSet::operator=(const CommonAlleleSet& cas)
{
  if(this != &cas)
  {
    my_model = cas.my_model;
    my_state = cas.my_state;
  
    my_valid_alleles[0] = cas.my_valid_alleles[0];
    my_valid_alleles[1] = cas.my_valid_alleles[1];
  
    my_allele_genotype_counts = cas.my_allele_genotype_counts;
    
    my_allele_pair_count = cas.my_allele_pair_count;
  
    my_is_dirty = cas.my_is_dirty;
  }

  return *this;
}

inline
size_t 
CommonAlleleSet::get_max_allele_count() const
{
  return my_allele_genotype_counts.size();
}

inline
CommonAlleleSet::StateEnum 
CommonAlleleSet::get_status() const
{
  if(my_is_dirty)
  {
    find_and_set_valid_alleles();
  }

  return my_state;
}

inline
CommonAlleleSet::AlleleID 
CommonAlleleSet::get_allele_id(size_t a) const
{
  if(my_is_dirty)
  {
    find_and_set_valid_alleles();
  }

  return my_valid_alleles[a];
}

inline
size_t 
CommonAlleleSet::get_allele_pair_count() const
{
  return my_allele_pair_count;
}

inline
size_t 
CommonAlleleSet::get_allele_pair_count_for_allele(AlleleID a) const
{
  return my_allele_genotype_counts[a.id()];
}

inline
bool 
CommonAlleleSet::is_dirty() const
{
  return my_is_dirty;
}

inline
bool 
CommonAlleleSet::add_allele_pair(AlleleID a1, AlleleID a2)
{
  // We first check for basic badness.  If the alleles are invalid or too large,
  // things are not good

  if(!a1.is_valid() || !a2.is_valid()      ||
     get_max_allele_count() <= a1.id() ||
     get_max_allele_count() <= a2.id()     )
  {
    SAGE_internal_error();
  }

  // Store the original state for later comparison
  
  StateEnum original_state = my_state;
  
  // We increment the allele genotype counts and the allele pair count.
  
  ++my_allele_pair_count;
  
  ++my_allele_genotype_counts[a1.id()];
  
  // We only increment the alleles once, even if there are two copies (homozygous)

  if(a1 != a2)
  {
    ++my_allele_genotype_counts[a2.id()];
  }
  
  // If we're already potentially changed, we can just return here.  We
  // assume that the status will be updated properly when requested.
  
  if(my_is_dirty) return true;
  
  // Depending on state, we do different things
  
  switch(my_state)
  {
    case EMPTY     :
      // If empty, we know that a1 and a2 are ok

      if(a1 != a2)
        set_to_two_valid(a1, a2);
      else
        set_to_one_valid(a1);

      break;

    case TWO_VALID :
    case ONE_VALID :
      // For both two and one valid, we have to check the alleles for
      // continuing validity
      determine_and_set_valid_alleles(a1,a2);
      break;

    case INVALID   :
      // In the invalid case, adding a pair will never make anything
      // become valid.
      break;
  }

  // Return if the state changed.  Note that since we're not dirty, we have
  // actually changed the state, but also the underlying variables, so we don't
  // need to set the dirty flag.
  
  return (my_state != original_state);
}

inline
bool 
CommonAlleleSet::remove_allele_pair(AlleleID a1, AlleleID a2)
{
  // We first check for basic badness.  If the alleles are invalid or too large,
  // things are not good

  if(!a1.is_valid() || !a2.is_valid()      ||
     get_max_allele_count() <= a1.id() ||
     get_max_allele_count() <= a2.id()     )
  {
    SAGE_internal_error();
  }

  // Store the original state for later comparison
  
  StateEnum original_state = my_state;
  
  // We decrement the allele genotype counts and the allele pair count.

  --my_allele_pair_count;
  
  --my_allele_genotype_counts[a1.id()];
  
  // We only decrement the alleles once, even if there are two copies (homozygous)

  if(a1 != a2)
  {
    --my_allele_genotype_counts[a2.id()];
  }
  
  // If we're dirty, we only check for no alleles, which means we're empty. 
  // This is to save time rather than generating state info when it's often
  // not required.

  if(is_dirty())
  {
    if(my_allele_pair_count == 0)
    {
      set_to_empty();
      my_is_dirty=false;
    }
  }
  else
  {
    // Because we're not dirty, we know our state is accurate, so we can
    // do a lot more with some cases.

    // Depending on state, we do different things
    
    switch(my_state)
    {
      case EMPTY     :
        // If empty, we shouldn't be doing a remove at all, so this is a major error
        SAGE_internal_error();
        break;
        
      case TWO_VALID :
        // Unless the number of allele pairs is now 0, the two alleles that were
        // valid will continue to be valid after the removal, so we only change
        // if we're now empty
        
        if(my_allele_pair_count == 0)
          set_to_empty();
        break;

      case ONE_VALID :
        // Unless the number of allele pairs is now 0, the allele we have now
        // will remain valid.  However, another allele might have also become valid.
        // So we check for empty status.  IF it's not empty, we don't know
        // what the valid alleles are, so we set the dirty flag.

        if(my_allele_pair_count == 0)
          set_to_empty();
        else
          my_is_dirty = true;
        break;
        
      case INVALID   :
        // In the invalid cases, we may have alleles that were
        // not valid which now are (due to decrease in allele pair count).
        // Since we don't know the alleles, we set the dirty flag
        
        my_is_dirty = true;
        break;
    }
  }

  // Return if the state changed or the dirty flag has been set
  
  return (my_state != original_state) || my_is_dirty;
}

inline
bool 
CommonAlleleSet::operator==(const CommonAlleleSet& cas) const
{
  return get_allele_id(0) == cas.get_allele_id(0) &&
         get_allele_id(1) == cas.get_allele_id(1);
}

inline
bool 
CommonAlleleSet::operator!=(const CommonAlleleSet& cas) const
{
  return get_allele_id(0) != cas.get_allele_id(0) ||
         get_allele_id(1) != cas.get_allele_id(1);
}

inline
void 
CommonAlleleSet::determine_and_set_valid_alleles(AlleleID a1, AlleleID a2)
{
  // Determine which of these alleles (if any) we'll keep.  We only keep
  // those alleles which have been in every allele pair added so far.

  bool keep_a1 = (my_allele_genotype_counts[a1.id()] == get_allele_pair_count());
  bool keep_a2 = (my_allele_genotype_counts[a2.id()] == get_allele_pair_count());

  // Determine state and set valid alleles.  This is done based on the 
  // joint keep status of the two alleles.
  
  if(keep_a1)
  {
    // We're keeping a1.  Determine if we should keep a2 as well
    // We keep both if keep_a2 is true and they're not the same allele

    if(keep_a2 && a1 != a2)
    {
      // Keeping both a1 and a2
      
      set_to_two_valid(a1, a2);
    }
    else
    {
      // Keeping only a1
      
      set_to_one_valid(a1);
    }
  }
  else
  {
     // We're not keeping a1, but we might keep a2
     if(keep_a2)
     {
       // Keeping only a2
       
       set_to_one_valid(a2);
     }
     else
     {
       // Keeping nothing

       set_to_invalid();
     }
  }
}

inline
void 
CommonAlleleSet::find_and_set_valid_alleles() const
{
  // Iterate through all the alleles, looking for alleles which have been in every
  // genotype pair.

  // We wish to keep a count of the alleles we've found.
  
  size_t valid_found = 0;
  
  // We only iterate until we've found at most 2 valid alleles, since there
  // can't be more than that.
  
  for(MLOCUS::allele_iterator i = my_model.allele_begin();
      valid_found < 2 && i != my_model.allele_end(); ++i)
  {
    // If the current allele i has been in every allele pair, it is valid
    
    if(get_allele_pair_count_for_allele(*i) == get_allele_pair_count())
    {
      my_valid_alleles[valid_found] = *i;
      
      ++valid_found;
    }
  }

  // At this point, we've found all the valid alleles, and know how many.  We
  // have to clear out the rest of the my_valid_alleles structure, though

  for(size_t i = valid_found; i < 2; ++i)
    my_valid_alleles[i] = MLOCUS::allele();

  // Wrapping up, we need to set our state.  The choice of values in the
  // StateEnum makes this easy, as it's converted directly from the number of
  // valid alleles
  
  my_state = (StateEnum) valid_found;

  my_is_dirty = false;
}

inline
void 
CommonAlleleSet::set_to_empty()
{
  my_state = EMPTY;

  my_valid_alleles[0] = my_valid_alleles[1] = MLOCUS::allele();
}

inline
void 
CommonAlleleSet::set_to_two_valid(AlleleID a1, AlleleID a2)
{
  my_state = TWO_VALID;

  // We want valid allele 0 to be the smaller, so we sort our source alleles.
  
  if(a2 < a1)
    std::swap(a1,a2);

  my_valid_alleles[0] = a1;
  my_valid_alleles[1] = a2;
}

inline
void 
CommonAlleleSet::set_to_one_valid(AlleleID a)
{
  my_state = ONE_VALID;

  my_valid_alleles[0] = a;
  my_valid_alleles[1] = MLOCUS::allele();
}

inline
void 
CommonAlleleSet::set_to_invalid()
{
  my_state = INVALID;

  my_valid_alleles[0] = my_valid_alleles[1] = MLOCUS::allele();
}

} // End MCMC
} // End SAGE
