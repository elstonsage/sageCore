#include "lvec/inheritance_vector.h"

namespace SAGE
{

inheritance_vector& inheritance_vector::operator=
    (const inheritance_vector& i)
{
  mm    = i.mm;
  found = i.found;
  nonf  = i.nonf;
  eq    = i.eq;
  iter  = i.iter;

  return *this;
}

inheritance_vector& inheritance_vector::operator() (const meiosis_map& m)
{
  found = nonf = eq = 0;

  mm = m;

  return *this;
}

inheritance_vector& inheritance_vector::operator++ ()
{
  // If we are iterating over the equivalence class
  if(iter == EQ)
  {
    nonf ^= inc_founder();                // Update the founder and mask nonf
  }
  else
  {
    if(++nonf & ~mm.nonfounder_mask() || nonf == 0)
    {
      storage_type new_mask = inc_founder();   // Increment the founders

      nonf = 0;

      // Set the equiv. class
      eq ^= ((((storage_type) 1) << mm.nonfounder_meiosis_count()) - 1) ^ new_mask;
    }
    else
    {
      // Determine the number of bits that changed in the non-founder.
      unsigned int i = 0;
      for( ; i < mm.nonfounder_meiosis_count() && !(nonf & (1 << i)); ++i);

      eq ^= (((storage_type) 1) << i+1) - 1;
    }
  }
  return *this;
}

inheritance_vector& inheritance_vector::operator-- ()
{
  // If we are iterating over the Equivalence class
  if(iter == EQ)
  {
    nonf ^= dec_founder();
  }
  else
  {
    if(nonf == 0)                            // Can't just decrement nonfounder
    {
      storage_type new_mask = dec_founder();  // decrement the founder bits

      nonf = (((storage_type) 1) << mm.nonfounder_meiosis_count()) - 1;

      eq  ^= nonf ^ new_mask;
    }
    else
    {
      --nonf;                            // Otherwise, decrement nonfounder

      // Determine number of bits that changed in non founder
      unsigned int i = 0;
      for( ; i < mm.nonfounder_meiosis_count() && (nonf & (1 << i)); ++i) ;

      eq ^= ((storage_type) 1 << i+1) - 1;            // Adjust the equiv. class
    }
  }
  return *this;
}

// inc_founder and dec_founder each iterate through the founder bits that
// change, generating a new mask for the equivalence class or non-founder
// bits as they go.

inheritance_vector::storage_type inheritance_vector::inc_founder()
{
  unsigned int i = 0;
  storage_type new_mask = 0;
  for( ; i < mm.founder_count() && (found & (1 << i)); ++i) 
    new_mask |= mm.mask(i);

  if(i < mm.founder_count())
    new_mask |= mm.mask(i);

  ++found;

  return new_mask;
}

inheritance_vector::storage_type inheritance_vector::dec_founder()
{
  unsigned int i = 0;
  storage_type new_mask = 0;
  for( ; i < mm.founder_count() && !(found & (1 << i)); ++i)
  {
    new_mask |= mm.mask(i);
  }

  if(i < mm.founder_count())
    new_mask |= mm.mask(i);

  --found;

  return new_mask;
}

} // End namespace SAGE
