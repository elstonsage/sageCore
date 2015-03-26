#ifndef FIXED_BIT_CALCULATOR_H
#include "lvec/fixed_bit_calculator.h"
#endif

namespace SAGE
{
  
// ===================
// fixed_bit_container
// ===================

inline
fixed_bit_container::fixed_bit_container(const FPED::Subpedigree& p)
   : my_individual_data(p.member_count()),
     my_family_data(p.family_count()),
     my_ind_sib_indices(p.member_count()),
     my_ped(&p)
{
  for(size_t i = 0; i < my_ped->family_count(); ++i)
  {
    const FPED::Family& fam = my_ped->family_index(i);

    size_t sib_count = fam.offspring_count();
    
    my_family_data[i].my_mother_sync_bits.resize(sib_count * sib_count, false);
    my_family_data[i].my_father_sync_bits.resize(sib_count * sib_count, false);

    FPED::OffspringConstIterator o = fam.offspring_begin();

    size_t fs = 0;

    for( ; o != fam.offspring_end(); ++o, ++fs )
    {
      my_ind_sib_indices[o->subindex()] = fs;
    }
  }
}

inline
fixed_bit_container::fixed_bit_container(const self& p)
  : my_individual_data(p.my_individual_data),
    my_family_data(p.my_family_data),
    my_ind_sib_indices(p.my_ind_sib_indices),
    my_ped(p.my_ped)
{ }

inline
fixed_bit_container& fixed_bit_container::operator=(const self& p)
{
  if(this != &p)
  {
    my_individual_data = p.my_individual_data;
    my_family_data = p.my_family_data;
    my_ind_sib_indices = p.my_ind_sib_indices;
  
    my_ped = p.my_ped;
  }

  return *this;
}

inline
fixed_bit_container::~fixed_bit_container()
{ }

inline
const FPED::Subpedigree& fixed_bit_container::get_pedigree() const
{
  return *my_ped;
}

inline
bool fixed_bit_container::mother_fixed(const FPED::Member& i) const
{
  return my_individual_data[i.subindex()].mother_bit_type == fixed;
}

inline
bool fixed_bit_container::father_fixed(const FPED::Member& i) const
{
  return my_individual_data[i.subindex()].father_bit_type == fixed;
}

inline
bool fixed_bit_container::mother_dont_care(const FPED::Member& i) const
{
  return my_individual_data[i.subindex()].mother_bit_type == dont_care;
}

inline
bool fixed_bit_container::father_dont_care(const FPED::Member& i) const
{
  return my_individual_data[i.subindex()].father_bit_type == dont_care;
}


inline
bool fixed_bit_container::mother_synchronized(const FPED::Member& i, const FPED::Member& j) const
{
  FPED::FamilyConstPointer fam = i.family(); 
  
  assert(fam && fam == j.family());
  
  size_t index_i = my_ind_sib_indices[i.subindex()];
  size_t index_j = my_ind_sib_indices[j.subindex()];
  
  return my_family_data[fam->subindex()].my_mother_sync_bits
    [index_i * fam->offspring_count() + index_j];
}

inline
bool fixed_bit_container::father_synchronized(const FPED::Member& i, const FPED::Member& j) const
{
  FPED::FamilyConstPointer fam = i.family(); 
  
  assert(fam && fam == j.family());
  
  size_t index_i = my_ind_sib_indices[i.subindex()];
  size_t index_j = my_ind_sib_indices[j.subindex()];
  
  return my_family_data[fam->subindex()].my_father_sync_bits
    [index_i * fam->offspring_count() + index_j];
}

inline
void fixed_bit_container::set_parent_bit_type(bool ptype, const FPED::Member& i, bit_type t)
{
  if(ptype)
    my_individual_data[i.subindex()].father_bit_type = t;
  else
    my_individual_data[i.subindex()].mother_bit_type = t;
}

inline
void fixed_bit_container::synchronize_parent(bool ptype, const FPED::Member& i, const FPED::Member& j)
{
  FPED::FamilyConstPointer fam = i.family(); 
  
  assert(fam && fam == j.family());
  
  size_t index_i = my_ind_sib_indices[i.subindex()];
  size_t index_j = my_ind_sib_indices[j.subindex()];
  
  if(ptype)
  {
    my_family_data[fam->subindex()].my_father_sync_bits
      [index_i * fam->offspring_count() + index_j] = true;
    my_family_data[fam->subindex()].my_father_sync_bits
      [index_j * fam->offspring_count() + index_i] = true;
  }
  else
  {
    my_family_data[fam->subindex()].my_mother_sync_bits
      [index_i * fam->offspring_count() + index_j] = true;
    my_family_data[fam->subindex()].my_mother_sync_bits
      [index_j * fam->offspring_count() + index_i] = true;
  }
}

inline
fixed_bit_container::individual_data::individual_data()
  : mother_bit_type(unknown),
    father_bit_type(unknown)
{ }

inline
fixed_bit_container::family_data::family_data()
  : my_mother_sync_bits(0),
    my_father_sync_bits(0)
{ }

// ====================
// fixed_bit_calculator
// ====================

inline
fixed_bit_calculator::~fixed_bit_calculator()
{ }
  
inline
const fixed_bit_container& fixed_bit_calculator::get_fixed_bit_container() const
{
  return my_data;
}

}

