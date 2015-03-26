// ============================
// Meiosis Map inline functions
// ============================

inline
const FPED::Multipedigree&
meiosis_map::get_multipedigree() const
{
  return *my_filtered_mped;
}

inline
meiosis_map::pedigree_const_pointer
meiosis_map::get_pedigree() const
{
  return my_filtered_subped->pedigree();
}

inline
meiosis_map::subpedigree_const_pointer
meiosis_map::get_subpedigree() const
{
  return my_filtered_subped;
}

inline
bool
meiosis_map::built() const
{
  return my_built;
}

inline
meiosis_map::index
meiosis_map::mother_meiosis(member_const_pointer id) const
{
  if( id )
    return p_meioses[id->subindex()].moth;

  return index_err;
}

inline
meiosis_map::index
meiosis_map::father_meiosis(member_const_pointer id) const
{
  if( id )
    return p_meioses[id->subindex()].fath;

  return index_err;
}

inline
meiosis_map::index
meiosis_map::mother_meiosis(size_type i) const
{
  return p_meioses[i].moth;
}

inline
meiosis_map::index
meiosis_map::father_meiosis(size_type i) const
{
  return p_meioses[i].fath;
}

inline
const meiosis_map::mask_type&
meiosis_map::founder_mask() const
{
  return fm;
}

inline
const meiosis_map::mask_type&
meiosis_map::nonfounder_mask() const
{
  return nfm;
}

inline
meiosis_map::mask_type
meiosis_map::mask(index i) const
{
  if(i < masks.size()) return masks[i];
 
  if(meiosis_bits <= i && i < 2 * meiosis_bits)
    return ((storage_type) 1) << (i - meiosis_bits);

  return 0;
}

inline
meiosis_map::mask_type
meiosis_map::mother_mask(size_type i) const
{
  return mask(mother_meiosis(i));
}

inline
meiosis_map::mask_type
meiosis_map::father_mask(size_type i) const
{
  return mask(father_meiosis(i));
}

inline
meiosis_map::mask_type
meiosis_map::mother_mask(member_const_pointer i) const
{
  return mask(mother_meiosis(i));
}

inline
meiosis_map::mask_type
meiosis_map::father_mask(member_const_pointer i) const
{
  return mask(father_meiosis(i));
}

inline
meiosis_map::size_type
meiosis_map::meiosis_count() const
{
  return founder_meiosis_count() + nonfounder_meiosis_count();
}

inline
meiosis_map::size_type
meiosis_map::founder_meiosis_count() const
{
  return fbits;
}

inline
meiosis_map::size_type
meiosis_map::nonfounder_meiosis_count() const
{
  return nfbits;
}

inline
meiosis_map::member_const_pointer
meiosis_map::member(size_type i) const
{
  const member_type& mem = my_filtered_subped->member_index(i);

  return &mem;
}

inline
meiosis_map::member_const_pointer
meiosis_map::operator[](size_type i) const
{
  const member_type& mem = my_filtered_subped->member_index(i);

  return &mem;
}

inline
bool
meiosis_map::founder(member_const_pointer id) const
{
  const member_type& mem = my_filtered_subped->member_index(id->subindex());

  if( mem.parent1() || mem.parent2() )
    return false;

  return true;
}

inline
bool
meiosis_map::founder(size_type i) const
{
  const member_type& mem = my_filtered_subped->member_index(i);

  if( mem.parent1() || mem.parent2() )
    return false;

  return true;
}

inline
meiosis_map::member_const_pointer
meiosis_map::mother(member_const_pointer id) const
{
  if( founder(id) )
    return NULL;

  const member_type& mem = my_filtered_subped->member_index(id->subindex());

  if(mem.parent1()->is_female() && !mem.parent2()->is_female())
    return mem.parent1();

  else if(mem.parent2()->is_female() && !mem.parent1()->is_female())
    return mem.parent2();

  return max(mem.parent1(), mem.parent2());
}

inline
meiosis_map::member_const_pointer
meiosis_map::father(member_const_pointer id) const
{
  if( founder(id) )
    return NULL;

  const member_type& mem = my_filtered_subped->member_index(id->subindex());

  if(mem.parent1()->is_male() && !mem.parent2()->is_male())
    return mem.parent1();

  else if(mem.parent2()->is_male() && !mem.parent1()->is_male())
    return mem.parent2();

  return min(mem.parent1(), mem.parent2());
}

inline
meiosis_map::member_const_pointer
meiosis_map::mother(size_type i) const
{
  if( founder(i) )
    return NULL;

  const member_type& mem = my_filtered_subped->member_index(i);

  if(mem.parent1()->is_female() && !mem.parent2()->is_female())
    return mem.parent1();

  else if(mem.parent2()->is_female() && !mem.parent1()->is_female())
    return mem.parent2();

  return max(mem.parent1(), mem.parent2());
}

inline
meiosis_map::member_const_pointer
meiosis_map::father(size_type i) const
{
  if( founder(i) )
    return NULL;

  const member_type& mem = my_filtered_subped->member_index(i);

  if(mem.parent1()->is_male() && !mem.parent2()->is_male())
    return mem.parent1();

  else if(mem.parent2()->is_male() && !mem.parent1()->is_male())
    return mem.parent2();

  return min(mem.parent1(), mem.parent2());
}

inline
meiosis_map::size_type
meiosis_map::mother_index(member_const_pointer id) const
{
  if( mother(id) == NULL )
    return index_err;

  return mother(id)->subindex();
}

inline
meiosis_map::size_type
meiosis_map::father_index(member_const_pointer id) const
{
  if( father(id) == NULL )
    return index_err;

  return father(id)->subindex();
}
                                                      
inline
meiosis_map::size_type
meiosis_map::mother_index(size_type i) const
{
  if( mother(i) == NULL )
    return index_err;

  return mother(i)->subindex();
}

inline
meiosis_map::size_type
meiosis_map::father_index(size_type i) const
{
  if( father(i) == NULL )
    return index_err;

  return father(i)->subindex();
}

inline
meiosis_map::size_type
meiosis_map::founder_count() const
{
  return my_founder_count;
}

inline
meiosis_map::size_type
meiosis_map::nonfounder_count() const
{
  return my_nonfounder_count;
}

inline
size_t
meiosis_map::family_count() const
{
  return get_subpedigree()->family_count();
}

inline
size_t
meiosis_map::individual_count() const
{
  return get_subpedigree()->member_count();
}

inline
meiosis_map::size_type
meiosis_map::bit_count() const
{
  return 2 * nonfounder_count() - founder_count();
}

inline
bool
meiosis_map::is_x_linked() const
{
  return my_x_linked;
}

inline
bool
meiosis_map::is_father_bit(index p) const
{
  for(size_t i = 0; i < p_meioses.size(); ++i )
  {
    if( mask(p_meioses[i].fath) == p )
      return true;
  }

  return false;
}
