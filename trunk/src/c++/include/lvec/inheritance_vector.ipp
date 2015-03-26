// =============================
// iv_reference inline functions
// =============================

iv_reference::operator bool() const
{
  return !(!(*p &mask));
}

iv_reference::reference& iv_reference::operator=(bool x)
{
  if(x != bool(*this))
  {
    *eq ^= eqmask; 
    if(x) *p |= mask;
    else  *p &= ~mask;
  }
  return *this;
}

iv_reference::reference& iv_reference::operator=(const reference& x)
{
  return *this = bool(x);
}

bool iv_reference::operator==(const reference& x) const
{
  return bool(*this) == bool(x);
}

bool iv_reference::operator<(const reference& x) const
{
  return bool(*this) < bool(x);
}

void iv_reference::flip()
{
  *p ^= mask; *eq ^= eqmask;
}

//====================================
// Inheritance Vector inline functions
//====================================

inline inheritance_vector::storage_type
    inheritance_vector::get_founders() const
{
  return found;
}

inline inheritance_vector::storage_type inheritance_vector::get_nonfounders() const
{
  return nonf; 
}

inline inheritance_vector::equivalence_class inheritance_vector::get_equivalence_class() const
{
  return eq;
}

inline inheritance_vector::equivalence_class
    inheritance_vector::set_equivalence_class(equivalence_class s) 
{
  nonf = eq = s & mm.nonfounder_mask();

  found = 0;

  return eq;
}

inline inheritance_vector::size_type
inheritance_vector::num_equivalence_classes() const
{
  return 1 << mm.nonfounder_meiosis_count();
}

inline inheritance_vector::reference inheritance_vector::meiosis(index i)
{
  meiosis_map::mask_type msk = mm.mask(i);

  if(i < meiosis_map::meiosis_bits)
    return reference(&found, &eq, 1 << i, msk);
  else
    return reference(&nonf,  &eq, 1 << (i - meiosis_map::meiosis_bits), msk);
}

inline inheritance_vector::const_reference inheritance_vector::meiosis(index i) const
{
  if(i < meiosis_map::meiosis_bits)
    return found & (storage_type) (1 << i);
  else
    return nonf  & (storage_type) (1 << (i - meiosis_map::meiosis_bits));
}

inline inheritance_vector::reference inheritance_vector::operator[] (index i)
{
  return meiosis(i);
}

inline inheritance_vector::const_reference
    inheritance_vector::operator[] (index i) const
{
  return meiosis(i);
}

inline inheritance_vector::reference
    inheritance_vector::mother_meiosis(member_const_pointer i)
{
  return meiosis(mm.mother_meiosis(i));
}

inline inheritance_vector::const_reference
    inheritance_vector::mother_meiosis(member_const_pointer i) const
{
  return meiosis(mm.mother_meiosis(i));
}

inline inheritance_vector::reference
    inheritance_vector::father_meiosis(member_const_pointer i)
{
  return meiosis(mm.father_meiosis(i));
}

inline inheritance_vector::const_reference
    inheritance_vector::father_meiosis(member_const_pointer i) const
{
  return meiosis(mm.father_meiosis(i));
}

inline inheritance_vector::iter_mthd inheritance_vector::get_iteration_method() const
{
  return iter;
}

inline inheritance_vector::iter_mthd inheritance_vector::set_iteration_method(iter_mthd i)
{
  return iter = i;
}

inline inheritance_vector inheritance_vector::operator++ (int)
{
  inheritance_vector v = *this;
  ++*this;
  return v;
}

inline inheritance_vector inheritance_vector::operator-- (int)
{
  inheritance_vector v = *this;
  --*this;
  return v;
}

inline bool inheritance_vector::isBegin() const
{
  if(iter == EQ)
    return !found;
  else
    return !(found || nonf);
}

inline bool inheritance_vector::isEnd() const
{
  return found >= ((storage_type) 1) << mm.founder_count();
}

