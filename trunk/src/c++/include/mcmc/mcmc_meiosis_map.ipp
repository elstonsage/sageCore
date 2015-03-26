// =================================
// MCMC Meiosis Map inline functions
// =================================

/// Returns the multipedigree
///
inline
const FPED::Multipedigree::multipedigree_type&
McmcMeiosisMap::get_multipedigree() const
{
  return *my_mped;
}

/// Returns the subpedigree
///
inline
const FPED::Subpedigree&
McmcMeiosisMap::get_subpedigree() const
{
  return *my_subped;
}

/// Returns the meiosis from the individual's mother
///
inline
size_t
McmcMeiosisMap::get_mother_meiosis(const FPED::Member& id) const
{
  return my_parent_meioses[id.subindex()].moth;
}

/// Returns the meiosis from the individual's father
///
inline
size_t
McmcMeiosisMap::get_father_meiosis(const FPED::Member& id) const
{
  return my_parent_meioses[id.subindex()].fath;
}

/// Returns the number of meioses
/// 
inline
size_t
McmcMeiosisMap::get_meiosis_count() const
{
  return my_meiosis_count;
}

/// Returns the mother of the individual if available, NULL otherwise.
///
inline
FPED::MemberConstPointer
McmcMeiosisMap::get_mother(const FPED::Member& id) const
{
  if( id.is_founder() )
    return NULL;

  if(    id.parent1()->is_female() )
    return id.parent1();

  return id.parent2();
}

/// Returns the father of the individual if available, NULL otherwise.
///
inline
FPED::MemberConstPointer
McmcMeiosisMap::get_father(const FPED::Member& id) const
{
  if( id.is_founder() )
    return NULL;

  if(    !id.parent1()->is_female() )
    return id.parent1();

  return id.parent2();
}

/// Returns the number of founders
///
inline
size_t
McmcMeiosisMap::get_founder_count() const
{
  return my_founder_count;
}

/// Returns the number of nonfounders
///
inline
size_t
McmcMeiosisMap::get_nonfounder_count() const
{
  return my_nonfounder_count;
}

/// Returns the member which is the child of the meiosis given.
///
inline
const FPED::Member&
McmcMeiosisMap::get_child_of_meiosis(size_t i) const
{
  return my_subped->member_index(my_individual_positions[i / 2]);
}

/// Returns the number of individuals in the subpedigree
inline
size_t
McmcMeiosisMap::get_individual_count()     const
{
  return get_subpedigree().member_count();
}

/// Returns the number of nuclear families in the subpedigree
inline
size_t
McmcMeiosisMap::get_family_count() const
{
  return get_subpedigree().family_count();
}

/// Returns \c true if x-linked, \c false otherwise
inline
bool
McmcMeiosisMap::is_x_linked() const
{
  return my_x_linked;
}

