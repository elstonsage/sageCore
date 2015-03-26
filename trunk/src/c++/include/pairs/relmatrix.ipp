// ---------------------------------------------------------------------------
// Inline implementation of RelationMatrix
// ---------------------------------------------------------------------------

inline
RelationType
RelationMatrix::get_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2) const
{
  if(!i1 || !i2 || i1->pedigree() != i2->pedigree())
    return RelationType();

  size_t ped = i1->pedigree()->index();
  size_t i   = i1->index();
  size_t j   = i2->index();

  return my_relmatrix[ped](i,j);
}

inline
RelationType
RelationMatrix::get_relationship(const pedigree_member_pair& p) const
{
  return get_relationship(p.first, p.second);
}

inline
RelationType
RelationMatrix::get_relationship(const const_pedigree_member_pair& p) const
{
  return get_relationship(p.first, p.second);
}

inline
CompleteRelationType
RelationMatrix::get_complete_relationship(const pedigree_member_pair& p) const
{
  pedigree_member_list a1;
  pedigree_member_list a2;
  
  RelationType rel = build_common_lineage(p.first, p.second, a1, a2);

  return CompleteRelationType(rel, a1, a2);
}

inline
CompleteRelationType
RelationMatrix::get_complete_relationship(const const_pedigree_member_pair& p) const
{
  pedigree_member_list a1;
  pedigree_member_list a2;
  
  RelationType rel = build_common_lineage(p.first, p.second, a1, a2);

  return CompleteRelationType(rel, a1, a2);
}

inline
CompleteRelationType
RelationMatrix::get_complete_relationship(const pedigree_member_type* i1,
                                          const pedigree_member_type* i2) const
{
  pedigree_member_list a1;
  pedigree_member_list a2;
  
  RelationType rel = build_common_lineage(i1, i2, a1, a2);

  return CompleteRelationType(rel, a1, a2);
}

inline
bool
RelationMatrix::equivalent(const pedigree_member_pair& p1, const pedigree_member_pair& p2) const
{
  return (get_relationship(p1.first, p1.second) == get_relationship(p2.first, p2.second));
}

inline
bool
RelationMatrix::equivalent(const const_pedigree_member_pair& p1, const const_pedigree_member_pair& p2) const
{
  return (get_relationship(p1.first, p1.second) == get_relationship(p2.first, p2.second));
}

inline
bool
RelationMatrix::gender_equivalent(const pedigree_member_pair& p1, const pedigree_member_pair& p2, bool strict) const
{
  return gender_equivalent( *reinterpret_cast<const const_pedigree_member_pair*>(&p1),
                            *reinterpret_cast<const const_pedigree_member_pair*>(&p2),
                             strict );
}

// end of RelationMatrix Implementation
