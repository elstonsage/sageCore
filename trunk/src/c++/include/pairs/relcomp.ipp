// ---------------------------------------------------------------------------
// Inline implementation of Relationship Comparison operations
// ---------------------------------------------------------------------------

inline
gender_equivalence::gender_equivalence(bool s)
{
  my_strict = s;
}

inline
bool 
gender_equivalence::operator()(const pedigree_member_type* p1,
                               const pedigree_member_type* p2) const
{
  if( !p1 || !p2 )
    return false;

  MPED::SexCode g1 = p1->get_detailed_sex();
  MPED::SexCode g2 = p2->get_detailed_sex();

  // NOTE: Unknown is eq. to unknown. 

  if(g1 == g2)
    return true;

  if(strict() && g1 != g2)
    return false;

  bool g1m = MPED::is_male(g1);
  bool g2m = MPED::is_male(g2);
  bool g1f = MPED::is_female(g1);
  bool g2f = MPED::is_female(g2);

  if(g1m && g1m == g2m)
    return true;

  if(g2f && g1f == g2f)
    return true;

  return false;
}

inline
bool
gender_equivalence::strict() const
{ 
  return my_strict; 
}

inline
void 
gender_equivalence::set_strict(bool s)
{ 
  my_strict = s;    
}

inline
gender_less::gender_less(bool s)
{
  my_strict = s;
}

inline
bool 
gender_less::operator()(const pedigree_member_type* p1,
                        const pedigree_member_type* p2) const
{
  if( !p1 || !p2 )
    return false;

  MPED::SexCode g1 = p1->get_detailed_sex();
  MPED::SexCode g2 = p2->get_detailed_sex();

  // NOTE: Unknown is eq. to unknown. 

  if( strict() )
    return g1 < g2;

  g1 = get_effective_sex(g1);
  g2 = get_effective_sex(g2);
  
  return g1 < g2;
}

inline
bool
gender_less::strict() const
{ 
  return my_strict; 
}

inline
void 
gender_less::set_strict(bool s)
{ 
  my_strict = s;    
}

//-----------------------------------------------------------------------

inline
RelationMateLess::RelationMateLess()
{ }

inline
bool
RelationMateLess::operator()(const RelationType& r1, const RelationType& r2 ) const
{
  if(  r1 == RelationType() && r2 == RelationType()
   &&  r1.has_offspring() && !r2.has_offspring() )
    return true;

  if( r1 < r2 )
    return true;

  return false;
}

inline
bool
RelationMateLess::operator()(const CompleteRelationType& r1, const CompleteRelationType& r2 ) const
{
  return (*this)(r1.relationship, r2.relationship);
}


//-----------------------------------------------------------------------

template <class MemberCompare>
CompleteRelationLess<MemberCompare>::CompleteRelationLess() 
{ }

template <class MemberCompare>
CompleteRelationLess<MemberCompare>::CompleteRelationLess(const MemberCompare& comp) 
       : member_compare(comp) 
{ }

template <class MemberCompare>
bool
CompleteRelationLess<MemberCompare>::operator()(const CompleteRelationType& r1,
                                                const CompleteRelationType& r2 ) const
{
  if( r1.relationship < r2.relationship )
    return true;

  if( r1.relationship > r2.relationship )
    return false;

  pedigree_member_list const* r1_l1 = &r1.lineage1;
  pedigree_member_list const* r1_l2 = &r1.lineage2;
  pedigree_member_list const* r2_l1 = &r2.lineage1;
  pedigree_member_list const* r2_l2 = &r2.lineage2;


  if( r1_l1->size() < r1_l2->size() )
    swap(r1_l1,r1_l2);
  else if( r1_l1->size() == r1_l2->size() 
        && std::lexicographical_compare(r1_l1->begin(), r1_l1->end(),
                                        r1_l2->begin(), r1_l2->end(),
                                        member_compare ) )
    swap(r1_l1,r1_l2);

  if( r2_l1->size() < r2_l2->size() )
    swap(r2_l1,r2_l2);
  else if( r2_l1->size() == r2_l2->size() 
        && std::lexicographical_compare(r2_l1->begin(), r2_l1->end(),
                                        r2_l2->begin(), r2_l2->end(),
                                        member_compare ) )
    swap(r2_l1,r2_l2);


  if( std::lexicographical_compare(r1_l1->begin(), r1_l1->end(),
                                   r2_l1->begin(), r2_l1->end(),
                                   member_compare ) )
    return true;

  if( std::lexicographical_compare(r2_l1->begin(), r2_l1->end(),
                                   r1_l1->begin(), r1_l1->end(),
                                   member_compare ) )
    return false;

  if( std::lexicographical_compare(r1_l2->begin(), r1_l2->end(),
                                   r2_l2->begin(), r2_l2->end(),
                                   member_compare ) )
    return true;

  return false;
}

//-----------------------------------------------------------------------

template <class MemberCompare>
CompleteRelationMateLess<MemberCompare>::CompleteRelationMateLess()
{ }

template <class MemberCompare>
bool
CompleteRelationMateLess<MemberCompare>::operator()(const CompleteRelationType& r1, const CompleteRelationType& r2 ) const
{
  if( mate_less(r1,r2) )
    return true;

  if( r1.relationship == RelationType() && r2.relationship == RelationType() 
   && r1.relationship.has_offspring()   != r2.relationship.has_offspring() )
    return false;

  if( base_relation_less(r1,r2) )
    return true;

  return false;
}

