//****************************************************************************
//* File:      reltype.ipp                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines relationtype class to define the     *
//*            relationship between two individuals of a pedigree using      *
//*            numbers.                                                      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

// ---------------------------------------------------------------------------
// Inline implementation of RelationType
// ---------------------------------------------------------------------------

inline
RelationType::RelationType()
{
    my_distance1     = 0;
    my_distance2     = 0;
    my_bridge        = NO_BRIDGE;
    my_phase1        = UNKNOWN_PPHASE;
    my_phase2        = UNKNOWN_PPHASE;
    my_has_offspring = false;
}

inline
RelationType::distance_type
RelationType::distance1() const
{
    return my_distance1;
}

inline
RelationType::distance_type
RelationType::distance2() const
{
    return my_distance2;
}

inline
bridge_type
RelationType::bridge() const
{
    return my_bridge;
}

inline
parental_phase_type
RelationType::phase1() const
{
    return my_phase1;
}

inline
parental_phase_type
RelationType::phase2() const
{
    return my_phase2;
}

inline
bool
RelationType::has_offspring() const
{
    return my_has_offspring;
}

inline
void
RelationType::set_distance1(distance_type d)
{
    my_distance1 = d;
}

inline
void
RelationType::set_distance2(distance_type d)
{
    my_distance2 = d;
}

inline
void
RelationType::set_bridge(bridge_type b)
{
    my_bridge = b;
}

inline
void
RelationType::set_phase1(parental_phase_type p)
{
    my_phase1 = p;
}

inline
void
RelationType::set_phase2(parental_phase_type p)
{
    my_phase2 = p;
}

inline
void
RelationType::set_has_offspring(bool o)
{
    my_has_offspring = o;
}

inline
void 
RelationType::swap_phase()
{
  std::swap(my_phase1,    my_phase2);
  std::swap(my_distance1, my_distance2);
}

inline
bool
RelationType::operator==(const RelationType& r) const
{
    // Handle the quick tests for equality first
    if( distance1() + distance2() != r.distance1() + r.distance2() )
      return false;

    if( bridge() != r.bridge() )
      return false;

    const distance_type   d1 = std::max(  distance1(),   distance2());
    const distance_type   d2 = std::min(  distance1(),   distance2());
    const distance_type r_d1 = std::max(r.distance1(), r.distance2());
    const distance_type r_d2 = std::min(r.distance1(), r.distance2());

    if( d1 != r_d1 || d2 != r_d2 || bridge() != r.bridge() )
      return false;

    return true;
}

inline
bool
RelationType::operator<(const RelationType& r) const
{
    size_t td1 =   distance1() +   distance2();
    size_t td2 = r.distance1() + r.distance2();

    if( !td1 &&   my_bridge == NO_BRIDGE )
      td1 = (distance_type)-1;

    if( !td2 && r.my_bridge == NO_BRIDGE )
      td2 = (distance_type)-1;
    
    if( td1 < td2 )
      return true;

    if( td1 > td2 )
      return false;

    const distance_type   d1 = std::max(  distance1(),   distance2());
    const distance_type   d2 = std::min(  distance1(),   distance2());
    const distance_type r_d1 = std::max(r.distance1(), r.distance2());
    const distance_type r_d2 = std::min(r.distance1(), r.distance2());

    if( d1 < r_d1 )
      return true;

    if( d1 > r_d1 )
      return false;

    if( d2 < r_d2 )
      return true;

    if( d2 > r_d2 )
      return false;

    if( my_bridge < r.my_bridge )
        return true;

    return false;
}

inline
bool
RelationType::operator!=(const RelationType& r) const
{
    return !(*this == r);
}

inline
bool
RelationType::operator>(const RelationType& r) const
{
  return r < *this;
}

inline
bool
RelationType::operator<=(const RelationType& r) const
{
  return !(*this > r);
}

inline
bool
RelationType::operator>=(const RelationType& r) const
{
  return !(*this < r);
}

//-------------------------------------------------------------------

inline
CompleteRelationType::CompleteRelationType() { }

inline
CompleteRelationType::CompleteRelationType(const RelationType &r, const pedigree_member_list& l1,
                                                                  const pedigree_member_list& l2)
        : relationship(r), lineage1(l1), lineage2(l2) { }

inline
void
CompleteRelationType::swap_phase()
{
  relationship.swap_phase();
  lineage1.swap(lineage2);
}

template <class CMP>
void
CompleteRelationType::normalize(const_pedigree_member_pair &p, CMP cmp)
{
  if( relationship.distance1() < relationship.distance2() )
  {
    swap_phase();
    std::swap(p.first, p.second);
  }
  else if( relationship.distance1() == relationship.distance2()
        && !std::lexicographical_compare( lineage1.begin(), lineage1.end() ,
                                          lineage2.begin(), lineage2.end() ,
                                          cmp ))
  {
    swap_phase();
    std::swap(p.first, p.second);
  }
}


// end of RelationType Implementation
