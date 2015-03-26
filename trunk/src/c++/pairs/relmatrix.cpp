//****************************************************************************
//* File:      relmatrix.cpp                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song & Kevin Jacobs                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     This source file defines a matrix of relationships            *
//*            individuals in a pedigree                                     *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <list>
#include <stdio.h>
#include "rped/rped.h"
#include "pairs/reldefs.h"
#include "pairs/reltype.h"
#include "pairs/reltypename.h"
#include "pairs/relmatrix.h"
#include "containers/ldatamap.h"

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     RelationMatrix                                               ~
// ~                                                                         ~
// ~ Purpose:   Generate & hold relationship in pedigrees.                   ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


// ---------------------------------------------------------------------------
// Out-of-line implementation of RelationMatrix
// ---------------------------------------------------------------------------

RelationMatrix::RelationMatrix(const RPED::RefMultiPedigree* mp)
{
  my_multiped = mp;
  generate_matrix();
}

void
RelationMatrix::generate_matrix()
{
  // Local declaration to simplify the code.
  const RPED::RefMultiPedigree* mp = my_multiped;

  size_t pedigree_count = mp->pedigree_count();
  my_relmatrix.resize(0);
  my_relmatrix.resize(pedigree_count);

  size_t i1;
  size_t i2;

  for( size_t ped = 0; ped < pedigree_count; ++ped )
  {
    const RPED::RefPedigree* p = &mp->pedigree_index(ped);

    if(!p)
      continue;

    size_t member_count = p->member_count();
    my_relmatrix[ped].resize(member_count);

    for( i1 = 0; i1 < member_count; ++i1 )
      for( i2 = 0; i2 <= i1; ++i2 )
        compute_relationship( &p->member_index(i1), &p->member_index(i2) );
  }
}

void
RelationMatrix::set_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2, RelationType rel)
{
  if(!i1 || !i2 || i1->pedigree() != i2->pedigree())
    return;

  size_t ped = i1->pedigree()->index();
  size_t i   = i1->index();
  size_t j   = i2->index();

  if( i < j )
    rel.swap_phase();

  my_relmatrix[ped](i,j) = rel;
}

void
RelationMatrix::compute_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2)
{
  RelationType rel;
  rel.set_distance1(0);
  rel.set_distance2(0);

  if( !i1 || !i2 || i1->pedigree() != i2->pedigree() )
    return;

  if( i1 == i2 )
  {
    rel.set_bridge(SELF_BRIDGE);
    set_relationship(i1, i2, rel);

    if( i1->parent1() && i1->parent2() )
    {
      rel = get_relationship(i1->parent1(), i1->parent2());
      rel.set_has_offspring(true);
      set_relationship(i1->parent1(), i1->parent2(), rel);
    }

    return;
  }

  // Check for parent-offspring
  if( i1->parent1() == i2 || i1->parent2() == i2 )
  {
    rel.set_distance1(1);

    if( i1->parent1() == i2 )
      rel.set_phase1(PARENT1_PHASE);
    else
      rel.set_phase1(PARENT2_PHASE);

    set_relationship(i1, i2, rel);

    return;
  }

  // Check for full-sibling type.
  if( i1->family() && i1->family() == i2->family() )
  {
    rel.set_distance1(1);
    rel.set_distance2(1);
    rel.set_bridge(FULL_SIB_BRIDGE);
    set_relationship(i1, i2, rel);
    return;
  }

  const pedigree_member_type* i1_p1 = std::min( i1->parent1(), i1->parent2() );
  const pedigree_member_type* i1_p2 = std::max( i1->parent1(), i1->parent2() );
  const pedigree_member_type* i2_p1 = std::min( i2->parent1(), i2->parent2() );
  const pedigree_member_type* i2_p2 = std::max( i2->parent1(), i2->parent2() );

  if( !i1_p2 && !i2_p2 )
  {
    // Set no relation
    return;
  }

  if( i1_p2 && i1_p1 == i2_p1 && i1_p2 == i2_p2 )
  {
    rel.set_distance1(1);
    rel.set_distance2(1);
    rel.set_bridge(FULL_SIB_BRIDGE);
    set_relationship(i1, i2, rel);
    return;
  }

  // Half sibs
  if(   (i1_p1 && i1_p1 == i2_p1) || (i1_p1 && i1_p1 == i2_p2)
     || (i1_p2 && i1_p2 == i2_p1) || (i1_p2 && i1_p2 == i2_p2) )
  {
    rel.set_distance1(1);
    rel.set_distance2(1);
    rel.set_bridge(HALF_SIB_BRIDGE);

    if( i1->parent1() && i1->parent1() == i2->parent1() )
    {
      rel.set_phase1(PARENT1_PHASE);
      rel.set_phase2(PARENT1_PHASE);
    }
    else if( i1->parent2() && i1->parent2() == i2->parent1() )
    {
      rel.set_phase1(PARENT2_PHASE);
      rel.set_phase2(PARENT1_PHASE);
    }
    else if( i1->parent1() && i1->parent1() == i2->parent2() )
    {
      rel.set_phase1(PARENT1_PHASE);
      rel.set_phase2(PARENT2_PHASE);
    }
    else if( i1->parent2() && i1->parent2() == i2->parent2() )
    {
      rel.set_phase1(PARENT2_PHASE);
      rel.set_phase2(PARENT2_PHASE);
    }
    set_relationship(i1, i2, rel);
    return;
  }

  // If none of above, then find relation type between i1's parents and i2.
  rel = infer_from_parental_relationship(i1, i2);
  set_relationship(i1, i2, rel);
}

RelationType
RelationMatrix::infer_from_parental_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2) const
{
  if(!i1 || !i2)
    return RelationType();

  RelationType p_rel1;
  RelationType p_rel2;
  RelationType no_relation;

  size_t i1_p1 = (i1->parent1())? i1->parent1()->index() : 0;
  size_t i1_p2 = (i1->parent2())? i1->parent2()->index() : 0;
  size_t i2_p1 = (i2->parent1())? i2->parent1()->index() : 0;
  size_t i2_p2 = (i2->parent2())? i2->parent2()->index() : 0;

  size_t mi1 = std::max(i1_p1, i1_p2);
  size_t mi2 = std::max(i2_p1, i2_p2);

  if( mi2 == std::max(mi1, mi2) )
    swap(i1, i2);

  p_rel1 = get_relationship(i1->parent1(), i2);
  p_rel2 = get_relationship(i1->parent2(), i2);

  p_rel1.set_has_offspring(false);
  p_rel2.set_has_offspring(false);

  if( p_rel1 != no_relation )
  {
    if( i1->parent1()->index() < i2->index() )
    {
      if(    (p_rel1.distance1() == 2 && p_rel1.distance2() == 0)
          && (p_rel2.distance1() == 0 && p_rel2.distance2() == 0) )
        p_rel1.set_bridge(HALF_SIB_BRIDGE);

      p_rel1.swap_phase();
    }

    p_rel1.set_phase1(PARENT1_PHASE);
    p_rel1.set_distance1( p_rel1.distance1() + 1 );
  }

  if( p_rel2 != no_relation )
  {
    if( i1->parent2()->index() < i2->index() )
    {
      if(    (p_rel2.distance1() == 2 && p_rel2.distance2() == 0)
          && (p_rel1.distance1() == 0 && p_rel1.distance2() == 0) ) 
        p_rel2.set_bridge(HALF_SIB_BRIDGE);

      p_rel2.swap_phase();
    }

    p_rel2.set_phase1(PARENT2_PHASE);
    p_rel2.set_distance1( p_rel2.distance1() + 1 );
  }

  if(    !p_rel1.distance1()
      && !p_rel1.distance2()
      && p_rel1.distance1() ==  p_rel2.distance1()
      && p_rel1.distance2() ==  p_rel2.distance2()
      && p_rel1.bridge()    ==  p_rel2.bridge()    )
  {
    p_rel1.set_phase1(UNKNOWN_PPHASE);
    p_rel2.set_phase1(UNKNOWN_PPHASE);
    p_rel1.set_bridge(NO_BRIDGE);
    p_rel2.set_bridge(NO_BRIDGE);
  }


  RelationType min_rel = std::min(p_rel1, p_rel2);

  if( mi2 == std::max(mi1, mi2) )
    min_rel.swap_phase();

  return min_rel;
}

bool
RelationMatrix::gender_equivalent(const pedigree_member_list& a1, const pedigree_member_list& a2, bool strict) const
{
  if( a1.size() != a2.size() )
    return false;

  return std::equal(a1.begin(), a1.end(), a2.begin(), gender_equivalence(strict) );
}

bool
RelationMatrix::gender_equivalent(const const_pedigree_member_pair& p1, const const_pedigree_member_pair& p2, bool strict) const
{
  // First we must be non-gender equivalent -- if not, we can't be gender equiv.
  if( !equivalent(p1, p2) )
    return false;

  pedigree_member_list i11_ans;
  pedigree_member_list i12_ans;
  pedigree_member_list i21_ans;
  pedigree_member_list i22_ans;

  // We include the person in their ancestor list
  build_common_lineage(p1.first, p1.second, i11_ans, i12_ans);
  build_common_lineage(p2.first, p2.second, i21_ans, i22_ans);

  if( gender_equivalent(i11_ans, i21_ans, strict) && gender_equivalent(i12_ans, i22_ans, strict) )
    return true;

  if( gender_equivalent(i11_ans, i22_ans, strict) && gender_equivalent(i12_ans, i21_ans, strict) )
    return true;

  return false;
}

RelationType
RelationMatrix::build_common_lineage(const pedigree_member_type* i1,
                                     const pedigree_member_type* i2,
                                           pedigree_member_list& a1,
                                           pedigree_member_list& a2) const
{
  RelationType rel = get_relationship(i1, i2);

  if( i1 == i2 )

  {
    a1.push_back( const_cast<pedigree_member_type*>(i1) );
    return rel;
  }

  if( !rel.distance1() && rel.distance1() == rel.distance2() )
  {
    a1.push_back( const_cast<pedigree_member_type*>(i1) );
    a2.push_back( const_cast<pedigree_member_type*>(i2) );
    return rel;
  }

  // Build ancestor list for phase 1
  build_lineage_list(i1, i2, a1, rel);

  rel.swap_phase();

  // Build ancestor list for phase 2
  build_lineage_list(i2, i1, a2, rel);

  rel.swap_phase();
  return rel;
}

void
RelationMatrix::build_lineage_list(const pedigree_member_type* i1,
                                   const pedigree_member_type* i2,
                                         pedigree_member_list& a,
                                         RelationType rel) const
{
  pedigree_member_type* m = const_cast<pedigree_member_type*>(i1);

  a.clear();

  if( !rel.distance1() && rel.phase1() == UNKNOWN_PPHASE )
    return;

  // Keep going so long as we have a ladder
  while( m )
  {
    // Store the current step
    a.push_back(m);

    // Pick the right direction
    if     ( rel.phase1() == PARENT1_PHASE )
      m = m->parent1();
    else if( rel.phase1() == PARENT2_PHASE )
      m = m->parent2();
    else
      break;

    // Climb up to the next rung
    rel = get_relationship(m, i2);

    // Normalize phase (pick a consistent direction)
    if( m->index() < i2->index() )
      rel.swap_phase();
  }
}

void
RelationMatrix::view(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw) const
{
  for( size_t k = 0; k < my_relmatrix.size(); ++k )
  {
    size_t rank = my_relmatrix[k].size();

    out << std::endl << std::left
        << "--------------------------------------------------------------------"
        << std::endl
        << "Pedigree: " << mp.pedigree_index(k).name()
        << std::endl << std::endl;

    for( size_t i = 0; i < rank ; ++i )
    {
      if( i % 8 == 0 )
      {
        // Handle boundry conditions
        size_t stop = i;
        if( !i || (i+8) >= rank)
          stop = rank;

        out << std::setw(11) << "";
        for( size_t j = 0; j < stop ; ++j )
          out << std::setw(11) << mp.pedigree_index(k).member_index(j).name();
        out << std::endl;
      }

      out << "  "
          << std::setw(5) << mp.pedigree_index(k).member_index(i).name()
          << ": |";

      for( size_t j = 0; j <= i; ++j )
        if(raw)
          out << setw(10) << my_relmatrix[k](i,j).str()              << "|";
        else
          out << setw(10) << relationship_name(my_relmatrix[k](i,j)) << "|";
      out << std::endl;
    }
  }
}

void
RelationMatrix::view_pairs(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw) const
{
  typedef   LinearDatamap<RelationType,
                          const_pedigree_member_pair,
                          RelationMateLess>          relation_map;

  relation_map relpairs;

  for( size_t k = 0; k < my_relmatrix.size(); ++k )
  {
    size_t rank = my_relmatrix[k].size();

    for( size_t i = 0; i < rank ; ++i )
      for( size_t j = 0; j <= i; ++j )
      {
        const_pedigree_member_pair p( &mp.pedigree_index(k).member_index(i),
                                      &mp.pedigree_index(k).member_index(j) );

        RelationType c = get_relationship(p);

        relpairs.insert( make_pair(c, p) );
      }
  }

  relation_map::type_const_iterator rel_type = relpairs.type_begin();
  relation_map::type_const_iterator rel_end  = relpairs.type_end();

  for( ; rel_type != rel_end; ++rel_type )
  {
    out << setw(12) << relationship_name(rel_type->first) << " [n="
        << setw(3)  << std::right << rel_type->second.size() << "]: " << std::left;

    relation_map::value_const_iterator        p = rel_type->second.begin();
    relation_map::value_const_iterator pair_end = rel_type->second.end();

    for( size_t k = 0; p != pair_end; ++p, ++k )
    {
      if( k && k % 10 == 0 )
        out << endl << setw(22) << "";
      out << '(' << setw(2) << p->first->name()
          << ',' << setw(2) << p->second->name() << ") ";
    }
    out << endl;
  }
}

void
RelationMatrix::view_pairs_gender(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw) const
{
  typedef   LinearDatamap<CompleteRelationType,
                         const_pedigree_member_pair,
                         CompleteRelationMateLess<> >  relation_map;

  relation_map relpairs;

  for( size_t k = 0; k < my_relmatrix.size(); ++k )
  {
    size_t rank = my_relmatrix[k].size();

    for( size_t i = 0; i < rank ; ++i )
      for( size_t j = 0; j <= i; ++j )
      {
        const_pedigree_member_pair p( &mp.pedigree_index(k).member_index(i),
                                      &mp.pedigree_index(k).member_index(j) );

        CompleteRelationType c = get_complete_relationship(p);

        c.normalize(p, gender_less());

        relpairs.insert(make_pair(c, p));
      }
  }

  relation_map::type_const_iterator rel_type = relpairs.type_begin();
  relation_map::type_const_iterator rel_end  = relpairs.type_end();

  for( ; rel_type != rel_end; ++rel_type )
  {
    out << setw(12)
        << relationship_name(rel_type->first.relationship) << " ";

    std::pair<string, string> gs = gender_strings(rel_type->first);

    string gender_list;
    if( gs.first.size() && gs.second.size() )
      gender_list = gs.first + "," + gs.second;
    else
      gender_list = gs.first + gs.second;

    out << setw(9) << gender_list;

    out << " [n=" << setw(3) << std::right << rel_type->second.size() << "]: " << std::left;

    relation_map::value_const_iterator        p = rel_type->second.begin();
    relation_map::value_const_iterator pair_end = rel_type->second.end();

    for( size_t k = 0; p != pair_end; ++p, ++k )
    {
      if( k && k % 10 == 0 )
        out << endl << setw(32) << "";
      out << '(' << setw(2) << p->first->name()
          << ',' << setw(2) << p->second->name()
          << ") ";
    }
    out << endl;
  }
}

// end of RelationMatrix Implementation

} // end of namespace SAGE
