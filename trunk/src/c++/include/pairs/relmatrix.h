#ifndef RELATIONMATRIX_H
#define RELATIONMATRIX_H

//****************************************************************************
//* File:      relmatrix.h                                                   *
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

#include <algorithm>
#include <cmath>
#include <list>
#include "rped/rped.h"
#include "pairs/reldefs.h"
#include "pairs/reltype.h"
#include "pairs/relcomp.h"
#include "numerics/trimatrix.h"

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     RelationMatrix                                               ~
// ~                                                                         ~
// ~ Purpose:   Generate & hold relationship in pedigrees.                   ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class RelationMatrix
{
  public:

    RelationMatrix(const RPED::RefMultiPedigree* mp);

    // Test for relationship equivalance between two pairs of relatives.
    //   The first version ignores the gender of all individuals.

    bool                        equivalent(const pedigree_member_pair& p1, const pedigree_member_pair& p2) const;
    bool                        equivalent(const const_pedigree_member_pair& p1, const const_pedigree_member_pair& p2) const;

    //   The second version if conditional on all the genders of the
    //   individuals as well as common ancestors (given by
    //   build_common_lineage).  If the strict flag is set to true, then
    //   assumed genders are ignored.  Unknown gender is always treated
    //   strictly in the sense that only it only matches other unknown
    //   genders.

    bool                   gender_equivalent(const pedigree_member_pair& p1, const pedigree_member_pair& p2, bool strict) const;
    bool                   gender_equivalent(const const_pedigree_member_pair& p1, const const_pedigree_member_pair& p2, bool strict) const;

    RelationType           get_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2)  const;
    RelationType           get_relationship(const pedigree_member_pair& p)  const;
    RelationType           get_relationship(const const_pedigree_member_pair& p)  const;

    CompleteRelationType   get_complete_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2)  const;
    CompleteRelationType   get_complete_relationship(const pedigree_member_pair& p)  const;
    CompleteRelationType   get_complete_relationship(const const_pedigree_member_pair& p)  const;

    // Build the list of each person and their closest common ancestors.
    // For unilineal pairs with full-sib bridge, this does not include the
    // two ancestoral parents.  For bilineal pairs (other than simple full
    // sibs) these lists approximate the shortest path to a common ancestor.


    void                       view(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw = false) const;
    void                 view_pairs(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw = false) const;
    void          view_pairs_gender(std::ostream& out, const RPED::RefMultiPedigree& mp, bool raw = false) const;

  private:

    const RPED::RefMultiPedigree*                               my_multiped;
    vector< TriangleMatrix<RelationType> >                my_relmatrix;


    void                generate_matrix();
    void                compute_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2);
    void                    set_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2, RelationType r);

    RelationType        infer_from_parental_relationship(const pedigree_member_type* i1, const pedigree_member_type* i2) const;

    bool                gender_equivalent(const pedigree_member_list& a1, const pedigree_member_list& a2, bool strict) const;

    RelationType        build_common_lineage(const pedigree_member_type* i1,
                                             const pedigree_member_type* i2,
                                             pedigree_member_list& a1,
                                             pedigree_member_list& a2)    const;

    void                build_lineage_list(const pedigree_member_type* i1,
                                           const pedigree_member_type* i2,
                                           pedigree_member_list& a,
                                           RelationType rel)              const;

}; // end of class definition

#include "pairs/relmatrix.ipp"

} // end of namespace SAGE

#endif
