#ifndef RELATIONCOMPARE_H
#define RELATIONCOMPARE_H

//****************************************************************************
//* File:      relcomp.h                                                     *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines comparison operations between        *
//*            relationships.                                                *
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

namespace SAGE
{

class gender_equivalence : std::binary_function<const pedigree_member_type*,
                                                const pedigree_member_type*,
                                                bool>
{
public:

  explicit      gender_equivalence(bool strict = false);

  bool                    operator()(const RPED::RefPedigree::member_type*,
                                     const RPED::RefPedigree::member_type*)  const;

  bool                      strict()                    const;
  void                  set_strict(bool s);

private:

  bool     my_strict;
};

class gender_less :
     std::binary_function<const pedigree_member_type*,
                          const pedigree_member_type*,
                          bool>
{
public:

  explicit      gender_less(bool strict = false);

  bool                    operator()(const pedigree_member_type*,
                                     const pedigree_member_type*)  const;

  bool                      strict()                    const;
  void                  set_strict(bool s);

private:

  bool     my_strict;
};

struct RelationMateLess : std::binary_function<const RelationType&,
                                               const RelationType&,
                                               bool>
{
  RelationMateLess();

  bool operator()(const RelationType& r1,
                  const RelationType& r2 ) const;

  bool operator()(const CompleteRelationType& r1,
                  const CompleteRelationType& r2 ) const;
};

template <class MemberCompare = gender_less>
struct CompleteRelationLess : std::binary_function<const CompleteRelationType&,
                                                   const CompleteRelationType&,
                                                   bool>
{
  CompleteRelationLess();
  CompleteRelationLess(const MemberCompare& comp);

  bool operator()(const CompleteRelationType& r1,
                  const CompleteRelationType& r2 ) const;

  MemberCompare   member_compare;
};

template <class MemberCompare = gender_less>
struct CompleteRelationMateLess : std::binary_function<const CompleteRelationType&,
                                                       const CompleteRelationType&,
                                                       bool>
{
  CompleteRelationMateLess();
  CompleteRelationMateLess(const MemberCompare& comp);

  bool operator()(const CompleteRelationType& r1,
                  const CompleteRelationType& r2 ) const;

  CompleteRelationLess<MemberCompare>     base_relation_less;
  RelationMateLess                        mate_less;
};

#include "pairs/relcomp.ipp"

} // end of namespace SAGE

#endif
