#ifndef RELATIONDEFS_H
#define RELATIONDEFS_H

//****************************************************************************
//* File:      reldefs.h                                                     *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     This header file defines various types used in computing      *
//*            relationships between individuals in a reference pedigree.    *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <utility>
#include <list>
#include "rped/rped.h"

namespace SAGE
{

typedef RPED::RefPedigree::member_type                        pedigree_member_type;
typedef RPED::RefPedigree::member_pointer                     pedigree_member_pointer;

typedef std::pair<pedigree_member_pointer,
                  pedigree_member_pointer>              pedigree_member_pair;

typedef std::pair<const pedigree_member_type*,
                  const pedigree_member_type*>          const_pedigree_member_pair;

typedef std::list<pedigree_member_pointer>              pedigree_member_list;

typedef std::pair<pedigree_member_list,
                  pedigree_member_list>                 pedigree_member_list_pair;

} // end of namespace SAGE

#endif
