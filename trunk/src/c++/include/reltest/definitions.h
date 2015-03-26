#ifndef RELTEST_DEFS_H
#define RELTEST_DEFS_H

//==========================================================================
//  File:      definitions.h
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                              Jul. 03
//
//  Notes:     This header file contains various type definition statements
//             used through out the reltest files.
//
//  Copyright (c) 2003 R.C. Elston
//    All Rights Reserved
//==========================================================================

#include "numerics/print_util.h"
#include "numerics/histogram.h"
#include "util/dots.h"
#include "rped/genome_description.h"
#include "ibd/exact_ibd_analysis.h"

namespace SAGE
{

namespace RELTEST
{

typedef RPED::RefMultiPedigree::subpedigree_type             subpedigree;
typedef RPED::RefMultiPedigree::subpedigree_const_pointer    subpedigree_const_pointer;
typedef const RPED::RefMultiPedigree::member_type            ind_type;
typedef RPED::RefMultiPedigree::member_const_pointer         ind_id;
typedef std::pair<ind_id, ind_id>                      id_pair;

typedef RPED::genome_description::region_type                region;

typedef Likelihood_Vector                              lvector;
typedef lvector::size_type                             size_type;
typedef meiosis_map                                    mmap;

enum putative_type{SIB              = 0, HSIB    = 1, MZTWIN    = 2,
                   PARENT_OFFSPRING = 3, MARITAL = 4, UNRELATED = 5 };

//Cu:   unrelated
//Ch:   half sib
//Cm:   MZtwin
//Cp:   parent offspring
//Call: total cutpoints
enum cutpoint_type { Cu = 0, Ch = 1, Cm = 2, Cp = 3, Call = 4 };

struct solution_type
{
  solution_type() { function_value = std::numeric_limits<double>::quiet_NaN();
                    mean           = std::numeric_limits<double>::quiet_NaN();
                    variance       = std::numeric_limits<double>::quiet_NaN(); }

  double function_value;
  double mean;
  double variance;
};

inline bool test_location(double d)
{
  d = d - (int) d;
  
  return d < 0.01 || d > 0.99;
}

inline bool is_sib(const ind_id i1, const ind_id i2)
{
  if(i1 == i2) return false;
  if( !i1->parent1() || !i1->parent2() || !i2->parent1() || !i2->parent2() )
    return false;
  
  return (i1->parent1() == i2->parent1() || i1->parent1() == i2->parent2()
       || i1->parent2() == i2->parent1() || i1->parent2() == i2->parent2() );
}

inline bool is_sib(const id_pair i)
{
  return is_sib(i.first, i.second);
}

inline bool is_fsib(const ind_id i1, const ind_id i2)
{
  if(i1 == i2) return false;
  if( !is_sib(i1,i2) ) return false;

  return ((i1->parent1() == i2->parent1() && i1->parent2() == i2->parent2())
       || (i1->parent1() == i2->parent2() && i1->parent2() == i2->parent1()));
}

inline bool is_fsib(const id_pair i)
{
  return is_fsib(i.first, i.second);
}

inline bool is_hsib(const ind_id i1, const ind_id i2)
{
  if(i1 == i2) return false;

  return is_sib(i1,i2) && !is_fsib(i1,i2);
}

inline bool is_hsib(const id_pair i)
{
  return is_hsib(i.first, i.second);
}

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
