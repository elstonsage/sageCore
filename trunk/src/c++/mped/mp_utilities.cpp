//===================================================================
//
//  File:	mp_utilities.cpp
//
//  Author:	Stephen Gross
//
//  History:	sag Initial implementation		Aug 01 2001
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//===================================================================

#include <algorithm>
#include <functional>
#include <boost/bind.hpp>
#include "mped/mp_utilities.h"

using namespace std;

namespace SAGE {
namespace MPED {

//===================================================================
// Loop testing
//===================================================================
bool 
mp_utilities::has_loops( const subpedigree_base& sped)
{
  size_t number_of_founders = 
      count_if(sped.member_begin(), sped.member_end(),
               mem_fun_ref(&member_base::is_founder));

  return (sped.family_count() != (number_of_founders-1));
}

bool
mp_utilities::no_loops( const subpedigree_base& sped)
{
  return !has_loops(sped);
}

//===================================================================
// Chain testing
//===================================================================
bool
mp_utilities::has_chains(const subpedigree_base& sped)
{
  for(member_const_iterator i  = sped.member_begin();
                            i != sped.member_end();   ++i)
  {
    if(i->mate_count() > 1)
    {
      for(mate_const_iterator m = i->mate_begin(); m != i->mate_end(); ++m)
      {
        if(m->mate().mate_count() > 1) return true;
      } // End of mate loop
    }
  } // End of member loop
  return false;
}

bool
mp_utilities::no_chains(const subpedigree_base& sped)
{
  return !has_chains(sped);
}

//===================================================================
// Max cluster size
//===================================================================
int
mp_utilities::max_cluster_size(const subpedigree_base& sped)
{
  size_t my_max = max_element(sped.member_begin(), sped.member_end(),
                              boost::bind(less<uint>(),
                                          boost::bind(&member_base::mate_count,_1),
                                          boost::bind(&member_base::mate_count,_2)))->mate_count();

  return my_max;
}

//===================================================================
// getMaxSibshipSize(...)
//===================================================================
int
mp_utilities::getMaxSibshipSize(const multipedigree_base & RMP)
{
  int max_sibship_size = 0;

  for(pedigree_base::pedigree_const_iterator pedigree_itr  = RMP.pedigree_begin ();
                                             pedigree_itr != RMP.pedigree_end   (); ++pedigree_itr)
  {
    for(family_base::family_const_iterator fam = pedigree_itr->family_begin(); 
                                           fam !=  pedigree_itr->family_end(); ++fam)
    {
      if((int) fam->offspring_count() > max_sibship_size)
        max_sibship_size = fam->offspring_count();
    }
  }

  return max_sibship_size;
}

} // End namespace MPED
} // End namespace SAGE
