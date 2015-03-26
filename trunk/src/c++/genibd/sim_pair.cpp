//==========================================================================
//  File:    sim_pair.cpp
//
//  Author:  Qing Sun
//           Yeunjoo Song
//
//  History: Took simulation pair structure from old MP_mcmc.h
//           & updated to new libraries.                        - yjs Sep 04
//
//  Notes:   This file implements the basic relative pair data structures for
//           ibd simulation method.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_pair.h"
  
namespace SAGE
{

namespace GENIBD
{

sim_relative_pair::sim_relative_pair(size_type m)
{ 
  my_data.resize(m);

  my_member_one = NULL;
  my_member_two = NULL;
}

sim_relative_pair::sim_relative_pair(size_type m,
                                     fmember_const_pointer m1, fmember_const_pointer m2,
                                     pair_type pt)
{
  my_data.resize(m);

  my_member_one = m1;
  my_member_two = m2;

  my_pair_type  = pt;
}

sim_relative_pair::~sim_relative_pair() {}

sim_relative_pair&
sim_relative_pair::operator=(const sim_relative_pair& r)
{
  my_data = r.my_data;

  my_member_one = r.my_member_one;
  my_member_two = r.my_member_two;

  my_pair_type = r.my_pair_type;
  
  return *this;
}

} // end of namespace GENIBD

} // end of namespace SAGE
