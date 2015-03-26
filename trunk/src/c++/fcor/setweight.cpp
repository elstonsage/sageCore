//****************************************************************************
//* File:      setweight.cpp                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file defines functions to calculate & set new     *
//*            weight to each pair of pairsetdata.                           *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/setweight.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SetWeight                                                    ~
// ~                                                                         ~
// ~ Purpose:   Defines functions to set weight of pairsetdata.              ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of SetWeight
// ---------------------------------------------------------------------------

SetWeight::SetWeight()
{}

void
SetWeight::set_weight(pairset_by_pedigree_type& pairsets, weight_type w)
{
  if( w == PAIR_WISE )
  {
    pairset_by_pedigree_iterator value     = pairsets.begin();
    pairset_by_pedigree_iterator value_end = pairsets.end();
  
    for( ; value != value_end; ++value )
      set_pair_weight(*value);

    return;
  }
     
  if( w == UNIFORM )
  {
    pairset_by_pedigree_iterator value     = pairsets.begin();
    pairset_by_pedigree_iterator value_end = pairsets.end();
  
    for( ; value != value_end; ++value )
      set_uniform_weight(*value);

    return;
  }

  if( w == MEAN )
  {
    pairset_by_pedigree_iterator value     = pairsets.begin();
    pairset_by_pedigree_iterator value_end = pairsets.end();
  
    for( ; value != value_end; ++value )
      set_mean_weight(*value);

    return;
  }

  return;
}

void
SetWeight::set_pair_weight(pairset_type& pairset)
{
  pairset_iterator value     = pairset.begin();
  pairset_iterator value_end = pairset.end();

  for( ; value != value_end; ++value )
    value->pair_weight[PAIR_WISE] = 1.0;
}  

void
SetWeight::set_uniform_weight(pairset_type& pairset)
{
  size_t m = pairset.size();
  
  pairset_iterator value     = pairset.begin();
  pairset_iterator value_end = pairset.end();

  for( ; value != value_end; ++value )
  {
    if( m )
      value->pair_weight[UNIFORM] = 1.0 / double( m );
    else
      value->pair_weight[UNIFORM] = 0.0;
  }
}

void
SetWeight::set_mean_weight(pairset_type& pairset)
{
  pairset_iterator value     = pairset.begin();
  pairset_iterator value_end = pairset.end();

  for( ; value != value_end; ++value )
  {
    double pair = value->pair_weight[PAIR_WISE];
    double unif = value->pair_weight[UNIFORM];

    if( !isnan(pair) && !isnan(unif) )
      value->pair_weight[MEAN] = (pair + unif) / 2.0;
    else
      value->pair_weight[MEAN] = 0.0;
  }
}  

// end of SetWeight Implementation

} // end of namespace FCOR
} // end of namespace SAGE
