#ifndef SETWEIGHT_H
#define SETWEIGHT_H

//****************************************************************************
//* File:      setweight.h                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines functions to calculate & set new     *
//*            weight to each pair of pairsetdata.                           *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/definitions.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SetWeight                                                    ~
// ~                                                                         ~
// ~ Purpose:   Defines functions to set weight of pairsetdata.              ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SetWeight
{
  public:
  
    SetWeight();

    void         set_weight(pairset_by_pedigree_type& pairset, weight_type w);

  private:
    
    void         set_pair_weight   (pairset_type& pairset);
    void         set_uniform_weight(pairset_type& pairset);
    void         set_mean_weight   (pairset_type& pairset);
};    
    
} // end of namespace FCOR
} // end of namespace SAGE

#endif
