#ifndef GENIBD_TRANS_CALCULATOR_H
#define GENIBD_TRANS_CALCULATOR_H

//============================================================================
// File:      trans_calculator.h
//
// Author:    Dan Baechle
//
// History:   10/14/2 - created.                                 djb
//            Modified for GENIBD                                yjs Apr 2004
//
// Notes:     defines a class to calculate transition probabilities.
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

#include "genibd/definitions.h"

namespace SAGE
{

namespace GENIBD
{

class trans_calculator
{
  public:

    trans_calculator();
  
    double transition(const ind_genotype& mom, 
                      const ind_genotype& dad,
                      const ind_genotype& kid ) const;

  private:

    static bool  kd(const allele& one, const allele& two);    // Kronecker delta.

    double transmission(const ind_genotype& p, MPED::SexCode s, const allele& ka) const;

    trans_calculator(const trans_calculator& other);
    trans_calculator&  operator=(const trans_calculator& other);
  
    //Data members.
};

#include "genibd/trans_calculator.ipp"

} // end of GENIBD namespace

} // end of SAGE namespace

#endif
