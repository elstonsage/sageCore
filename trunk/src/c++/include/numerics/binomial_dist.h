#ifndef BINOMIAL_DIST_H
#define BINOMIAL_DIST_H

//========================================================================
//
//  File:	binomial_dist.h
//
//  Author:	Stephen Gross
//
//  Copyright (c) 2001, R. C. Elston
//
//========================================================================

#include <cstddef>
#include <cmath>
#include "numerics/cephes.h"

namespace SAGE     {
namespace NUMERICS {

static double bin_prob(size_t q, double p, size_t x) 
{
  double bin_prob     = 0.0;
  double chooses_term = 0.0;

  if(x > q) 
    bin_prob = 0;
  else
  {
    chooses_term = choose(q,x);
    bin_prob     = chooses_term * pow((double)p,(double)x) *
                   pow((float)(1-p),(float)(q-x));
  }
  return bin_prob;
}

} // End namespace NUMERICS
} // End namespace SAGE

#endif
