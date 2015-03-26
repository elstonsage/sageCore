//============================================================================
// File:      maxfun_utilities.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/16/3 - created.                       -djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "maxfun/maxfun_utilities.h"

namespace SAGE
{

// - Return the total number of independent parameters, both those involved in
//   a functional relationship and those not involved in a functional relationship
//
size_t
ind_param_count(Maxfun& max)
{
  size_t count = 0;
  for(size_t i = 0; i < (size_t) max.nt(); ++i)
  {
    if(max.istin(i) == 1 || max.istin(i) == 2)
    {
      ++count;
    }
  }
  
  return  count;
}

}

