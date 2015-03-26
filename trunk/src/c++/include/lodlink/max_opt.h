#ifndef LODLINK_MAX_OPT_H
#define LODLINK_MAX_OPT_H
//============================================================================
// File:      max_opt.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/10/2 - created.                                   djb
//                                                                          
// Notes:     Recipe for optimal likelihood maximization.  Copied from
//            maxfun/maxex_opt.cpp.  See max_opt.cpp for history of orig-
//            inal code
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>
#include "maxfunapi/maxfunapi.h"

using std::ostream;

namespace SAGE
{

namespace LODLINK
{
  
void  set_maxfun_configuration(MAXFUN::SequenceCfg& configuration);

}
}

#endif

