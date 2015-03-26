//============================================================================
// File:      max_opt.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 12/10/2.                                                   
//                                                                          
// Notes:     Recipe for optimal likelihood maximization.  Copied from
//            maxfun/maxex_opt.cpp and modified.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

/* MAXEX-EASY:  SAMPLE DRIVER PROGRAM FOR MAXFUN
 *
 * Last Modified:
 *  1-NOV-2000 - Geoff Wedig  - Using the 'optimum' strategy.
 *  7-JUN-1999 - Kevin Jacobs - Updated for new Fortran runtime
 * 13-APR-1998 - Kevin Jacobs - C++ conversion
 * 24-JAN-1996 - Kevin Jacobs - updated
 *  7-DEC-1995 - Kevin Jacobs - created
*/

#include "lodlink/max_opt.h"

namespace SAGE
{

namespace LODLINK
{

void
set_maxfun_configuration(MAXFUN::SequenceCfg& configuration)
{
  int  failure = configuration.addRunCfg(MAXFUN::RunCfg::DIRECT_WITHOUT, 3);
  assert(! failure);
  MAXFUN::RunCfg&  first_RunCfg = configuration.getLatestRunCfg();
  first_RunCfg.epsilon1 = 1e-2;
  first_RunCfg.epsilon2 = 1e-3;
  
  failure = configuration.addRunCfg(MAXFUN::RunCfg::VAR_METRIC_ESTIMATE, 20);
  assert(! failure);
  MAXFUN::RunCfg&  second_RunCfg = configuration.getLatestRunCfg();
  second_RunCfg.epsilon1 = 1e-4;
  second_RunCfg.epsilon2 = 1e-7;
  second_RunCfg.var_cov = MAXFUN::RunCfg::FINAL;
}

}
}

