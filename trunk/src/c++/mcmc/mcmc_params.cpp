//==========================================================================
//  File:       mcmc_params.cpp
//
//  Author:     Qing Sun
//              Geoff Wedig
//              Yeunjoo Song
//
//  History:    Initial implementation.                        Aug. 28, 1998
//              1.0 Completely rewritten so makes sense        Feb   4, 2000
//              2.0 Updated to new libraries.                  yjs May. 2004
//
//  Notes:      This file implements a parameter for mcmc analysis.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/mcmc_params.h"

namespace SAGE
{

namespace MCMC
{

mcmc_parameters::mcmc_parameters()
{
  my_multipoint            = true;

  my_max_tunnel            = DEFAULT_TT_MAX_TRANS;
  my_T[0]                  = DEFAULT_T0_WEIGHT;
  my_T[1]                  = DEFAULT_T1_WEIGHT;
  my_T[2]                  = DEFAULT_T2_WEIGHT;
  my_single_marker         = 0.0;
  my_local_marker          = DEFAULT_LOCAL_WEIGHT;
  my_local_individual      = 0.0; // This needs testing!

  my_dememorization_step   = DEFAULT_DEMEMO_STEPS;
  my_simulation_step       = DEFAULT_MCMC_STEPS;
  my_batch_count           = DEFAULT_BATCH_COUNT;

  my_use_factor            = true;

  my_dememorization_factor = 15;
  my_simulation_factor     = 150;
  my_batch_factor          = 30;
  
  my_random_seed           = 0;
}

mcmc_parameters::mcmc_parameters(const mcmc_parameters& p)
{
  my_multipoint            = p.my_multipoint;

  my_max_tunnel            = p.my_max_tunnel;
  my_T[0]                  = p.my_T[0];
  my_T[1]                  = p.my_T[1];
  my_T[2]                  = p.my_T[2];
  my_single_marker         = p.my_single_marker;
  my_local_marker          = p.my_local_marker;
  my_local_individual      = p.my_local_individual;

  my_dememorization_step   = p.my_dememorization_step;
  my_simulation_step       = p.my_simulation_step;
  my_batch_count           = p.my_batch_count;

  my_use_factor            = p.my_use_factor;

  my_dememorization_factor = p.my_dememorization_factor;
  my_simulation_factor     = p.my_simulation_factor;
  my_batch_factor          = p.my_batch_factor;
  
  my_random_seed           = p.my_random_seed;
}

mcmc_parameters&
mcmc_parameters::operator=(const mcmc_parameters& p)
{
  my_multipoint            = p.my_multipoint;

  my_max_tunnel            = p.my_max_tunnel;
  my_T[0]                  = p.my_T[0];
  my_T[1]                  = p.my_T[1];
  my_T[2]                  = p.my_T[2];
  my_single_marker         = p.my_single_marker;
  my_local_marker          = p.my_local_marker;
  my_local_individual      = p.my_local_individual;

  my_dememorization_step   = p.my_dememorization_step;
  my_simulation_step       = p.my_simulation_step;
  my_batch_count           = p.my_batch_count;

  my_use_factor            = p.my_use_factor;

  my_dememorization_factor = p.my_dememorization_factor;
  my_simulation_factor     = p.my_simulation_factor;
  my_batch_factor          = p.my_batch_factor;
  
  my_random_seed           = p.my_random_seed;

  return *this;
}

mcmc_parameters::~mcmc_parameters()
{}

void
mcmc_parameters::dump_parameters(ostream& out) const
{
  SAGE::cerrorstream o(out);
  
  o.prefix("MCMC : ");

  if( is_multipoint() )
    o << "Multipoint";
  else
    o << "Singlepoint";

  o << " MCMC analysis ";

  o << ':' << endl;

  o.prefix("  * ");

  o << "max_tunnel:            " << my_max_tunnel << endl;
  o << "T[0]:                  " << my_T[0] << endl;
  o << "T[1]:                  " << my_T[1] << endl;
  o << "T[2]:                  " << my_T[2] << endl;
  o << "single_marker:         " << my_single_marker << endl;
  o << "local_marker:          " << my_local_marker << endl;
  o << "local_individual:      " << my_local_individual << endl;

  o << "dememorization_step:   " << my_dememorization_step << endl;
  o << "simulation_step:       " << my_simulation_step << endl;
  o << "batch_count:           " << my_batch_count << endl;

  o << "use_factor:            " << my_use_factor << endl;

  o << "dememorization_factor: " << my_dememorization_factor << endl;
  o << "simulation_factor:     " << my_simulation_factor << endl;
  o << "batch_factor:          " << my_batch_factor << endl;
  
  if(my_random_seed)
    o << "random seed:           " << my_random_seed << endl;

  o << endl;
}

} // end of namespace MCMC

} // end of namespace SAGE

