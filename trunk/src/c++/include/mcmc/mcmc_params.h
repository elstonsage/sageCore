#ifndef  MCMC_PARAM_H
#define  MCMC_PARAM_H

//==========================================================================
//  File:       mcmc_params.h
//
//  Author:     Qing Sun
//              Geoff Wedig
//
//  History:    Version 0.90                                   Aug. 28, 1998
//              1.0 Completely rewritten so makes sense        Feb   4, 2000
//              2.0 Updated to new libraries.                  yjs Aug. 2004
//
//  Notes:      This header gives basic parameters that are used in both
//              mcmc_ibd simulation program and mcmc haplotyping program.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/definitions.h"

namespace SAGE
{

namespace MCMC
{

class mcmc_parameters
{
  public:  

    mcmc_parameters();

    mcmc_parameters(const mcmc_parameters&);

    mcmc_parameters& operator=(const mcmc_parameters&);

    ~mcmc_parameters();

    void      dump_parameters(ostream& out = cout) const;

    // Return fuctions

    bool      is_multipoint               ()          const;

    size_t    get_max_tunnel              ()          const;
    double    get_single_marker_tunneling ()          const;
    double    get_local_marker            ()          const;
    double    get_local_individual        ()          const;
    double    get_transition_weight       (size_t)    const;

    long      get_dememorization_step     ()          const;
    long      get_simulation_step         ()          const;
    long      get_batch_count             ()          const;
                                                  
    bool      get_use_factor              ()          const;
                                                  
    double    get_dememorization_factor   ()          const;
    double    get_simulation_factor       ()          const;
    double    get_batch_factor            ()          const;

    unsigned long get_random_seed         ()          const;

    // Modification functions

    bool set_multipoint              (bool m);

    bool set_max_tunnel              (size_t p);
    bool set_single_marker_tunneling (double p);
    bool set_local_marker            (double p);
    bool set_local_individual        (double p);
    bool set_transition_weight       (size_t t, double p);

    bool set_dememorization_step     (long p);
    bool set_simulation_step         (long p);
    bool set_batch_count             (long p);
                                     
    bool set_use_factor              (bool b = true);

    bool set_dememorization_factor   (double p); 
    bool set_simulation_factor       (double p);
    bool set_batch_factor            (double p);

    bool set_random_seed             (unsigned long p);

    void normalize_weights           ();

  protected:

    bool      my_multipoint;

    size_t    my_max_tunnel;
    double    my_single_marker;
    double    my_local_marker;
    double    my_local_individual;
    double    my_T[3];

    long      my_dememorization_step;  
    long      my_simulation_step;
    long      my_batch_count;

    bool      my_use_factor;

    double    my_dememorization_factor;
    double    my_simulation_factor;
    double    my_batch_factor;
    
    unsigned long my_random_seed;
};

#include "mcmc/mcmc_params.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
