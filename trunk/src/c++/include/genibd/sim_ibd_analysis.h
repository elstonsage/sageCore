#ifndef SIM_IBD_ANALYSIS_H
#define SIM_IBD_ANALYSIS_H

//==========================================================================
//  File:    sim_ibd_analysis.h
//
//  Author:  Qing Sun
//
//  History: Version 0.90
//           Maternal & paternal bit split for sib pair done.   - yjs Jun 02
//           Updated to new libraries.                          - yjs Sep 04
//
//  Notes:   This header defines detailed ibd_file implementations related
//           to MCMC ibd simulation program, and provide interfaces to
//           MCMC simulation program.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_storage_ibd.h"
#include "genibd/sim_simulator.h"

namespace SAGE
{

namespace GENIBD
{

class sim_ibd_analysis
{
   public:

    sim_ibd_analysis(const mcmc_parameters& p, cerrorstream& = sage_cerr);

   ~sim_ibd_analysis(); 

    cerrorstream get_errors()  const;
    cerrorstream set_errors(cerrorstream& s);

    bool         built()       const;
    bool         valid()       const;

    IBD*         ibd_adaptor() const;

    //build() must be called before calling compute. It builds up necessary
    //objects for analysis. Only one pedigree is used at a time. If the
    // boolean variable analysis_type=true, do sinlgepoint analysis, 
    //otherwise do multipoint analysis.
    bool         build();

    bool         set_pedigree(const meiosis_map& mmap, const pedigree_region& pr);

    bool         build_ibds();

    size_t       add_pair(mem_pointer i1, mem_pointer i2, pair_type pt);

    bool         compute(ostream& info);

   private:

    //members
    //
    mcmc_parameters             my_params;

    meiosis_map                 my_meiosis_map;

    region_type                 my_region;
    pedigree_region             my_ped_region;

    vector<sim_relative_pair>   my_relative_pairs;

    ibd_mcmc_simulator*         my_simulator;

    IBD*                        my_ibds;

    bool                        my_built;
    bool                        my_valid;
    bool                        my_verbose;
    bool                        my_multipoint;

    cerrorstream                errors;
};

#include "genibd/sim_ibd_analysis.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
