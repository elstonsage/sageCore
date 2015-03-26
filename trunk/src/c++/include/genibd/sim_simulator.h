#ifndef  GENIBD_MCMC_SIMULATOR_H
#define  GENIBD_MCMC_SIMULATOR_H

//==========================================================================
//  File:       sim_simulator.h
//
//  Author:     Qing Sun
//
//  History:    Version 0.90
//              Updated to new libraries.                       - yjs Oct 04
//
//  Notes:      This header defines interfaces to IBD simulation object.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/sim_pair.h"

namespace SAGE
{

namespace GENIBD
{

class ibd_mcmc_simulator : public mcmc_simulator
{
  public:

    ibd_mcmc_simulator(const mcmc_meiosis_map& ped, const pedigree_region& pr,
                       mcmc_parameters*        par, cerrorstream&          err);

   ~ibd_mcmc_simulator();

    virtual bool start(ostream& info);

    vector<sim_relative_pair>*  set_pairs(vector<sim_relative_pair>* rp);
    vector<sim_relative_pair>*  get_pairs() const;
    
  private:

    struct tree_data
    {
      bit_field         pattern;
      size_t            last_used;
      vector<size_t>    tree;
    };
    
    int             set_init_ibd_sc();

    void            singlepoint_ibd_sharing(int offset);
    void            multipoint_ibd_sharing(int offset);

    void            ibd_sharing(int offset, int mid);

    void            do_multipoint_analysis(cerrorstream& o);
    void            do_singlepoint_analysis(cerrorstream& o);

    vector<sim_relative_pair>*   my_relative_pairs;

    vector<tree_data>            my_trees;
};

#include "genibd/sim_simulator.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
