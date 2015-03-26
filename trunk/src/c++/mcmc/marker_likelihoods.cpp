//==========================================================================
//  File:    marker_likelihoods.cpp
//
//  Author:  Geoff Wedig
//
//  History: Version 0.90
//           1.0 Updated to new libraries                        yjs May. 04
//
//  Notes:
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/marker_likelihoods.h"

namespace SAGE
{

namespace MCMC
{

#define DEBUG(x)

// =============================
// marker_likelihood_calculator 
// =============================

marker_likelihood_calculator::marker_likelihood_calculator(const pedigree_region&  pr,
                                                           const McmcMeiosisMap& mmap,
                                                           mcmc_data_accessor&     dat)
  : my_cache_data(pr.get_region().locus_count()),
    my_graphs(pr.get_region().locus_count()),
    my_region(pr),
    my_data(&dat)
{
  // Build the initial trees.

  size_t marker_count = pr.get_region().locus_count();

  for( size_t m = 0; m < marker_count; ++m )
  {
    const MLOCUS::inheritance_model& imodel = pr[m];

    my_graphs[m].setup(mmap, imodel);

    if( !my_data->is_valid_locus(m) ) continue;
  }
}

void
marker_likelihood_calculator::dump_graphs(ostream& o) const
{
  for(size_t i = 0; i < my_region.get_region().locus_count(); ++i)
  {
     o << "Marker " << i << endl;

     my_graphs[i].dump_graph(o);
  }
}

} // end of namespace MCMC

} // end of namespace SAGE
