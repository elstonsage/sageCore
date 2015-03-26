//==========================================================================
//  File:    mcmc_data_access.cpp
//
//  Author:  Geoff Wedig
//           Yeunjoo Song
//
//  History: Version 0.90
//           1.0 Changed the class name MCMCDataStrategy
//               to mcmc_data_access & updated to new libraries  yjs May. 04
//
//  Notes:   This file implements accesses to individual's raw data
//           for MCMC IBD simulatior.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/mcmc_data_accessor.h"

namespace SAGE
{

namespace MCMC
{

mcmc_data_accessor::mcmc_data_accessor(const McmcMeiosisMap& mmap, size_t lc)
                  : my_mcmc_meiosis_map(mmap), my_locus_count(lc)
{ 
  my_valid_loci.resize(locus_count(), false);

  my_indicator.resize(locus_count(), bit_field(my_mcmc_meiosis_map.get_nonfounder_count() * 2));
}

void
mcmc_data_accessor::dump_accessor(ostream& o) const
{
  o << endl << "mcmc_data_accessor dump :" << endl;
  o << "locus_count        = "  << my_locus_count << endl;
  o << "indicator size     = "  << my_indicator.size() << endl;
  o << "valid_loci size    = "  << my_valid_loci.size() << endl << endl;
}

} // end of namespace MCMC

} // end of namespace SAGE
