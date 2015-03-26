#ifndef  MCMC_ANALYSIS_H
#define  MCMC_ANALYSIS_H

//==========================================================================
//  File:       test_mcmc_analysis.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Notes:      This class is the main driver performing test_mcmc_analysis.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/test_mcmc_parser.h"
#include "mcmc/marker_likelihoods.h"

namespace SAGE
{

namespace MCMC
{

class test_mcmc_analysis
{
  public:

    test_mcmc_analysis(const RPED::RefMultiPedigree& mp,
                       const test_mcmc_parameters&   params,
                       RPED::genome_description*     genome,
                       ostream&                      info,
                       cerrorstream&                 err = sage_cerr);

    ~test_mcmc_analysis();

    bool run_analysis();

    const test_mcmc_parameters*     get_parameters()    const;
    const RPED::RefMultiPedigree*   get_multipedigree() const;

  private:

    // ANALYSIS
    //
    bool  do_analysis();

    void  process_pedigree(const string&      out,
                           const RPED::RefPedigree& rped,  const RPED::genome_description::region_type& r);

    void  process_subpedigree(const string& out,
                              const string&  name,  const RPED::genome_description::region_type& r,
                              const McmcMeiosisMap& mmap);

    void  test_mcmc_engine_components(const RPED::genome_description::region_type& r,
                                      const McmcMeiosisMap& mmap);

    //
    // MEMBER DATA
    //
    const RPED::RefMultiPedigree*     my_multipedigree;
    const test_mcmc_parameters*       my_parameters;
    RPED::genome_description*         my_genome;

    pedigree_region                   my_ped_region;
    
    ostream&                          info;
    cerrorstream                      errors;
};

#include "mcmc/test_mcmc_analysis.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
