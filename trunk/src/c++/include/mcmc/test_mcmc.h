#ifndef TEST_MCMC_H
#define TEST_MCMC_H

//============================================================================
// File:      test_mcmc.h
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation                                   May. 04
//
// Notes:
//
// Copyright (c) 2004 R.C. Elston
// All Rights Reserved
//============================================================================

#include "app/SAGEapp.h"
#include "mcmc/test_mcmc_analysis.h"

namespace SAGE
{

namespace MCMC
{

//----------------------------------------------------------------------------
//  Class:    test_mcmc
//                                                                          
//  Purpose:  Application class for the test_MCMC program.
//                                                                          
//----------------------------------------------------------------------------
//
class test_mcmc : public APP::SAGEapp
{
  public:

    // Constructor/destructor.
    test_mcmc(int argc = 0, char** argv = NULL);
    ~test_mcmc();
    
    // Print program information.
    virtual void  print_help(std::ostream&);

    // Run the application.
    virtual int main();

  protected:
  
    LSF_ptr<RPED::genome_description> build_genome_description(const RPED::RefMPedInfo& minfo,
                                                         LSFBase*           regions,
                                                         bool               multipoint,
                                                         double             distance,
                                                         cerrorstream&      err = sage_cerr);

    void print_analysis_table(test_mcmc_parameters& params, ostream& o);

};

} // end of namespace MCMC

} // end of namespace SAGE

#endif
