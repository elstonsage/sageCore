#ifndef MCMC_INPUT_H
#define MCMC_INPUT_H

//==========================================================================
//  File:     test_mcmc_input.h
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.                                May. 04
//
//  Notes:
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "data/SAGEdata.h"

namespace SAGE
{

namespace MCMC
{

class test_mcmc_data : public APP::SAGE_Data
{
  public:

    test_mcmc_data(const string& program_name, bool debug);

    ~test_mcmc_data();

    bool input(int argc, char** argv);

    virtual bool read_analysis();
    
    const vector<LSF_ptr<LSFBase> >& get_analysis() const;

  private:

    vector<LSF_ptr<LSFBase> >    my_analysis;

};

// ================
// Inline Functions
// ================

inline
const vector<LSF_ptr<LSFBase> >& test_mcmc_data::get_analysis() const
{
  return my_analysis;
}


} // end of namespace MCMC

} // end of namespace SAGE

#endif
