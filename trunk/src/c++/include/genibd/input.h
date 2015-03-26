#ifndef GENIBD_INPUT_H
#define GENIBD_INPUT_H

//==========================================================================
//  File:     input.h
//
//  Author:   Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:  0.01 gcw Initial Interface design                 Jun 15, 1998
//            0.1  gcw Release version                          Jul  8, 1998
//            2.0  yjs Updated to new libraries.                Nov.    2004
//
//  Notes:    GENIBD input routines
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "data/SAGEdata.h"

namespace SAGE
{

namespace GENIBD
{

class genibd_data : public APP::SAGE_Data
{
  public:

    genibd_data(const string& program_name, bool debug);

    ~genibd_data();

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
const vector<LSF_ptr<LSFBase> >& genibd_data::get_analysis() const
{
  return my_analysis;
}


} // end of namespace GENIBD

} // end of namespace SAGE

#endif
