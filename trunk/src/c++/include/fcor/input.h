#ifndef FCORINPUT_H
#define FCORINPUT_H

//****************************************************************************
//* File:      input.h                                                       *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Aug 01 *
//*                                                                          *
//* Notes:     This header file defines fcor_data derived from SAGEdata.     *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "data/SAGEdata.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     FcorData                                                     ~
// ~                                                                         ~
// ~ Purpose:   Defines fcor input data object derived from SAGEdata.        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FcorData : public APP::SAGE_Data
{
  public:

    FcorData(const string& program_name, bool debug = false);
    ~FcorData();
   
    void input(int argc, char** argv);
   
    // Read fcor_analysis block
    virtual bool read_analysis();

    const vector<LSF_ptr<LSFBase> >& analysis() { return my_analysis; }

  protected:

     vector<LSF_ptr<LSFBase> >  my_analysis;
};

} // end of namespace FCOR
} // end of namespace SAGE

#endif
