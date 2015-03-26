#ifndef LODPAL_INPUT_H
#define LODPAL_INPUT_H

//****************************************************************************
//* File:      lodpal_input.h                                                *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Aug 01 *
//*                                                                          *
//* Notes:     This header file defines lodpal_data derived from SAGEdata.   *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "data/SAGEdata.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     lodpal_data                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines lodpal input data object derived from SAGEdata.      ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class lodpal_data : public APP::SAGE_Simple_Data
{
  public:

    lodpal_data(const string& program_name, bool debug);
    ~lodpal_data();
   
    void input(int argc, char** argv);
   
    // Read lodpal_analysis block
    virtual bool read_analysis();

    const vector<LSF_ptr<LSFBase> >& analysis() { return my_analysis; }

    bool  dump_pairs() { return my_dump_pairs; }

  protected:

     vector<LSF_ptr<LSFBase> >  my_analysis;

     bool                       my_dump_pairs;
};

} // end of namespace LODPAL
} // end of namespace SAGE 

#endif
