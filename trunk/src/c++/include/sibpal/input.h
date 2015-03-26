#ifndef SIBPAL_INPUT_H
#define SIBPAL_INPUT_H

//****************************************************************************
//* File:      sibpal_input.h                                                *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Aug 01 *
//*                                                                          *
//* Notes:     This header file defines sibpal_data derived from SAGEdata.   *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "data/SAGEdata.h"

namespace SAGE   {
namespace SIBPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     sibpal_data                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines sibpal input data object derived from SAGEdata.      ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class sibpal_data : public APP::SAGE_Simple_Data
{
  public:

    sibpal_data(const string& program_name, bool debug);
    ~sibpal_data();
   
    void input(int argc, char** argv);
   
    // Read sibpal_analysis block
    virtual bool read_analysis();

    const vector< pair<LSF_ptr<LSFBase>, string> >& analysis()      { return my_analysis; }
    const vector< pair<LSF_ptr<LSFBase>, string> >& treg_analysis() { return my_treg_analysis; }
    const vector< pair<LSF_ptr<LSFBase>, string> >& mean_analysis() { return my_mean_analysis; }

    bool  dump_pairs() { return my_dump_pairs; }

  protected:

     vector< pair<LSF_ptr<LSFBase>, string> >  my_analysis;
     vector< pair<LSF_ptr<LSFBase>, string> >  my_treg_analysis;
     vector< pair<LSF_ptr<LSFBase>, string> >  my_mean_analysis;

     bool                                      my_dump_pairs;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
