#ifndef RELPAL_INPUT_H
#define RELPAL_INPUT_H

//=============================================================================
// File:      input.h
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation                             yjs Mar 07 *
//
// Notes:     This header file defines RelpalData derived from SAGEdata.
//
// Copyright (c) 2007 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "data/SAGEdata.h"

namespace SAGE   {
namespace RELPAL {

class relpal_data : public APP::SAGE_Simple_Data
{
  public:

    relpal_data(const string& program_name, bool debug);
    ~relpal_data();
   
    void input(int argc, char** argv);
   
    // Read relpal_analysis block
    virtual bool read_analysis();

    const vector< LSF_ptr<LSFBase> >& analysis() { return my_analysis; }

  protected:

     vector< LSF_ptr<LSFBase> >  my_analysis;
};

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
