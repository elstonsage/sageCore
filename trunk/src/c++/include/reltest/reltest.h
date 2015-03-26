#ifndef RELTEST_H
#define RELTEST_H

//============================================================================
// File:      reltest.h
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation                                   Jul 03
//
// Notes:
//
// Copyright (c) 2003 R.C. Elston
// All Rights Reserved
//============================================================================

#include "app/SAGEapp.h"

namespace SAGE
{

namespace RELTEST
{

//----------------------------------------------------------------------------
//  Class:    reltest
//                                                                          
//  Purpose:  Application class for the Reltest program.
//                                                                          
//----------------------------------------------------------------------------
//
class reltest : public APP::SAGEapp
{
  public:

    // Constructor/destructor.
    reltest(int argc = 0, char** argv = NULL);
    ~reltest();
    
    // Run the application.
    virtual int main();

  protected:
  
};

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
