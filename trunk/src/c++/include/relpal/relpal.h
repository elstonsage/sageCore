#ifndef RELPAL_H
#define RELPAL_H

//=============================================================================
// File:    relpal.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                   yjs Mar 07
//
// Notes:   This file contains definition for following data structures.
//            class Relpal : public APP::SAGEapp
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "app/SAGEapp.h"

namespace SAGE   {
namespace RELPAL {

class Relpal : public APP::SAGEapp
{
  public:

    Relpal(int argc=0, char **argv=NULL);
     
    // Run the application
    virtual int main();

  protected:

};

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
