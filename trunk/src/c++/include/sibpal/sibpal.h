#ifndef SIBPAL_H
#define SIBPAL_H

//=============================================================================
// File:    sibpal.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//            class Sibpal : public APP::SAGEapp
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include <limits>
#include <iostream>
#include "app/SAGEapp.h"

namespace SAGE   {
namespace SIBPAL {

class Sibpal : public APP::SAGEapp
{
  public:

    Sibpal(int argc=0, char **argv=NULL);
     
    // Run the application
    virtual int main();

  protected:

    LSF_ptr<SymbolTable>     syms;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
