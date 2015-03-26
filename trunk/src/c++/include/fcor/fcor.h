#ifndef FCOR_H
#define FCOR_H

//****************************************************************************
//* File:      fcor.h                                                        *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                  Jan 31 01 *
//*                                                                          *
//* Notes:     This header file defines fcor app. derived from SAGEapp.      *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "app/SAGEapp.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     Fcor                                                         ~
// ~                                                                         ~
// ~ Purpose:   Defines fcor application derived from SAGEapp.               ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Fcor : public APP::SAGEapp
{
  public:

    Fcor(int argc=0, char **argv=NULL);
   
    // Run the application
    virtual int main();

  protected:
  
};

} // end of namespace FCOR
} // end of namespace SAGE

#endif
