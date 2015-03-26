#ifndef LODPAL_H
#define LODPAL_H

//****************************************************************************
//* File:      lodpal.h                                                      *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation.                yjs Jan 01 *
//*                                                                          *
//* Notes:     This header file defines lodpal app. derived from SAGEapp.    *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "app/aparser.h"
#include "app/SAGEapp.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     Lodpal                                                       ~
// ~                                                                         ~
// ~ Purpose:   Defines lodpal application derived from SAGEapp.             ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Lodpal : public APP::SAGEapp
{
  public:
  
    Lodpal(const string& program_name, int argc=0, char **argv=NULL);
     
    // Run the application
    virtual int main();

  protected:

    LSF_ptr<SymbolTable>  syms;
};

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
