#ifndef TEST_OUTPUT_H
#define TEST_OUTPUT_H
//============================================================================
// File:      test_output.h
//
// Author:    Yeunjoo Song
//
// History:   8/07/01 - created.                            yjs
//
// Notes:     Defines test_output application object.
//
// Copyright (c) 2001 R.C. Elston
// All Rights Reserved
//============================================================================


#include <limits>
#include <iostream>
#include "app/SAGEapp.h"
#include "app/aparser.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "LSF/LSFsymbol.h"

namespace SAGE
{
namespace SEGREG
{

//----------------------------------------------------------------------------
//  Class:     segreg_test_output
//
//  Purpose:   Defines segreg_test_output application derived from SAGEapp.
//
//----------------------------------------------------------------------------

class test_output : public SAGEapp
{
  public:

    test_output(const string& program_name, int argc=0, char **argv=NULL);
     
    // Print program information
    virtual void print_help(std::ostream& );
    
    // Run the application
    virtual int main();

    bool init_output_streams();

  protected:

    LSF_ptr<SymbolTable>  syms;

    ofstream                 info_file;
    SAGE::cerrormultistream  errors;
};

} // end of namespace REGX
} // end of namespace SAGE

#endif
