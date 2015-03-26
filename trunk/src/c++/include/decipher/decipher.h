#ifndef DECIPHER_H
#define DECIPHER_H

//============================================================================
// File:      decipher.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/11/4 - created.                                   djb
//            4/25/5 - changed program name to decipher           djb
//                                                                          
// Notes:      
//
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "LSF/LSFinit.h"
#include "rped/rped.h"
#include "app/SAGEapp.h"
#include "app/SAGEapp_version_bank.h"
#include "mped/mp_utilities.h"
#include "mlocus/penmodel.h"
#include "decipher/parser.h"
#include "decipher/analysis.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace RPED;
using namespace APP;
using namespace MPED;
using namespace MLOCUS;


class decipher : public SAGEapp
{
  public:

    decipher(int argc=0, char **argv=NULL);
   
    // Run the application
    virtual int main();
    
  private:
    
};

}
} 

#endif

