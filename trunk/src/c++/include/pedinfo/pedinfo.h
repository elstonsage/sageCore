#ifndef PEDINFO_H
#define PEDINFO_H
//============================================================================
// File:     pedinfo.h
//                                                                          
// Author:   
//
// Program                                                                          
// History:   /  /2 - Rewritten almost entirely.                       -djb
//           6/19/3 - Added listing of individuals w. multiple mates.  -djb
//           6/23/3 - Added listing of consanguineous mating pairs.    -djb
//           6/25/3 - Inheritance vector bits to ind. pedigree stats.  -djb
//           6/25/3 - Changed output so that number of generations
//                    for pedigrees w. non-marriage loops is given
//                    as 'undet.'                                      -djb
//           6/30/3 - fixed 'off by one' errors in log histogram.      -djb
//           9/18/3 - added several stats to trait specific output.    -djb
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"

namespace SAGE
{


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     pedinfo                                                      ~
// ~                                                                         ~
// ~ Purpose:   Defines pedinfo application derived from SAGEapp.            ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class pedinfo : public APP::SAGEapp
{
  public:

    pedinfo(int argc=0, char **argv=NULL);
   
    // Run the application
    virtual int main();
};


} // end of namespace SAGE 

#endif
