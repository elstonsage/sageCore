#ifndef LODLINK_H
#define LODLINK_H
//============================================================================
// File:      lodlink.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 9/20/2
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <iostream>
#include "LSF/LSFinit.h"
#include "rped/rped.h"
#include "app/SAGEapp.h"
#include "app/SAGEapp_version_bank.h"
#include "mped/mp_utilities.h"
#include "lodlink/input.h"
#include "lodlink/parser.h"
#include "lodlink/analysis.h"
#include "lodlink/results.h"

namespace SAGE
{

namespace LODLINK
{

void build_headers();

class lodlink : public APP::SAGEapp
{
  public:

    lodlink(int argc=0, char **argv=NULL);
   
    // Print program information
    virtual void print_title(std::ostream& );
    
    // Run the application
    virtual int main();
    
  private:
    static void  check_for_loops(const RPED::RefMultiPedigree& mped, cerrorstream& errors);
    static void  check_allele_freqs(const RPED::RefMultiPedigree& mped, cerrorstream& errors);
};

}
} 

#endif
