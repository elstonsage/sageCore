#ifndef RPEDTEST_INPUT_H
#define RPEDTEST_INPUT_H

//============================================================================
// File:      rpedtest_input.h
//                                                                          
// Author:    
//                                                                          
// History:   3-27-01 modified to add test of genome_description class.  - djb                                                   
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <iostream>
#include <fstream>
#include <list>
#include "mlocus/imodel.h"
#include "mlocus/mfile.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "rped/genome_description.h"
#include "error/errorstream.h"
#include "LSF/LSFsymbol.h"

namespace SAGE {
namespace RPED {

class rptest_data
{
public:

  rptest_data() { }

  ~rptest_data();

  bool input(int argc, char* argv[]);

// Get the data members

  SymbolTable*        parameters()       const;
  AttrVal             parameter (string) const;
 
  const RefMultiPedigree*   pedigrees()        const;

  cerrormultistream errors;
  cerrormultistream information;

  ofstream messages;
  
private:

  // File Readers

  void read_par();                     // Read Parameter File
  void read_ldf(char sep);             // Read Locus Description File
  void read_fdf();                     // Read Family Data File
  void read_gdf();                     // Read Genome Description File

  int    argc;
  char** argv;
  
  LSF_ptr<SymbolTable>         params;

  RefMultiPedigree             pedigree_data;
};

} // End namespace RPED
} // End namespace SAGE

#include "rped/rpedtest_input.ipp"

#endif
