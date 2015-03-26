#ifndef GELIMTEST_INPUT_H
#define GELIMTEST_INPUT_H

#include <iostream>
#include <fstream>
#include <list>
#include "mlocus/imodel.h"
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "LSF/LSFsymbol.h"

namespace SAGE {

class segreg_data
{
public:

  segreg_data() { }

  ~segreg_data();

  bool input(int argc, char* argv[]);

// Get the data members

  SymbolTable*        parameters()       const;
  AttrVal             parameter (string) const;
 
  const RPED::RefMultiPedigree*   pedigrees()        const;

  SAGE::cerrormultistream errors;
  SAGE::cerrormultistream information;

  ofstream messages;
  
private:

  // File Readers

  void read_par ();                     // Read Parameter File
  void read_ldf (char sep);             // Read Locus Description
  void read_fdf ();                     // Read Family Data File

  int    argc;
  char** argv;
  
  LSF_ptr<SymbolTable>         params;

  RPED::RefMultiPedigree             pedigree_data;
};

#include "segreg/test_segreg_input.ipp"

} // End sage namespace

#endif
