#ifndef TESTPEELER_INPUT_H
#define TESTPEELER_INPUT_H

#include <iostream>
#include <fstream>
#include <list>
#include "rped/rped.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "LSF/LSFsymbol.h"

namespace SAGE
{

class test_peeler_data
{
public:

  test_peeler_data() { }

  ~test_peeler_data();

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

#include "peeling/test_peeler_input.ipp"

}

#endif
