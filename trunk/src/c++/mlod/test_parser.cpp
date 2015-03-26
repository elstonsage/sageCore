//============================================================================
// File:      test_parser.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 5/28/92.                                                   
//                                                                          
// Notes:     Tests the mlod parser class.
//    
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================
#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "mlod/parser.h"
#include "mlod/data.h"

using namespace std;
using namespace SAGE;
using namespace MLOD;

int main(int argc, char** argv)
{
  if (argc != 6)
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <trait> <locus> <genome>"
         << endl << endl
         << "Command line parameters:"             << endl
         << "  parameters   - Parameter File"      << endl
         << "  pedigree     - Pedigree Data File"  << endl
         << "  trait        - Trait Locus Description File" << endl
         << "  locus        - Marker Locus Description File" << endl
         << "  genome       - Genome Description File" << endl
         << endl << endl;
    exit(EXIT_FAILURE);
  }
  
  LSFInit();
  
  Data  mlod_data("Test_Parser", false);
  
  mlod_data.input(argc, argv);
  
  // - Parse segreg_analysis blocks.
  //
  for(int i = 0; i < (int)mlod_data.get_analyses().size(); ++i)
  {
    mlod_data.get_analyses()[i].print_analysis_table(mlod_data.screen());
  }

  exit(EXIT_SUCCESS);
}
