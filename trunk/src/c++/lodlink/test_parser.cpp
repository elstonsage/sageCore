//============================================================================
// File:      test_parser.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 9/4/2                                                   
//                                                                          
// Notes:     Tests the lodlink parser class.
//    
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
#include "lodlink/parser.h"
#include "data/SAGEdata.h"

using namespace std;
using namespace SAGE;
using namespace LODLINK;

int main(int argc, char* argv[])
{
  // ============  The preliminaries.  =============
  //
  if (! (argc == 4 || argc == 5))
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> <locus> <trait>\n\n"
         << "Command line parameters:\n"
         << "  parameters   - Parameter File\n"
         << "  pedigree     - Pedigree Data File\n"       
         << "  trait        - Trait File\n"
         << "\n" << endl;
    exit(EXIT_FAILURE);
  }
  
  LSFInit();  
  
  // - Read data.
  //
  lodlink_data  my_data("testparser", false);

  my_data.input(argc, argv);
  
  
  const RPED::RefMultiPedigree*  my_mped = &(my_data.pedigrees());
  assert(my_mped != 0);

  cout << "Generating statistics....................." << flush;

  parser  my_parser(my_mped, my_data, cout, my_data.errors());
  
  ofstream out_file;
  out_file.open("testparser.out");
  
  if(! out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testparser.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }  
 
  // - Parse each analysis block.
  //
  for(size_t a = 0; a < my_data.analyses().size(); ++a)
  {
    LSF_ptr<LSFBase> analysis_ptr = my_data.analyses()[a];
    assert(analysis_ptr != 0);

    size_t  id = my_parser.analysis_id() + 1;
  
    my_data.errors() << priority(error) << endl;
    my_data.errors() << priority(error) << "---------- " << id << " New Analysis Block " << id << " ---------- " << endl;
    my_data.errors() << priority(error) << endl;
    
        
    out_file << endl;
    out_file << "---------- " << id << " New Analysis Block " << id << " ---------- " << endl;
    out_file << endl;
    
    my_parser.parse(analysis_ptr);
    
    out_file << my_parser.user_instructions() << endl;
  }
}
