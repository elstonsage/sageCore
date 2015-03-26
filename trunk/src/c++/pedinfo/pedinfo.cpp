//============================================================================
// File:      pedinfo.cpp
//                                                                          
// Author:    
//                                                                          
// History:   10/00 Underlying files for generating and viewing statistics
//                  revised or totally rewritten.  Statistics for traits added. 
//                  - Dan Baechle                                                   
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "pedinfo/stats.h"
#include "pedinfo/stats_view.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "pedinfo/pedparse.h"
#include "pedinfo/pedinfo.h"
#include "pedinfo/pedinfo_input.h"

using namespace std;
using namespace SAGE;

pedinfo::pedinfo(int argc, char** argv)
       : APP::SAGEapp(APP::APP_PEDINFO, true, argc, argv)
{
  LSFInit();
}

int pedinfo::main()
{
  // Create input data.
  //
  pedinfo_data data(name, debug());

  print_title(data.info());

  data.input(argc, argv);

  cerrormultistream errors;

  errors.insert(data.errors());

  // - Calculate and write results.
  //
  cout << "Generating statistics....................." << flush;
  
  bool analysis_specified = false;
  
  for( size_t a = 0; a < data.analyses().size(); ++a )
  {
    LSF_ptr<LSFBase> i = data.analyses()[a];
    
    if( !*i ) continue;

    pedinfo_parser  parser;

    analysis_specified = true;
    
    RPED::MP_stats                    mps(errors);
    parser.parse(i, errors, &data.pedigrees());
    pedinfo_parser::trait_list  traits = parser.traits();
    size_t                      trait_count = traits.size();
    
    if(trait_count == 0)
    {
      mps.compute(&data.pedigrees());
    }
    else if(trait_count == 1)
    {
      mps.compute(&data.pedigrees(), traits[0]);
    }
    else
    {
      mps.compute(&data.pedigrees(), traits);
    }
    
    ofstream out_file;  
    out_file.open(parser.file_name().c_str());
    if(!out_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file: pedinfo.out.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    
    print_title(out_file);
    RPED::mp_stats_viewer mp_view(out_file, mps);
    
    // - General stats may be suppressed only if there is a valid trait.
    //
    mp_view.view(parser.show_each_pedigree(), 
                 parser.suppress_general() && trait_count);
  }

  if( !analysis_specified )
  {
    RPED::MP_stats     mps(&data.pedigrees(), errors);
    
    ofstream out_file;  
    out_file.open("pedinfo.out");
    if(!out_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file: pedinfo.out.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
      
    print_title(out_file);
    RPED::mp_stats_viewer mp_view(out_file, mps);
    mp_view.view();
  }

  //cout << "Writing results................." << flush;
  cout << "done." << endl;
  
  cout << endl << "Analysis complete!" << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{
  free(malloc(1));

  pedinfo pedinfo_inst(argc, argv);

  pedinfo_inst.main();

  exit(EXIT_SUCCESS);
}
