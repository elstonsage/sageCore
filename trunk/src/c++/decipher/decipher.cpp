//============================================================================
// File:      decipher.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/11/4  created.                             djb
//            9/22/6  changed file name to decipher.       djb
//                                                                          
// Notes:     decipher class represents the decipher program.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/decipher.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

//============================================================================
// IMPLEMENTATION:  decipher
//============================================================================
//
decipher::decipher(int argc, char** argv)
       : SAGEapp(APP_DECIPHER, true, argc, argv)
{
  LSFInit();
}

int decipher::main()
{
  // - Read data.
  //
  decipher_data  my_data(name, debug());

  print_title(my_data.info());
  my_data.input(argc, argv);  
  
  const RefMultiPedigree*  my_mped = &(my_data.pedigrees());
  assert(my_mped != 0);

  parser  my_parser(my_mped, my_data, cout, my_data.errors(), my_data.genome()); 
  
  // - Do the analysis specified in each parameter file analysis block.
  //
  bool  analysis_specified = false;
  for(size_t a = 0; a < my_data.analyses().size(); ++a)
  {
    LSF_ptr<LSFBase> analysis_ptr = my_data.analyses()[a];
    assert(analysis_ptr != 0);

    analysis_specified = true;

    my_parser.parse(analysis_ptr);
    if(! my_parser.user_instructions().valid)
    {
      continue;    // Parser has already informed user.
    }
    
    // - Open output files.
    //
    ofstream  summary_file; 
    string    summary_file_name = my_parser.user_instructions().file_name_root + ".sum"; 
    summary_file.open(summary_file_name.c_str());
    if(! summary_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file '" << summary_file_name << "'.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    
    print_title(summary_file);
    
    ofstream  detail_file; 
    string    detail_file_name = my_parser.user_instructions().file_name_root + ".det"; 
    detail_file.open(detail_file_name.c_str());
    if(! detail_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file '" << detail_file_name << "'.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    } 
    
    print_title(detail_file);    
    
    ofstream  dump_file; 
    string    dump_file_name = my_parser.user_instructions().file_name_root + ".dmp"; 
    dump_file.open(dump_file_name.c_str());
    if(! dump_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file '" << dump_file_name << "'.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }     
    
    print_title(dump_file);
    
    analysis  my_analysis(my_data.get_ostreams(), detail_file, summary_file, dump_file, *my_mped, 
                          my_data.genome(), my_parser.user_instructions());
    
    try
    {
      my_analysis.analyze();
    }
    catch(const bad_alloc&)
    {
      my_data.get_ostreams().errors() << priority(critical) << "Not enough memory to process "
            << "analysis '" << my_analysis.title() << "'.  Skipping analysis ..." << endl;
      continue;
    }
  }

  if(! analysis_specified)
  {
    my_data.errors() << priority(error) << "Parameter file contains no analysis." << endl;
  }

  print_inf_banner(cout);

  return  EXIT_SUCCESS;
}

}
}


//============================================================================
// DECIPHER, the program
//============================================================================
//
int main(int argc, char* argv[])
{
  SAGE::DECIPHER::decipher decipher_inst(argc, argv);

  decipher_inst.main();
  
  exit(EXIT_SUCCESS);
}


