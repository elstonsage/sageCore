//============================================================================
// File:      data.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   5/30/02 - created.                         djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              data 
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "mlod/data.h"

//using namespace std;

namespace SAGE
{

namespace MLOD
{


//============================================================================
// IMPLEMENTATION:  data
//============================================================================
//
// - Read and store the data needed for the mlod application.
//
void
Data::input(int argc, char** const argv)
{
  parse_cmdline(argc, argv);
  
  //lint --e{534} Ignoring returns
  // Deal with parameter file and symbol table (derived from parameter file)

  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);
  
  // Read in the locus models
  
  read_locus_models(argv);
  
  copy_locus_models_to_multipedigree();

  //read pedigree data file

  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], false, true, true, false, false);

  // - Genome description file.
  //
  read_genome_description_file(my_parsed_arguments.get_arguments(APP::GENOME_FILE)[0], true, 2.0);
  
//  build_trait_genome();

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // - Find and store pointers to mlod analysis blocks.
  //
  read_analysis();
}

bool
Data::read_analysis()
{
  // Create a Parser object
  Parser  mlod_parser(*this, messages(), errors());

  // Iterate through the parameter file parameters, looking for
  // MLOC analysis blocks.
  LSFList::const_iterator i;
  for(i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
  {
    if(! *i) 
    {
      continue;
    }

    if(toUpper((*i)->name()) == "MLOD"          || 
       toUpper((*i)->name()) == "MLOD_ANALYSIS"   )
    {
      mlod_parser.parse(*i);

      bool valid_analysis = mlod_parser.parameters().is_valid();
      
      if(valid_analysis)
        my_analyses.push_back(mlod_parser.parameters());
    }
  }

  return true;
}

void Data::read_locus_models(char* const _argv[])
{
  // Read the marker loci
  cout << "Reading Marker Locus Description File....." << flush;

  read_locus_description_file_to_map(my_marker_loci, my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);

  cout << "done." << endl;

  // Read the trait loci
  cout << "Reading Trait Model Description File......" << flush;

  read_locus_description_file_to_map(my_trait_loci, my_parsed_arguments.get_arguments(APP::TRAIT_MODEL_FILE)[0]);

  cout << "done." << endl;
}


void
Data::read_locus_description_file_to_map
    (MLOCUS::inheritance_model_map& imap, const string& fname)
{
  MLOCUS::InheritanceModelFile marker_reader(errors());

  if( !marker_reader.input(imap, fname,
                           my_pedigrees.info().get_pheno_reader_info()) )
  {
    cout << endl;
    errors() << priority(fatal)
             << "Error reading model descriptions from file '" << fname
             << "'."
             << endl;
    
    exit(EXIT_FAILURE);
  }
}

void Data::copy_locus_models_to_multipedigree()
{
  my_pedigrees.info().markers().insert(my_marker_loci.index_begin(),
                                       my_marker_loci.index_end());
  my_pedigrees.info().markers().insert(my_trait_loci.index_begin(),
                                       my_trait_loci.index_end());
}

}
}

