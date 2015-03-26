//
//  File:     input.cpp
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston

#include <cassert>
#include "markerinfo/input.h"

using namespace std;
using namespace SAGE;

markerinfo_data::markerinfo_data(const string& program_name, bool debug)
               : APP::SAGE_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,  APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE,     APP::ArgumentRuleset::ZERO_OR_ONE));
}

markerinfo_data::~markerinfo_data()
{}

bool
markerinfo_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);
  
  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  if( my_parsed_arguments.argument_specified(APP::LOCUS_FILE) )
  {
    // Read in allele delimiter & marker locus file.
    //
    read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);
  }

  // Create RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0],
                        false, true, false, false, true);

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read markerinfo_analysis block.
  //
  read_analysis();

  // Build genome.
  //
  build_genome();

  return true;
}

bool
markerinfo_data::read_analysis()
{
  my_analysis.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if(    toUpper((*i)->name() ) == "MARKERINFO"
        || toUpper((*i)->name() ) == "MARKERINFO_ANALYSIS" )
      my_analysis.push_back(*i);
  }

  return true;
}

void
markerinfo_data::build_genome()
{
  my_genome = new RPED::genome_description(my_pedigrees.info(), errors());

  my_genome->add_region("REGION");

  for(size_t i = 0; i < markers().size(); ++i)
    my_genome->add_locus(i, 0.0);

  // Doesn't matter, but needs to be there for correct building
  my_genome->set_mapping_function(RPED::genome_description::MapTypeShPtr(new RPED::Haldane()));

  my_genome->build();

  my_genome->freeze();
}
