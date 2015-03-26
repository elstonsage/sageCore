//
//  File:     input.cpp
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston

#include <cassert>
#include "lvec/test_ivg_input.h"

using namespace std;
using namespace SAGE;

test_ivg_data::test_ivg_data(const string& program_name, bool debug)
               : APP::SAGE_Data(program_name, debug)
{}

test_ivg_data::~test_ivg_data()
{}

bool
test_ivg_data::input(int argc, char** argv)
{
  // Load parameter file & check for validity.
  //
  read_parameter_file(argv[1]);

  if( argc > 3 && argv[3] )
  {
    // Read in allele delimiter & marker locus file.
    //
    read_locus_description_file(argv[3]);
  }

  // Create RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  read_family_data_file(argv[2], false, true, true, false, true);

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read test_ivg_analysis block.
  //
  read_analysis();

  // Build genome.
  //
  build_genome();

  return true;
}

bool
test_ivg_data::read_analysis()
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
test_ivg_data::build_genome()
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
