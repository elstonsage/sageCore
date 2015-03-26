//==========================================================================
//  File:     input.cpp
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.                                Nov. 03
//
//  Copyright (c) 2003 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/input.h"

namespace SAGE
{

namespace GENIBD
{

genibd_data::genibd_data(const string& program_name, bool debug)
           : APP::SAGE_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,  APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE,     APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::GENOME_FILE,     APP::ArgumentRuleset::ZERO_OR_ONE));
}

genibd_data::~genibd_data()
{}

bool
genibd_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);
  
  //read parameter file
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Read in allele delimiter & marker locus file.
  read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);

  //read pedigree data file
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], false, true, true, false, false);

  //read genome file
  //
  if( my_parsed_arguments.argument_specified(APP::GENOME_FILE) )
    read_genome_description_file(my_parsed_arguments.get_arguments(APP::GENOME_FILE)[0]);

  //read genibd analysis
  //
  read_analysis();

  return true;
}

bool
genibd_data::read_analysis()
{
  my_analysis.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if( toUpper((*i)->name() ) == "GENIBD_ANALYSIS" ||
        toUpper((*i)->name() ) == "GENIBD" ||
        toUpper((*i)->name() ) == "IBD_ANALYSIS" )
      my_analysis.push_back(*i);
  }

  return true;
}

} // end of namespace GENIBD

} // end of namespace SAGE
