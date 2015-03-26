//****************************************************************************
//* File:      lodpal_input.cpp                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Aug 01 *
//*                                                                          *
//* Notes:     This source file implements lodpal_input.h                    *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_input.h"

namespace SAGE   {
namespace LODPAL {

lodpal_data::lodpal_data(const string& program_name, bool debug)
           : APP::SAGE_Simple_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::IBD_FILE, APP::ArgumentRuleset::ONE));
}

lodpal_data::~lodpal_data()
{}

void
lodpal_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, false, false, true, false);

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read lodpal_analysis block.
  // Start lodpal analysis.
  //
  read_analysis();
}

bool
lodpal_data::read_analysis()
{
  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if(    toUpper((*i)->name() ) == "LODPAL_ANALYSIS"
        || toUpper((*i)->name() ) == "LODPAL"
        || toUpper((*i)->name() ) == "WIDE_OUT" )
      my_analysis.push_back(*i);
  }

  return true;
}

} // end of namespace LODPAL
} // end of namespace SAGE
