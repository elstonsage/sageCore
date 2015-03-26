//****************************************************************************
//* File:      input.cpp                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Aug 01 *
//*            1.0 Added reading of marker locus file.            yjs Oct 01 *
//*                                                                          *
//* Notes:     This source file implements input.h                           *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/input.h"

using namespace std;

namespace SAGE {
namespace FCOR {

FcorData::FcorData(const string& program_name, bool debug)
        : SAGE_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE, APP::ArgumentRuleset::ZERO_OR_ONE));
}

FcorData::~FcorData()
{}

void
FcorData::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  bool dump_marker = false;

  if( my_parsed_arguments.argument_specified(APP::LOCUS_FILE) )
  {
    read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);

    // Print marker only when it is read.
    //
    dump_marker = true;
  }

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, dump_marker, false, false, true);

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read fcor_analysis block.
  // Start fcor analysis.
  //
  read_analysis();
}

bool
FcorData::read_analysis()
{
  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if(    toUpper((*i)->name() ) == "FCOR"
        || toUpper((*i)->name() ) == "FCOR2"
        || toUpper((*i)->name() ) == "FCOR_ANALYSIS"
        || toUpper((*i)->name() ) == "FCOR2_ANALYSIS" )
      my_analysis.push_back(*i);
  }

  return true;
}

} // end of namespace FCOR
} // end of namespace SAGE
