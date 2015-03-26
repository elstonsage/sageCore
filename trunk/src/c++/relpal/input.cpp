//=============================================================================
// File:      input.cpp
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation                                 yjs Mar 01
//
// Notes:     This source file implements input.h
//
// Copyright (c) 2007 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/input.h"

using namespace std;

namespace SAGE   {
namespace RELPAL {

relpal_data::relpal_data(const string& program_name, bool debug)
           : APP::SAGE_Simple_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::IBD_FILE, APP::ArgumentRuleset::ZERO_OR_MORE));
}

relpal_data::~relpal_data()
{}

void
relpal_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, true, false, false, true);
            
  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read relpal_analysis block.
  //
  read_analysis();
}

bool
relpal_data::read_analysis()
{
  my_analysis.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;
    
    if(    toUpper((*i)->name()) == "RELPAL_ANALYSIS"
        || toUpper((*i)->name()) == "RELPAL" )
    {
      my_analysis.push_back(*i);
    }
  }

  return true;
}

} // end of namespace RELPAL
} // end of namespace SAGE
