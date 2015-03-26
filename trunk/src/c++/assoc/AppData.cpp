//========================================================================
//
//  File:	AppData.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2001 R. C. Elston
//========================================================================

#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "assoc/AppData.h"

namespace SAGE  {
namespace ASSOC {

//========================================================================
//
//  AppData(...) CONSTRUCTOR
//
//========================================================================
AppData::AppData(const string& program_name, bool debug_on) :
	APP::SAGE_Simple_Data(program_name, debug_on)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,  APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE,     APP::ArgumentRuleset::ZERO_OR_ONE));
}

//========================================================================
//
//  process_input(...)
//
//========================================================================
void
AppData::process_input(int argc, char** argv)
{
  parse_cmdline(argc, argv);
  
  //lint --e{818) Pointer 'argv' could be const
  //lint --e{534} Ignoring return value

  read_parameter_file   ( my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0] );

  if( my_parsed_arguments.argument_specified(APP::LOCUS_FILE) )
    //lint -e(534) Ignoring return value
    read_locus_description_file( my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0] );

  read_family_data_file ( my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0] , true, false, false, false, true);
  evaluate_functions    ();
  read_analysis         ();
}

//========================================================================
//
//  read_analysis()
//
//========================================================================
bool
AppData::read_analysis()
{
  my_analyses.clear();

  for(LSFList::const_iterator i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
  {
    if(!*i) continue;

    if(toUpper((*i)->name()) == "ASSOC_ANALYSIS" || toUpper((*i)->name()) == "ASSOC")
    {
      try
      {
        my_analyses.push_back(Parser::parseAssocBlock(my_analyses.size(), *i, pedigrees(), my_output.info(), errors()));
      }
      catch(const std::exception &)
      {
        // Do nothing - analysis invalid
      }
    }
  }

  return true;
}

} // End namespace ASSOC
} // End namespace SAGE
