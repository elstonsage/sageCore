//=======================================================================
//
//  File:	AppData.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2001 R. C. Elston
//=======================================================================

#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "ageon/AppData.h"
#include "ageon/Parser.h"

namespace SAGE {
namespace AO {


//=======================================================================
//
//  AppData(...) CONSTRUCTOR
//
//=======================================================================
AppData::AppData(const string& program_name, bool debug) :
	APP::SAGE_Simple_Data(program_name, debug)
{
   my_cmdline_rules.add_rule(
      APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE)
   );
   my_cmdline_rules.add_rule(
      APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE)
   );
   return;
}



//=======================================================================
//
//  process_input(...)
//
//=======================================================================
void
AppData::process_input(int argc, char** argv)
{
  parse_cmdline(argc, argv);
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  read_family_data_file(
    my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0],
    true,    // dump_trait
    false,    // dump_marker
    false,   // skip_traits
    true,   // skip_markers
    true     // dynamic_markers
  );

  evaluate_functions    ();
  read_analysis         ();
}


//=======================================================================
//
//  read_analysis()
//
//=======================================================================
bool
AppData::read_analysis()
{
  // 0. Set up local variables:

	Model  model;
	Parser parser(pedigrees(), my_output.messages(), errors());

  // 1. Clear out the analysis vector:

	my_analyses.clear();

  // 2. Iterate through all analysis blocks:

	for(LSFList::const_iterator i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
	{
	  if(!*i) continue;

	  if(toUpper((*i)->name()) == "AGEON")
	  {

  // 2.0.1. Parse the analysis block:

	    parser.parse_test_parameter_section(*i);

  // 2.0.2. Fetch the model from the parser:

	    model = parser.get_model();

  // 2.0.3. Make sure it's valid:

	    if(model.is_valid())
	    {
	      my_analyses.push_back(model);
	    }
	  }
	}

  // 3. Return true:

	return true;
}

} // End namespace AO
} // End namespace SAGE
