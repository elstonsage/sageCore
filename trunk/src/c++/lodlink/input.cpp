//============================================================================
// File:      input.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 9/20/2                                            
//                                                                          
// Notes:     implementation of class, lodlink_data.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/input.h"

using namespace std;

namespace SAGE
{

namespace LODLINK
{

lodlink_data::lodlink_data(const string& program_name, bool debug)
            : APP::SAGE_Data(program_name, debug), my_true_marker_count(0)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));  
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::TRAIT_MODEL_FILE, APP::ArgumentRuleset::ZERO_OR_ONE));  
}

size_t
lodlink_data::true_marker_count() const
{
  return  my_true_marker_count;
}

// - Read input files whose names are specified on the command line.
//
void
lodlink_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);
  read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);

  my_true_marker_count = my_pedigrees.info().marker_count();
  
  if(my_parsed_arguments.argument_specified(APP::TRAIT_MODEL_FILE))
  {
    read_locus_description_file(my_parsed_arguments.get_arguments(APP::TRAIT_MODEL_FILE)[0]);
  }
  
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], false, true, true, false, false);
  filter_x_linked_markers();
  filter_markers();
  
  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  read_analysis();
}

bool
lodlink_data::read_analysis()
{
  my_analyses.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    assert(*i != 0);

    if( toUpper((*i)->name() ) == "LODLINK_ANALYSIS" || 
        toUpper((*i)->name() ) == "LODLINK" )
      my_analyses.push_back(*i);
  }

  return true;
}

}
}
