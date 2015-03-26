//============================================================================
// File:      input.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 2/12/4                                            
//                                                                          
// Notes:     implementation of class, decipher_data.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/input.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

decipher_data::decipher_data(const string& program_name, bool debug)
            : SAGE_Data(program_name, debug), my_true_marker_count(0)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));  
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE, APP::ArgumentRuleset::ZERO_OR_ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::GENOME_FILE, APP::ArgumentRuleset::ZERO_OR_ONE));  
}

size_t
decipher_data::true_marker_count() const
{
  return  my_true_marker_count;
}

void
decipher_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);
  
  if(my_parsed_arguments.argument_specified(APP::LOCUS_FILE))
  {
    read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);
  }  
  
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, true, false, false, true);
  
  filter_markers();
  filter_single_allelic_markers();

  if(my_parsed_arguments.argument_specified(APP::GENOME_FILE))
  {
    read_genome_description_file(my_parsed_arguments.get_arguments(APP::GENOME_FILE)[0], true, 2.0);
  }
  
  evaluate_functions();

  read_analysis();
}

bool
decipher_data::read_analysis()
{
  my_analyses.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    assert(*i != 0);

    if( toUpper((*i)->name() ) == "DECIPHER_ANALYSIS" || 
        toUpper((*i)->name() ) == "DECIPHER" )
      my_analyses.push_back(*i);
  }

  return true;
}

}
}
