//****************************************************************************
//* File:      sibpal_input.cpp                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Aug 01 *
//*                                                                          *
//* Notes:     This source file implements sibpal_input.h                    *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "sibpal/input.h"

using namespace std;

namespace SAGE   {
namespace SIBPAL {

sibpal_data::sibpal_data(const string& program_name, bool debug)
           : APP::SAGE_Simple_Data(program_name, debug)
{
  my_dump_pairs = false;

  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::IBD_FILE, APP::ArgumentRuleset::ONE_OR_MORE));
}

sibpal_data::~sibpal_data()
{}

void
sibpal_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.  Check pedigrees.
  //
  my_dump_pairs = read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, true, false, false, true);
            
  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read sibpal_analysis block.
  // Start sibpal analysis.
  //
  read_analysis();
}

bool
sibpal_data::read_analysis()
{
  my_analysis.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;
    
    LSFBase* param = *i;

    string name = toUpper((*i)->name());

    if( name == "TRAIT_REGRESSION" )
    {
      string output_file_name = "";

      AttrVal v = attr_value(param, "out");
      if( !v.has_value() || !v.String().size() )
        v = attr_value(param, "output");
      if( v.has_value() && v.String().size() )
      {
        output_file_name = v.String();
      }

      my_treg_analysis.push_back(make_pair(*i, output_file_name));
    }
    else if(    name == "MEAN_TEST"         || name == "MEANS_TEST"
             || name == "MARKER_REGRESSION" || name == "MEAN_REGRESSION" )
    {
      string output_file_name = "";

      AttrVal v = attr_value(param, "out");
      if( !v.has_value() || !v.String().size() )
        v = attr_value(param, "output");
      if( v.has_value() && v.String().size() )
      {
        output_file_name = v.String();
      }

      my_mean_analysis.push_back(make_pair(*i, output_file_name));
    }
    else if( name == "SIBPAL_ANALYSIS" || name == "SIBPAL" )
    {
      string output_file_name = "";

      AttrVal v = attr_value(param, "out");
      if( !v.has_value() || !v.String().size() )
        v = attr_value(param, "output");
      if( v.has_value() && v.String().size() )
      {
        output_file_name = v.String();
      }

      my_analysis.push_back(make_pair(*i, output_file_name));
    }
    else if( name == "TRAIT_REGRESSION_DEFAULT" )
    {
      string output_file_name = "";
      my_analysis.push_back(make_pair(*i, output_file_name));
    }
  }

  return true;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
