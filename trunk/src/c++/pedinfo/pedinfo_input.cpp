#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "pedinfo/pedinfo_input.h"
#include "data/ArgumentRuleset.h"

using namespace std;
using namespace SAGE;

pedinfo_data::pedinfo_data(const string& program_name, bool debug)
            : APP::SAGE_Simple_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE, APP::ArgumentRuleset::ONE));  
}

pedinfo_data::~pedinfo_data()
{}

void
pedinfo_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);

  // Load parameter file & check for validity.
  //
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  // Create RPED::RefMultiPedigree & read pedigree file & check for validity.
  // Check for the errors in pedigrees.
  // Sort pedigrees.
  //
  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, false, false, false, true);

  // - Create traits specified by function blocks.  -djb 5/2/01
  //
  evaluate_functions();

  // Read pedinfo_analysis block.
  //
  read_analysis();
}

bool
pedinfo_data::read_analysis()
{
  my_analyses.resize(0);

  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if( toUpper((*i)->name() ) == "PEDINFO_ANALYSIS" || 
        toUpper((*i)->name() ) == "PEDINFO" )
      my_analyses.push_back(*i);
  }

  return true;
}
