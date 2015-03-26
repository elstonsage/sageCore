#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "segreg/segreg_input.h"
#include "segreg/parser.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{


segreg_data::segreg_data(const string& program_name, bool debug)
  : APP::SAGE_Simple_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,  APP::ArgumentRuleset::ONE));
}

segreg_data::~segreg_data()
{}

void
segreg_data::input(int argc, char** argv)
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

  read_analysis();
}

bool
segreg_data::read_analysis()
{
  my_analyses.resize(0);

  locus_indic_vec.clear(); // new member of segreg_data class added by JA

  SEGREG::parser p(&pedigrees(), my_output.messages(), errors());
  
  LSFList::const_iterator i;
  for( i = my_params->List()->begin(); i != my_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if(    toUpper((*i)->name() ) == "SEGREG_ANALYSIS"
        || toUpper((*i)->name() ) == "SEGREG" )
    {
      p.zero_poly_loci = false;
      p.parse_test_parameter_section(*i);

      if(p.get_model().get_model_class() != model_INVALID)
      {
         my_analyses.push_back(p.get_model()); // original code

           if (p.zero_poly_loci){locus_indic_vec.push_back(1);}
            else {locus_indic_vec.push_back(0);}
           like_cutoff = p.like_cutoff; 
       }
      }
    }
  return true;
}

}
}
