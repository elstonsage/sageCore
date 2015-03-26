#include "data/ArgumentRuleset.h"
#include "LSF/LSFinit.h"
#include "freq/freq.h"

namespace SAGE {
namespace FREQ {

//=========================
//
//  CONSTRUCTOR
//
//=========================
freq::freq(int argc, char** argv) : APP::SAGEapp(APP::APP_FREQ, true, argc, argv)
{
  // Initialize LSF:
  LSFInit();

  // Set error output prefix:
  sage_cerr << prefix("%%FREQ-%P: ");
}

freq_data::freq_data(const std::string& program_name, bool debug)
      : APP::SAGE_Simple_Data(program_name, debug)
{
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,  APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE,     APP::ArgumentRuleset::ZERO_OR_ONE));
}

void
freq_data::input(int argc, char** argv)
{
  parse_cmdline(argc, argv);
    
  read_parameter_file(my_parsed_arguments.get_arguments(APP::PARAMETER_FILE)[0]);

  if(my_parsed_arguments.argument_specified(APP::LOCUS_FILE))
  {
    read_locus_description_file(my_parsed_arguments.get_arguments(APP::LOCUS_FILE)[0]);
  }

  read_family_data_file(my_parsed_arguments.get_arguments(APP::PEDIGREE_FILE)[0], true, true, false, false, true);

  evaluate_functions();    
}

//===========================================
//
//  main()
//
//===========================================
int freq::main()
{
  // Create freq_data instance and read in input files:
  //
  freq_data data(name, debug());
  print_title(data.info());
  data.input(argc, argv);

  // Construct filtered multipedigree:
  //
  FPED::FilteredMultipedigree fped(data.pedigrees());
 
  FPED::MPFilterer::add_multipedigree_filtered_by_subpedigrees(fped, data.pedigrees(), boost::bind(&MPED::mp_utilities::no_loops, _1));
  FPED::MPFilterer::add_multipedigree_filtered_by_unconnecteds(fped, data.pedigrees(), FPED::always_keep());

  fped.construct();   

  // Parse parameter file into a list of Configuration's:

  std::vector<Configuration> configs;

  for(LSFList::iterator i = data.parameters()->List()->begin(); i != data.parameters()->List()->end(); ++i)
  {
    if(*i && toUpper((*i)->name()) == "FREQ" && (*i)->List())
    {
      configs.push_back(Parser::parseFreqBlock(configs.size(), *i, fped));
    }
  }

  if(!configs.size())
  {
    configs.push_back(Parser::createDefaultConfig(fped));
  }

  // Execute Configuration's:

  for(std::vector<Configuration>::const_iterator config_itr = configs.begin(); config_itr != configs.end(); ++config_itr)
  {
    Output::generateOutput(Estimator::runAnalysis(*config_itr, fped), fped);
  }
    
  print_inf_banner(cout);

  // Return success:
  
  return EXIT_SUCCESS;
}

} // End namespace FREQ
} // End namespace SAGE

//================================
//          main(...)
//================================
int main(int argc, char* argv[])
{
    return SAGE::FREQ::freq(argc, argv).main();
}

