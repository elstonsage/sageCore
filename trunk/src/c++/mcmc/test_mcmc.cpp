//==========================================================================
//  File:      test_mcmc.cpp
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                               May. 04
//
//  Notes:     Tests the components of mcmc.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "LSF/LSFinit.h"
#include "mcmc/test_mcmc.h"
#include "mcmc/test_mcmc_input.h"

using namespace std;

namespace SAGE
{

namespace MCMC
{

test_mcmc::test_mcmc(int argc, char** argv)
           : APP::SAGEapp(APP::TEST_PROGRAM, false, argc, argv)
{
  if( arg_count < 3 )
  {
    print_help(cerr);
    exit(EXIT_FAILURE);
  }

  LSFInit();
}

test_mcmc::~test_mcmc()
{}
    
void  
test_mcmc::print_help(std::ostream& o)
{
  o << "usage: " << argv[0] << " <parameters> <pedigree> <locus> [map]" 
    << endl << endl
    << "Command line parameters:"                << endl
    << "  parameters   - Parameter File"         << endl
    << "  pedigree     - Pedigree Data File"     << endl
    << "  locus        - Locus Description File" << endl
    << "  map          - Genome Map File (optional for single point analysis)"
    << endl << endl;
}

int
test_mcmc::main()
{
  // Create input data.
  //
  test_mcmc_data  analysis_data(name, debug());

  print_title(analysis_data.info());

  // Get the error stream and make it available locally
  //
  cerrorstream& errors = analysis_data.errors();

  analysis_data.input(argc, argv);

  // Start test_mcmc analysis
  //
  cout << endl << "  No analyses specified."
       << endl << "Performing MCMC_ENGINE default analysis..." << endl << endl << flush;

  test_mcmc_parser g_parser(errors);

  test_mcmc_parameters& analysis_param = g_parser.get_parameters();
  
//  print_analysis_table(analysis_param, cout);

  LSF_ptr<RPED::genome_description> genome
    = build_genome_description(analysis_data.pedigrees().info(),
                               analysis_data.regions(),
                               analysis_param.is_multipoint(),
                               2.0,
                               errors);

  analysis_param.build_analysis_region(genome, g_parser.get_regions(), errors);

  test_mcmc_analysis  g_analysis(analysis_data.pedigrees(),
                                 analysis_param,
                                 genome,
                                 analysis_data.info(),
                                 errors);

  g_analysis.run_analysis();

  cout << endl << "Analysis complete!" << endl << endl;

  return EXIT_SUCCESS;
}

LSF_ptr<RPED::genome_description>
test_mcmc::build_genome_description(const RPED::RefMPedInfo& minfo,
                                    LSFBase*           regions,
                                    bool               multipoint,
                                    double             distance,
                                    cerrorstream&      errors)
{
  cerrormultistream m;

  m.insert(errors);

  if( !regions || !regions->List() )
  {
    if( multipoint )
    {
      errors << priority(fatal)
             << "A Genome map must be provided when "
             << "performing multipoint analysis." << endl;

      exit(3);
    }
    else
    {
      RPED::genome_description* genome = new RPED::genome_description(minfo, m);
      genome->add_region("default");

      for( int i = 0; i < (int)minfo.marker_count(); ++i )
        genome->add_locus(i, 0.0);

      genome->set_mapping_function(RPED::genome_description::MapTypeShPtr(new RPED::Haldane()));
      genome->build();

      return genome;
    }
  }

  LSF_ptr<RPED::genome_description> genome = new RPED::LSFgenome_description(minfo, m, regions, multipoint);

  genome->set_scan_distance(distance);

  genome->build();

  return genome;
}

void
test_mcmc::print_analysis_table(test_mcmc_parameters& params, ostream& o)
{
  o << endl << "Analysis" << endl << "========" << endl << endl;

  params.dump_parameters(o);
}

} // end of namespace MCMC

} // end of namespace SAGE


int main(int argc, char** argv)
{
  free(malloc(1));

  boost::shared_ptr<SAGE::MCMC::test_mcmc> gen(new SAGE::MCMC::test_mcmc(argc, argv));

  if( !gen.get() )
    exit(EXIT_FAILURE);

  gen->main();

  exit(EXIT_SUCCESS);
}

