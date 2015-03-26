//==========================================================================
//  File:      genibd.cpp
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                               Nov. 03
//
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "LSF/LSFinit.h"
#include "genibd/input.h"
#include "genibd/parser.h"
#include "genibd/analysis.h"
#include "genibd/genibd.h"

using namespace std;

namespace SAGE
{

namespace GENIBD
{

genibd::genibd(int argc, char** argv)
      : APP::SAGEapp(APP::APP_GENIBD, true, argc, argv)
{
  LSFInit();
}

genibd::~genibd()
{}
    
int
genibd::main()
{
  // Create input data.
  //
  genibd_data  analysis_data(name, debug());

  print_title(analysis_data.info());

  // Get the error stream and make it available locally
  //
  cerrorstream& errors = analysis_data.errors();

  analysis_data.input(argc, argv);

  // Start genibd analysis
  //
  int g_anals = 0;

  for( size_t a = 0; a < analysis_data.get_analysis().size(); ++a )
  {
    LSF_ptr<LSFBase> i = analysis_data.get_analysis()[a];

    if( !*i ) continue;

    cout << endl << "GENIBD analysis......."
         << right << setw(3) << ++g_anals << endl << endl << flush;

    genibd_parser g_parser(errors);

    g_parser.parse_test_parameter_section(i);

    genibd_parameters& analysis_param = g_parser.get_parameters();

    LSF_ptr<genome_description> genome
      = build_genome_description(analysis_data.pedigrees().info(),
                                 analysis_data.regions(),
                                 analysis_param.is_multipoint(),
                                 analysis_param.interval_distance(),
                                 errors);

    analysis_param.build_analysis_region(genome, g_parser.get_regions(), errors);

    if( !analysis_param.title().size() )
    {
      string title = "Analysis " + doub2str(g_anals);

      errors << priority(information)
             << "Analysis title unspecified.  Will use '"
             << title << "'." << endl;

      analysis_param.set_title(title);
    }

#if 0
    analysis_param.dump_parameters(cout);

    print_analysis_table(analysis_param, cout);
    print_analysis_table(analysis_param, analysis_data.info());
#endif

    genibd_analysis  g_analysis(analysis_data.pedigrees(),
                                analysis_param,
                                genome,
                                analysis_data.info(),
                                errors);

    g_analysis.run_analysis();

#if 0
    g_analysis.dump_pairs(cout);
#endif

  }

  if( !g_anals )
  {
    cout << endl << "  No analyses specified."
         << endl << "Performing GENIBD default analysis..." << endl << endl << flush;

    genibd_parser g_parser(errors);

    genibd_parameters& analysis_param = g_parser.get_parameters();

    LSF_ptr<genome_description> genome
      = build_genome_description(analysis_data.pedigrees().info(),
                                 analysis_data.regions(),
                                 analysis_param.is_multipoint(),
                                 analysis_param.interval_distance(),
                                 errors);

    analysis_param.build_analysis_region(genome, g_parser.get_regions(), errors);

#if 0
    g_parser.get_parameters().dump_parameters(cout);

    print_analysis_table(analysis_param, cout);
    print_analysis_table(analysis_param, analysis_data.info());
#endif

    genibd_analysis  g_analysis(analysis_data.pedigrees(),
                                analysis_param,
                                genome,
                                analysis_data.info(),
                                errors);

    g_analysis.run_analysis();

#if 0
    g_analysis.dump_pairs(cout);
#endif
  }

  cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

LSF_ptr<genome_description>
genibd::build_genome_description(const RefMPedInfo& minfo,
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
      vector<size_t> autos;
      vector<size_t> x_links;
      for( size_t i = 0; i < minfo.marker_count(); ++i )
      {
        if( minfo.marker_info(i).is_x_linked() )
          x_links.push_back(i);
        else
          autos.push_back(i);
      }

      if( autos.size() && x_links.size() )
      {
        errors << priority(warning)
               << "There are markers from both autosomal and X chromosomes. "
               << "Markers will be splitted into two regions "
               << "to perform singlepoint analysis." << endl;
      }

      genome_description* genome = new genome_description(minfo, m);

      if( autos.size() )
      {
        genome->add_region("default");

        for( size_t i = 0; i < autos.size(); ++i )
          genome->add_locus(autos[i], 0.0);
      }

      if( x_links.size() )
      {
        genome->add_region("default_x", true);

        for( size_t i = 0; i < x_links.size(); ++i )
          genome->add_locus(x_links[i], 0.0);
      }

      genome->set_mapping_function(genome_description::MapTypeShPtr(new RPED::Haldane()));
      genome->build();

      return genome;
    }
  }

  LSF_ptr<genome_description> genome = new RPED::LSFgenome_description(minfo, m, regions, multipoint);

  genome->set_scan_distance(distance);

  genome->build();

  return genome;
}

void
genibd::print_analysis_table(genibd_parameters& params, ostream& o)
{
  o << endl << "Analysis" << endl << "========" << endl << endl;

  params.dump_parameters(o);
}

} // end of namespace GENIBD

} // end of namespace SAGE


int main(int argc, char** argv)
{
  free(malloc(1));

  boost::shared_ptr<SAGE::GENIBD::genibd> gen(new SAGE::GENIBD::genibd(argc, argv));

  if( !gen.get() )
    exit(EXIT_FAILURE);

  gen->main();

  exit(EXIT_SUCCESS);
}

