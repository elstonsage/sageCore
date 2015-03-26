//==========================================================================
//  File:      reltest.cpp
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                               Jul. 03
//
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "LSF/LSFinit.h"
#include "reltest/analysis.h"
#include "reltest/reltest.h"

using namespace std;

namespace SAGE
{

namespace RELTEST
{

reltest::reltest(int argc, char** argv)
       : APP::SAGEapp(APP::APP_RELTEST, true, argc, argv)
{
  LSFInit();
}

reltest::~reltest()
{}
    
int
reltest::main()
{
  // Create input data.
  //
  reltest_data  analy_data(name, debug());

  print_title(analy_data.info());

  // Get the error stream and make it available locally
  //
  cerrorstream& errors = analy_data.errors();

  analy_data.input(argc, argv);

  // Start reltest analysis
  //
  typedef boost::shared_ptr<ofstream> of_ptr;

  of_ptr  s_file;
  of_ptr  n_file;
  of_ptr  d_file;

  int r_anals = 0;

  for( size_t a = 0; a < analy_data.get_analysis().size(); ++a )
  {
    LSF_ptr<LSFBase> i = analy_data.get_analysis()[a];

    if( !*i ) continue;

    std::cout << endl << "RELTEST analysis......."
              << setw(3) << ++r_anals << endl << endl << flush;

    reltest_parser r_parser;

    r_parser.parse_reltest_analysis(analy_data, i, errors);

#if 0
    r_parser.view_parameter();
#endif

    ostream* sum_output = NULL;
    ostream* nuc_output = NULL;
    ostream* det_output = NULL;

    ofstream sum_file;
    ofstream nuc_file;
    ofstream det_file;

    string   sum_filename;
    string   nuc_filename;
    string   det_filename;

    AttrVal v = attr_value(i, "out");
    if( !v.has_value() || !v.String().size() )
      v = attr_value(i, "output");

    if( v.has_value() && v.String().size() )
    {
      string output = v.String();
      sum_filename = output + ".sum";
      sum_file.open( sum_filename.c_str() ); 

      if(sum_file)
        print_title(sum_file);
      sum_output = &sum_file;

      if( r_parser.generate_nucfam_output() )
      {
        nuc_filename = output + ".fam";
        nuc_file.open( nuc_filename.c_str() ); 
        if(nuc_file)
          print_title(nuc_file);
        nuc_output = &nuc_file;
      }

      if( r_parser.generate_detailed_output() )
      {
        det_filename = output + ".det";
        det_file.open( det_filename.c_str() ); 
        if(det_file)
          print_title(det_file);
        det_output = &det_file;
      }
    }
    else
    {
      sum_filename = "reltest.sum";
      if( s_file.get() == NULL )
      {
        s_file = of_ptr( new ofstream(sum_filename.c_str()) );
        if(*s_file)
          print_title(*s_file);
      }
      sum_output = s_file.get();

      if( r_parser.generate_nucfam_output() )
      {
        nuc_filename = "reltest.fam";
        if( n_file.get() == NULL )
        {
          n_file = of_ptr( new ofstream(nuc_filename.c_str()) );
          if(*n_file)
            print_title(*n_file);
        }
        nuc_output = n_file.get();
      }

      if( r_parser.generate_detailed_output() )
      {
        det_filename = "reltest.det";
        if( d_file.get() == NULL )
        {
          d_file = of_ptr( new ofstream(det_filename.c_str()) );
          if(*d_file)
            print_title(*d_file);
        }
        det_output = d_file.get();
      }
    }

    if( !sum_output || !*sum_output )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << sum_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( r_parser.generate_nucfam_output() && (!nuc_output || !*nuc_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << nuc_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( r_parser.generate_detailed_output() && (!det_output || !*det_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << det_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    reltest_analysis  rt_analysis(analy_data, r_parser, errors);

    rt_analysis.run_analysis(*sum_output, *nuc_output, *det_output);

#if 0
    rt_analysis.dump_pairs(cout);
#endif

  }

  if( !r_anals )
  {
    ostream* sum_output = NULL;
    ostream* nuc_output = NULL;
    ostream* det_output = NULL;

    ofstream sum_file;

    string   sum_filename = "reltest.sum";

    if( s_file.get() == NULL )
    {
      s_file = of_ptr( new ofstream(sum_filename.c_str()) );
      if(*s_file)
        print_title(*s_file);
    }
    sum_output = s_file.get();

    cout << endl << "  No analyses specified."
         << endl << "Performing RELTEST default analysis..." << endl << endl << flush;

    reltest_parser r_parser;

    r_parser.build_default_reltest_analysis(analy_data);

#if 0
    r_parser.view_parameter();
#endif

    reltest_analysis  rt_analysis(analy_data, r_parser, errors);

    rt_analysis.run_analysis(*sum_output, *nuc_output, *det_output);

#if 0
    rt_analysis.dump_pairs(cout);
#endif
  }

  cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

} // end of namespace RELTEST

} // end of namespace SAGE


int main(int argc, char** argv)
{
  free(malloc(1));

  boost::shared_ptr<SAGE::RELTEST::reltest> rel(new SAGE::RELTEST::reltest(argc, argv));

  if( !rel.get() )
    exit(EXIT_FAILURE);

  rel->main();

  exit(EXIT_SUCCESS);
}

