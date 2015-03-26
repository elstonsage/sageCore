//****************************************************************************
//* File:      fcor.cpp                                                      *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Oct 98 *
//*                    1.0 Modified as fcor object.               yjs Jan 01 *
//*                                                                          *
//* Notes:     This source file implements fcor app. derived from SAGEapp.   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "LSF/LSFinit.h"
#include "fcor/fcor.h"
#include "fcor/input.h"
#include "fcor/analysis.h"

using namespace std;

namespace SAGE {
namespace FCOR {

Fcor::Fcor(int argc, char** argv)
    : APP::SAGEapp(APP::APP_FCOR, true, argc, argv)
{
  LSFInit();
}

int Fcor::main()
{
  // Create input data.
  //
  FcorData fdata(name, debug());

  print_title(fdata.info());

  // Get the error stream and make it available locally

  cerrorstream& errors = fdata.errors();

  fdata.input(argc, argv);

  // No trait? No analysis!
  //
  if( !fdata.pedigrees().info().trait_count() )
  {
    errors << priority(error)
           << "No trait to be analyzed...  Exiting..." << endl;

    return EXIT_FAILURE;
  }

  // Start fcor analysis.
  //
  typedef boost::shared_ptr<ofstream> of_ptr;
  
  of_ptr  s_file;
  of_ptr  d_file;
  of_ptr  p_file;
  of_ptr  c_file;
  of_ptr  x_file;

  int f_anals = 0;

  for( size_t a = 0; a < fdata.analysis().size(); ++a )
  {
    LSF_ptr<LSFBase> i = fdata.analysis()[a];

    if( !i ) continue;

    std::cout << endl << "FCOR analysis........"
              << setw(3) << ++f_anals << endl << flush;

    FcorParser f_parser(fdata.pedigrees(), errors);

    f_parser.parse_fcor_analysis(i);

    if( !f_parser.get_trait_list().size() )
      continue;

    ostream* sum_output   = NULL;  
    ostream* det_output   = NULL;  
    ostream* pair_output  = NULL;  
    ostream* cov_output   = NULL;  
    ostream* xls_output   = NULL;  

    ofstream sum_file;
    ofstream det_file;
    ofstream pair_file;
    ofstream cov_file;
    ofstream xls_file;

    string   sum_filename;
    string   det_filename;
    string   pair_filename;
    string   cov_filename;
    string   xls_filename;

    AttrVal v = attr_value(i, "out");
    if( !v.has_value() || !v.String().size() )
      v = attr_value(i, "output");
    if(v.has_value() && v.String().size())
    {
      /*
      if( f_parser.get_pairset_type() != MAINTYPES )
      {
        sub_filename = v.String() + ".sub";
        sub_file.open( sub_filename.c_str() );

        if(sub_file)
          print_title(sub_file);
        sub_output = &sub_file;
      }

      if( f_parser.get_pairset_type() != SUBTYPES )
      {
        main_filename = v.String() + ".main";
        main_file.open( main_filename.c_str() );

        if(main_file)
          print_title(main_file);
        main_output = &main_file;
      }
      */

      sum_filename = v.String() + ".out";
      sum_file.open( sum_filename.c_str() );

      if(sum_file)
        print_title(sum_file);
      sum_output = &sum_file;

      if( f_parser.get_analysis_options().detailed_output ) 
      {
        det_filename = v.String() + ".det";
        det_file.open( det_filename.c_str() );

        if(det_file)
          print_title(det_file);
        det_output = &det_file;
      }

      if( f_parser.get_analysis_options().xls_output ) 
      {
        xls_filename = v.String() + ".alt";
        xls_file.open( xls_filename.c_str() );

        if(xls_file)
          print_title(xls_file);
        xls_output = &xls_file;
      }

      if( f_parser.get_analysis_options().pair_output ) 
      {
        pair_filename = v.String() + ".pair";
        pair_file.open( pair_filename.c_str() );

        if(pair_file)
          print_title(pair_file);
        pair_output = &pair_file;
      }

      if( f_parser.get_analysis_options().var_covs.size() ) 
      {
        cov_filename = v.String() + ".cov";
        cov_file.open( cov_filename.c_str() );

        if(cov_file)
          print_title(cov_file);
        cov_output = &cov_file;
      }
    }
    else
    {
      /*
      if( f_parser.get_pairset_type() != MAINTYPES )
      {
        sub_filename = "fcor.sub";
        if( s_file.get() == NULL )
        {
          s_file.reset( new ofstream(sub_filename.c_str()) );
          if(*s_file)
            print_title(*s_file);
        }
        sub_output = s_file.get();
      }

      if( f_parser.get_pairset_type() != SUBTYPES )
      {
        main_filename = "fcor.main";
        if( m_file.get() == NULL )
        {
          m_file.reset( new ofstream(main_filename.c_str()) );
          if(*m_file)
            print_title(*m_file);
        }
        main_output = m_file.get();
      }
      */

      sum_filename = "fcor.out";
      if( s_file.get() == NULL )
      {
        s_file.reset( new ofstream(sum_filename.c_str()) );
        if(*s_file)
          print_title(*s_file);
      }
      sum_output = s_file.get();

      if( f_parser.get_analysis_options().detailed_output ) 
      {
        det_filename = "fcor.det";
        if( d_file.get() == NULL )
        {
          d_file.reset( new ofstream(det_filename.c_str()) );
          if(*d_file)
            print_title(*d_file);
        }
        det_output = d_file.get();
      }

      if( f_parser.get_analysis_options().xls_output ) 
      {
        xls_filename = "fcor.alt";
        if( x_file.get() == NULL )
        {
          x_file.reset( new ofstream(xls_filename.c_str()) );
          if(*x_file)
            print_title(*x_file);
        }
        xls_output = x_file.get();
      }

      if( f_parser.get_analysis_options().pair_output ) 
      {
        pair_filename = "fcor.pair";
        if( p_file.get() == NULL )
        {
          p_file.reset( new ofstream(pair_filename.c_str()) );
          if(*p_file)
            print_title(*p_file);
        }
        pair_output = p_file.get();
      }

      if( f_parser.get_analysis_options().var_covs.size() ) 
      {
        cov_filename = "fcor.cov";
        if( c_file.get() == NULL )
        {
          c_file.reset( new ofstream(cov_filename.c_str()) );
          if(*c_file)
            print_title(*c_file);
        }
        cov_output = c_file.get();
      }
    }
    /*
    if( f_parser.get_pairset_type() != MAINTYPES && (!sub_output || !*sub_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << sub_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( f_parser.get_pairset_type() != SUBTYPES && (!main_output || !*main_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << main_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }
    */

    if( !sum_output || !*sum_output )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << sum_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( f_parser.get_analysis_options().detailed_output && (!det_output || !*det_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << det_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( f_parser.get_analysis_options().xls_output && (!xls_output || !*xls_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << xls_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( f_parser.get_analysis_options().pair_output && (!pair_output || !*pair_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << pair_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( f_parser.get_analysis_options().var_covs.size() && (!cov_output || !*cov_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << cov_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    FcorAnalysis fcor_analysis(f_parser, errors);

    fcor_analysis.run_analysis(*sum_output, *det_output, *pair_output,
                               *cov_output, *xls_output);
  }
 
  if( !f_anals )
  {
    ostream* sum_output   = NULL;  
    ostream* det_output   = NULL;  
    ostream* pair_output  = NULL;  
    ostream* cov_output   = NULL;  
    ostream* xls_output   = NULL;  

    ofstream sum_file;

    string   sum_filename = "fcor.out";

    if( s_file.get() == NULL )
    {
      s_file.reset( new ofstream(sum_filename.c_str()) );
      if(*s_file)
        print_title(*s_file);
    }
    sum_output = s_file.get();

    cout << endl << "  No analyses specified."
         << endl << "Performing FCOR default analysis..." << endl << flush;

    FcorParser f_parser(fdata.pedigrees(), errors);

    f_parser.build_default_fcor_analysis();

    FcorAnalysis fcor_analysis(f_parser, errors);

    fcor_analysis.run_analysis(*sum_output, *det_output, *pair_output,
                               *cov_output, *xls_output);
  }

  std::cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

} // end of namespace FCOR
} // end of namespace SAGE

int main(int argc, char* argv[])
{
  free(malloc(1));

  boost::shared_ptr<SAGE::FCOR::Fcor> fcor(new SAGE::FCOR::Fcor(argc, argv));

  if( !fcor.get() )
    exit(EXIT_FAILURE);

  fcor->main();

  exit(EXIT_SUCCESS);
}
