//=============================================================================
// File:    relpal.cpp
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                   yjs Mar 07
//
// Notes:   This source file implements relpal app. derived from SAGEapp.
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "LSF/LSFinit.h"
#include "relpal/input.h"
#include "relpal/analysis.h"
#include "relpal/relpal.h"

using namespace std;

namespace SAGE   {
namespace RELPAL {

Relpal::Relpal(int argc, char** argv)
      : APP::SAGEapp(APP::APP_RELPAL, true, argc,argv) 
{ 
  LSFInit();
}

int Relpal::main()
{
  // Create relpal input data.
  //
  relpal_data rdata(name, debug());

  print_title(rdata.info());

  rdata.input(argc, argv);

  cerrorstream& errors    = rdata.errors();
  ostream&      info_file = rdata.info();

  // Create filtered multipedigree.
  //
  FPED::Multipedigree fped(rdata.pedigrees());
  FPED::MPFilterer::add_multipedigree(fped, rdata.pedigrees());

  fped.construct();

  // Create RelativePairs.
  //
  relative_pairs rp;

  // Set the filtered multipedigree within relative_pairs.
  //
  rp.set_fped(fped);
  rp.prebuild();

  if( !rp.trait_count() )
  {
    errors << priority(critical)
           << "No traits to analyze in the pedigree data.  "
           << "Please have at least one valid 'trait' statement in the pedigree block. "
           << "Aborting." << endl;
    exit(EXIT_FAILURE);
  }

  size_t total_i_file = rdata.parsed_arguments().get_arguments(APP::IBD_FILE).size();
  bool multipoint_ibd = false;

  PALBASE::pair_analysis_ibd ibd(rp);

  if( total_i_file )
  {
    // Read IBD file.
    //
    cout << "Reading pairs............................." << flush;

    errors << priority(error);
    RefIBDReadFileStdIO ibd_reader(errors);
    ibd_reader.set_warn_invalid_pairs(true);

    for( size_t i_file = 0; i_file < total_i_file; ++i_file )
    {
      if( !ibd_reader.input( rdata.parsed_arguments().get_arguments(APP::IBD_FILE)[i_file], &ibd ) )
      {
        errors << priority(fatal)
               << "Error reading IBD file: "
               << rdata.parsed_arguments().get_arguments(APP::IBD_FILE)[i_file]
               << "...  Please use the correct ibd file from GENIBD without any modification. "
               << "Aborting." << endl;
        exit(EXIT_FAILURE);
      }
    }
    cout << "done." << endl;

    if( !rp.marker_count() )
    {
      errors << priority(critical)
             << "No markers to analyze... "
             << "Please make sure the ibd file has the information from the valid markers.  "
             << "Aborting." << endl;
      exit(EXIT_FAILURE);
    }

    // Sort pairs.
    //
    cout << "Sorting pairs............................." << flush;
    rp.build(false, true);
    cout << "done." << endl;

    if( !ibd.get_ibd_option().old_ibd_format )
    {
      if( toUpper(ibd.get_ibd_option().ibd_mode) == "MULTIPOINT" )
        multipoint_ibd = true;
      else
        multipoint_ibd = false;      
    }
    else
    {
      errors << priority(warning)
             << "There are no information whether the ibd values are singlepoint or multipoint.  "
             << "They will be processed as singlepoint ibd values."
             << endl;
    }
  }
  else
  {
    // Sort pairs.
    //
    cout << "Building pairs............................" << flush;
    rp.build_from_pedfile(rdata.pedigrees(), false, true);
    cout << "done." << endl;
  }

  if(!rp.valid())
  {
    errors << priority(fatal)
           << "Error building pair data structures!" << endl;
    exit(EXIT_FAILURE);
  }

  //rp.dump_pairs(cout);

  // Start relpal analysis.
  //
  typedef boost::shared_ptr<ofstream> of_ptr;

  of_ptr   s_file;
  of_ptr   d_file;
  of_ptr   e_file;

  int  r_anals = 0;
  bool run_sucess = false;

  for( size_t a = 0; a < rdata.analysis().size(); ++a )
  {
    LSF_ptr<LSFBase> i = rdata.analysis()[a];

    if( !i ) continue;

    cout << endl << "RELPAL analysis......."
         << setw(3) << ++r_anals << endl << flush;

    info_file << endl << "RELPAL analysis......."
              << setw(3) << r_anals << endl << endl << flush;

    relpal_parser r_parser(rp, errors);

    r_parser.parse_test_parameter_section(i);

    r_parser.dump_parser(info_file);

    if( !r_parser.get_traits().size() )
    {
      errors << priority(error)
             << "No Trait are specified in the analysis block!"
             << "  Please specify the trait(s) to be analyzed.  Skipping analysis..."
             << endl;
      
      continue;
    }

    if( !total_i_file && !r_parser.do_first_level_test() )
    {
      errors << priority(error)
             << "No IBD files are read!"
             << "  Please specify the ibd file(s) to do second level analysis.  Skipping analysis..."
             << endl;
      
      continue;
    }

    if(    r_parser.get_analysis_options().IBD_variance
        && r_parser.get_analysis_options().state_file_name.size() )
    {
      cout << "  Reading pair IBD states ................" << flush;

      errors << priority(error);
      RefIBDReadFileStdIO ibd_reader(errors);
      ibd_reader.set_warn_invalid_pairs(true);

      //PALBASE::pair_analysis_ibd ibd(rp);

      if( !ibd_reader.input_ibd_state(r_parser.get_analysis_options().state_file_name, &ibd) )
      {
        errors << priority(fatal)
               << "Error reading IBD state file: "
               << r_parser.get_analysis_options().state_file_name
               << "...  Please use the correct ibd state file from GENIBD without any modification. "
               << "Aborting." << endl;
        exit(EXIT_FAILURE);
      }

      cout << "done." << endl;
    }

    ostream* sum_output = NULL;
    ostream* det_output = NULL;
    ostream* exp_output = NULL;

    ofstream sum_file;
    ofstream det_file;
    ofstream exp_file;

    string   sum_filename;
    string   det_filename;
    string   exp_filename;

    AttrVal v = attr_value(i, "out");
    if( !v.has_value() || !v.String().size() )
      v = attr_value(i, "output");
    if( v.has_value() && v.String().size() )
    {
      sum_filename = v.String() + ".out";
      sum_file.open( sum_filename.c_str() );

      if( sum_file )
        print_title(sum_file);
      sum_output = &sum_file;

      if( r_parser.get_output_options().detailed_out )
      {
        det_filename = v.String() + ".det";
        det_file.open( det_filename.c_str() );

        if( det_file )
          print_title(det_file);
        det_output = &det_file;
      }

      if( r_parser.get_output_options().export_out )
      {
        exp_filename = v.String() + ".export";
        exp_file.open( exp_filename.c_str() );

        exp_output = &exp_file;
      }
    }
    else
    {
      sum_filename = "relpal.out";
      if( s_file.get() == NULL )
      {
        s_file.reset( new ofstream(sum_filename.c_str()) );
        if( *s_file )
          print_title(*s_file);
      }
      sum_output = s_file.get();

      if( r_parser.get_output_options().detailed_out )
      {
        det_filename = "relpal.det";
        if( d_file.get() == NULL )
        {
          d_file.reset( new ofstream(det_filename.c_str()) );
          if( *d_file )
            print_title(*d_file);
        }
        det_output = d_file.get();
      }

      if( r_parser.get_output_options().export_out )
      {
        exp_filename = "relpal.export";
        if( e_file.get() == NULL )
        {
          e_file.reset( new ofstream(exp_filename.c_str()) );
        }
        exp_output = e_file.get();
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

    if( r_parser.get_output_options().detailed_out && (!det_output || !*det_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open regression test output file: '" << det_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( r_parser.get_output_options().export_out && (!exp_output || !*exp_output) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open export output file: '" << exp_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    relpal_analysis r_analysis(errors);

    run_sucess = r_analysis.run_analysis(rp, r_parser, multipoint_ibd, *sum_output, *det_output, *exp_output);
  }

  if( !r_anals )
  {
    cout << endl << "  No analyses specified."
         << endl << "Please specify relpal analysis with at least one trait to perform an analysis."
         << endl << "Skipping..."
         << endl << endl << flush;
    
  }
  else if( !run_sucess )
    cout << endl << "Analysis failed!" << endl << endl;
  else
    cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

} // end of namespace RELPAL
} // end of namespace SAGE

int main(int argc, char* argv[])
{
  free(malloc(1));

  boost::shared_ptr<SAGE::RELPAL::Relpal> relpal(new SAGE::RELPAL::Relpal(argc, argv));

  if( !relpal.get() )
    exit(EXIT_FAILURE);

  relpal->main();

  exit(EXIT_SUCCESS);
}
