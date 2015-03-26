//****************************************************************************
//* File:      lodpal.cpp                                                    *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial Implementation.                yjs Mar 01 *
//*                                                                          *
//* Notes:     This source file implements lodpal app. derived from SAGEapp. *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "LSF/LSFinit.h"
#include "palbase/pal_ibd.h"
#include "palbase/pair_info_file.h"
#include "lodpal/lodpal.h"
#include "lodpal/lodpal_input.h"
#include "lodpal/lodpal_analysis.h"

using namespace std;

namespace SAGE   {
namespace LODPAL {

Lodpal::Lodpal(const string& program_name, int argc, char** argv)
      : APP::SAGEapp(APP::APP_LODPAL, true, argc, argv) 
{ 
  LSFInit();

  sage_cerr << prefix("%%LODPAL-%P: ");

  syms = new SymbolTable("LODPAL Symbol Table");

  syms->add( "wide_out",        "false");
  syms->add( "csv_out",         "false");
  syms->add( "write_ibd_file",  "false");
}

int Lodpal::main()
{
  // Create lodpal input data.
  //  
  lodpal_data ldata("lodpal", debug());

  print_title(ldata.info());

  cerrorstream& errors    = ldata.errors();
  ostream&      info_file = ldata.info();

  ldata.input(argc, argv);

  // Create filtered multipedigree.
  //
  FPED::Multipedigree fped(ldata.pedigrees());
  FPED::MPFilterer::add_multipedigree(fped, ldata.pedigrees());

  fped.construct();

  // Create RelativePairs & read pedigree file & check for validity.
  //
  RelativePairs rp;

  // Set the multipedigree with in RelativePairs.
  //
  rp.set_fped(fped);
  rp.prebuild();

  // Read IBD file.
  //
  cout << "Reading pairs............................." << flush;

  errors << priority(error);
  RefIBDReadFileStdIO ibd_reader(errors);
  ibd_reader.set_warn_invalid_pairs(true);

  PALBASE::pair_analysis_ibd ibd(rp);

  if( !ibd_reader.input( ldata.parsed_arguments().get_arguments(APP::IBD_FILE)[0], &ibd ) )
  {
    errors << priority(fatal) << "Error reading IBD file: "
           << ldata.parsed_arguments().get_arguments(APP::IBD_FILE)[0] << endl;
    exit(EXIT_FAILURE);
  }

  cout << "done." << endl;

  if( !rp.trait_count() )
  {
    errors << priority(critical) << "No traits to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  if( !rp.marker_count() )
  {
    errors << priority(critical) << "No markers to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  // Sort pairs.
  //
  cout << "Sorting pairs............................." << flush;
  rp.build(false, rp.x_linked_marker_exist());
  cout << "done." << endl;

  if(!rp.valid())
  {
    errors << priority(fatal) << "Error building pair data structures!" << endl;
    exit(EXIT_FAILURE);
  }
  
  // Start lodpal analysis.
  //
  typedef boost::shared_ptr<ofstream> of_ptr;

  of_ptr   lodpal_file;
  of_ptr   diagno_file;
  of_ptr   x_link_file;

  int lodpals  = 0;

  lodpal_parser  global_lodpal_parser(rp, syms, errors);

  for( size_t a = 0; a < ldata.analysis().size(); ++a )
  {
    LSF_ptr<LSFBase> i = ldata.analysis()[a];

    if( !i ) continue;

    if(    toUpper(i->name()) == "LODPAL_ANALYSIS"
        || toUpper(i->name()) == "LODPAL" )
    {
      global_lodpal_parser.parse_parameter(i);

      lodpal_parser local_lodpal_parser = global_lodpal_parser;

      cout << endl << "LODPAL analysis......."
           << setw(3) << ++lodpals << endl << endl << flush;

      rp.invalidate_pair_info();

      local_lodpal_parser.parse_test_parameter_section(i);

      if( !ibd.get_ibd_option().old_ibd_format )
      {
        if( toUpper(ibd.get_ibd_option().ibd_mode) == "MULTIPOINT" )
          local_lodpal_parser.set_multipoint(true);
        else
          local_lodpal_parser.set_multipoint(false);
      }
      // backward compatibility for old ibd file.
      else
      {
        if( has_attr(i, "multipoint") )
          local_lodpal_parser.set_multipoint(true);
        else if( has_attr(i, "singlepoint") )
          local_lodpal_parser.set_multipoint(false);
        else
        {
          cout << endl;
          sage_cerr << priority(error)
                    << "Please specify analysis type, multipoint or singlepoint."
                    << "  Skipping analysis..." << endl;
          continue;
        }
      }

      // Check if valid pair_info_file exist, then read the file.
      //
      rp.invalidate_pair_info();
      if( local_lodpal_parser.pair_info_file() )
        PALBASE::read_pair_info_file(local_lodpal_parser.get_pair_info_file(), rp, errors, info_file);

      ostream* lodpal_output = NULL;
      ostream* diagno_output = NULL;
      ostream* x_link_output = NULL;
      ofstream lodpal_outfile;
      ofstream diagno_outfile;
      ofstream x_link_outfile;
      string   lodpal_output_filename;
      string   diagno_output_filename;
      string   x_link_output_filename;

      AttrVal v = attr_value(i, "out");
      if( !v.has_value() || !v.String().size() )
        v = attr_value(i, "output");
      if(v.has_value() && v.String().size())
      {
        if(    (   local_lodpal_parser.parameters().marker_count() > 0
                && local_lodpal_parser.parameters().autosomal_marker_exist())
            || (   rp.marker_count() > 0
                && rp.autosomal_marker_exist()) )
        {
          lodpal_output_filename = v.String() + ".out";
          lodpal_outfile.open( lodpal_output_filename.c_str() );

          if(lodpal_outfile)
            print_title(lodpal_outfile);
          lodpal_output = &lodpal_outfile;
        }

        if(local_lodpal_parser.diagnostic())
        {
          diagno_output_filename = v.String() + ".lod";
          diagno_outfile.open( diagno_output_filename.c_str() );

          if(diagno_outfile)
            print_title(diagno_outfile);
          diagno_output = &diagno_outfile;
        }

        if(    (   local_lodpal_parser.parameters().marker_count() > 0
                && local_lodpal_parser.parameters().x_linked_marker_exist())
            || (   rp.marker_count() > 0
                && rp.x_linked_marker_exist()) )
        {
          x_link_output_filename = v.String() + ".xln";
          x_link_outfile.open( x_link_output_filename.c_str() );

          if(x_link_outfile)
            print_title(x_link_outfile);
          x_link_output = &x_link_outfile;
        }
      }
      else
      {
        if(    (   local_lodpal_parser.parameters().marker_count() > 0
                && local_lodpal_parser.parameters().autosomal_marker_exist())
            || (   rp.marker_count() > 0
                && rp.autosomal_marker_exist()) )
        {
          lodpal_output_filename = "lodpal.out";
          if( !lodpal_file )
          {
            lodpal_file.reset( new ofstream(lodpal_output_filename.c_str()) );
            if( *lodpal_file )
              print_title(*lodpal_file);
          }
          lodpal_output = lodpal_file.get();
        }

        if(local_lodpal_parser.diagnostic())
        {
          diagno_output_filename = "lodpal.lod";
          if( !diagno_file )
          {
            diagno_file.reset( new ofstream(diagno_output_filename.c_str()) );
            if( *diagno_file )
              print_title(*diagno_file);
          }
          diagno_output = diagno_file.get();
        }

        if(    (   local_lodpal_parser.parameters().marker_count() > 0
                && local_lodpal_parser.parameters().x_linked_marker_exist())
            || (   rp.marker_count() > 0
                && rp.x_linked_marker_exist()) )
        {
          x_link_output_filename = "lodpal.xln";
          if( !x_link_file )
          {
            x_link_file.reset( new ofstream(x_link_output_filename.c_str()) );
            if( *x_link_file )
              print_title(*x_link_file);
          }
          x_link_output = x_link_file.get();
        }
      }

      if(   (   (   local_lodpal_parser.parameters().marker_count() > 0
                 && local_lodpal_parser.parameters().autosomal_marker_exist())
             || (   rp.marker_count() > 0
                 && rp.autosomal_marker_exist()) )
          && (!lodpal_output || !*lodpal_output) )
      {
        cout << endl;
        sage_cerr << priority(error)
                  << "Cannot open output file: '" << lodpal_output_filename
                  << "'.  Skipping analysis..." << endl;
        continue;
      }

      if( local_lodpal_parser.diagnostic() && (!diagno_output || !*diagno_output) )
      {
        cout << endl;
        sage_cerr << priority(error)
                  << "Cannot open output file: '" << diagno_output_filename
                  << "'.  Skipping analysis..." << endl;
        continue;
      }

      if(   (   (   local_lodpal_parser.parameters().marker_count() > 0
                 && local_lodpal_parser.parameters().x_linked_marker_exist() )
             || (   rp.marker_count() > 0
                 && rp.x_linked_marker_exist()) )
          && (!x_link_output || !*x_link_output) )
      {
        cout << endl;
        sage_cerr << priority(error)
                  << "Cannot open output file: '" << x_link_output_filename
                  << "'.  Skipping analysis..." << endl;
        continue;
      }

      if( !local_lodpal_parser.multipoint() )
        cout << "Computing singlepoint lod scores.........." << flush;
      else
        cout << "Computing multipoint lod scores..........." << flush;

      lodpal_analysis lod_anal(errors);

      lod_anal.do_lodpal_test(rp, local_lodpal_parser, *lodpal_output, *diagno_output, *x_link_output);

      cout << "done." << endl;
    }
  }

  if( !lodpals )
  {
    ostream* lodpal_output = NULL;
    ostream* diagno_output = NULL;
    ostream* x_link_output = NULL;

    if( !lodpal_file )
    {
      if( rp.autosomal_marker_exist() )
      {
        lodpal_file.reset( new ofstream("lodpal.out") );
        if( *lodpal_file )
          print_title(*lodpal_file);
      }
      lodpal_output = lodpal_file.get();

      if( rp.x_linked_marker_exist() )
      {
        if( !x_link_file )
        {
          x_link_file.reset( new ofstream("lodpal.xln") );
          if( *x_link_file )
            print_title(*x_link_file);
        }
        x_link_output = x_link_file.get();
      }
    }
    cout << endl << "  No analyses specified."
         << endl << "Performing default lodpal analysis........" << flush;

    lodpal_analysis lod_anal(errors);
    lod_anal.do_lodpal_test(rp, global_lodpal_parser, *lodpal_output, *diagno_output, *x_link_output);

    cout << "done." << endl;
  }

  cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

} // end of namespace LODPAL
} // end of namespace SAGE

int main(int argc, char* argv[])
{
  free(malloc(1));

  boost::shared_ptr<SAGE::LODPAL::Lodpal> lodpal(new SAGE::LODPAL::Lodpal("lodpal", argc, argv));

  if( !lodpal.get() )
    exit(EXIT_FAILURE);

  lodpal->main();

  exit(EXIT_SUCCESS);
}
