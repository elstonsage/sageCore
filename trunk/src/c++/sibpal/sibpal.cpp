#include "LSF/LSFinit.h"
#include "sibpal/input.h"
#include "sibpal/parser.h"
#include "sibpal/analysis.h"
#include "sibpal/sibpal.h"

using namespace std;

namespace SAGE   {
namespace SIBPAL {

Sibpal::Sibpal(int argc, char** argv)
      : APP::SAGEapp(APP::APP_SIBPAL, true, argc,argv) 
{ 
  LSFInit();

  sage_cerr << prefix("%%SIBPAL-%P: ");

  syms = new SymbolTable("Sibpal Symbol Table");

  syms->add( "wide_out", "false" );
  syms->add( "csv_out",  "false" );
}

int Sibpal::main()
{
  // Create sibpal input data.
  //
  sibpal_data sdata("sibpal", debug());

  print_title(sdata.info());

  sdata.input(argc, argv);

  cerrorstream& errors    = sdata.errors();
  ostream&      info_file = sdata.info();

  // Create filtered multipedigree.
  //
  FPED::Multipedigree fped(sdata.pedigrees());
  FPED::MPFilterer::add_multipedigree(fped, sdata.pedigrees());

  fped.construct();

  // Create RelativePairs.
  //
  relative_pairs rp;

  // Set the multipedigree with in RelativePairs.
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

  size_t total_i_file = sdata.parsed_arguments().get_arguments(APP::IBD_FILE).size();

  if( total_i_file )
  {
    // Read IBD file.
    //
    cout << "Reading pairs............................." << flush;

    errors << priority(error);
    RefIBDReadFileStdIO ibd_reader(errors);
    ibd_reader.set_warn_invalid_pairs(true);

    PALBASE::pair_analysis_ibd ibd(rp);

    for( size_t i_file = 0; i_file < total_i_file; ++i_file )
    {
      if( !ibd_reader.input( sdata.parsed_arguments().get_arguments(APP::IBD_FILE)[i_file], &ibd ) )
      {
        errors << priority(fatal)
               << "Error reading IBD file: "
               << sdata.parsed_arguments().get_arguments(APP::IBD_FILE)[i_file]
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
    rp.build(true, rp.x_linked_marker_exist());
    cout << "done." << endl;
  }
  else
  {
    // Sort pairs.
    //
    cout << "Building pairs............................" << flush;
    rp.build_from_pedfile(sdata.pedigrees(), true, rp.x_linked_marker_exist());
    cout << "done." << endl;
  }

  if(!rp.valid())
  {
    errors << priority(fatal)
           << "Error building pair data structures!" << endl;
    exit(EXIT_FAILURE);
  }

  // Start sibpal analysis.
  //
  typedef boost::shared_ptr<ofstream> of_ptr;

  of_ptr   mean_file;
  of_ptr   mean_csv_file;

  of_ptr   reg_file;
  of_ptr   reg_csv_file;
  of_ptr   reg_det_file;

  int sibpals = 0;

  regression_type default_reg_type = SINGLE_MARKER;

  regression_parser global_trait_parser(rp, syms, errors);
  meantest_parser   global_mean_parser(rp, syms, errors);

  for( size_t anal = 0; anal < sdata.analysis().size(); ++anal )
  {
    LSF_ptr<LSFBase> params = sdata.analysis()[anal].first;
    string output_name_root = sdata.analysis()[anal].second;

    if( !params ) continue;

    AttrVal a = attr_value(params,"TRAIT_REGRESSION_DEFAULT", 0);
    if( a.has_value() )
    {
      if( toUpper(a.String()) == "MULTIPLE" || toUpper(a.String()) == "MULTIPLE_MARKER" )
        default_reg_type = MULTIPLE_MARKER;
      else if( toUpper(a.String()) == "ZERO" || toUpper(a.String()) == "ZERO_MARKER" )
        default_reg_type = ZERO_MARKER;
    }

    if(    toUpper((params)->name()) == "SIBPAL"
        || toUpper((params)->name()) == "SIBPAL_ANALYSIS" )
    {
      size_t inside_tregs = 0;
      size_t inside_means = 0;

      ostream* treg_output = NULL;
      ostream* treg_csvput = NULL;
      ostream* treg_detput = NULL;

      ofstream treg_outfile;
      ofstream treg_csvfile;
      ofstream treg_detfile;

      string   treg_output_filename;
      string   treg_csvput_filename;
      string   treg_detput_filename;

      ostream* mean_output = NULL;
      ostream* mean_csvput = NULL;

      ofstream mean_outfile;
      ofstream mean_csvfile;

      string   mean_output_filename;
      string   mean_csvput_filename;

      LSFList::const_iterator i;
      for( i = params->List()->begin(); i != params->List()->end(); ++i )
      {
        if( !*i ) continue;

        LSFBase* param = *i;

        string name = toUpper(param->name());

        if( name == "TRAIT_REGRESSION" )
        {
          regression_parser local_trait_parser = global_trait_parser;

          regression_type local_reg_type = default_reg_type;
          
          if( has_attr(param, "single") || has_attr(param, "single_marker") )
            local_reg_type = SINGLE_MARKER;
          else if( has_attr(param, "multiple") || has_attr(param, "multiple_marker") )
            local_reg_type = MULTIPLE_MARKER;
          else if( has_attr(param, "zero") || has_attr(param, "zero_marker") )
            local_reg_type = ZERO_MARKER;

          cout << endl << "SIBPAL analysis......."
               << setw(3) << ++sibpals  << endl << endl << flush;

          local_trait_parser.set_regression_type(local_reg_type);
          local_trait_parser.parse_test_parameter_section(param);

          // Check if valid pair_info_file exist, then read the file.
          //
          if( local_trait_parser.get_pair_info_file().size() )
            PALBASE::read_pair_info_file(local_trait_parser.get_pair_info_file(), rp, errors, info_file);

          //local_trait_parser.dump_parser(cout);

          if( output_name_root.size() )
          {
            treg_output_filename = output_name_root + ".treg";
            treg_detput_filename = output_name_root + ".treg_det";

            if( !inside_tregs )
            {
              treg_outfile.open( treg_output_filename.c_str() );
              treg_detfile.open( treg_detput_filename.c_str() );

              if( treg_outfile )
                print_title(treg_outfile);

              if( treg_detfile )
                print_title(treg_detfile);
            }
            else
            {
              treg_outfile.open( treg_output_filename.c_str(), ofstream::out | ofstream::app );
              treg_detfile.open( treg_detput_filename.c_str(), ofstream::out | ofstream::app );
            }

            treg_output = &treg_outfile;
            treg_detput = &treg_detfile;

            if( local_trait_parser.get_output_options().export_out )
            {
              treg_csvput_filename = output_name_root + ".treg_export";

              if( !inside_tregs )
                treg_csvfile.open( treg_csvput_filename.c_str() );
              else
                treg_csvfile.open( treg_csvput_filename.c_str(), ofstream::out | ofstream::app );

              treg_csvput = &treg_csvfile;
            }
          }
          else
          {
            treg_output_filename = "traits.out";
            treg_detput_filename = "traits.det";

            if( !reg_file )
            {
              reg_file.reset( new ofstream(treg_output_filename.c_str()) );
              if( *reg_file )
                print_title(*reg_file);
            }
            treg_output = reg_file.get();

            if( !reg_det_file )
            {
              reg_det_file.reset( new ofstream(treg_detput_filename.c_str()) );
              if( *reg_det_file )
                print_title(*reg_det_file);
            }
            treg_detput = reg_det_file.get();

            if( local_trait_parser.get_output_options().export_out )
            {
              treg_csvput_filename = "traits.export";

              if( !reg_csv_file )
              {
                reg_csv_file.reset( new ofstream(treg_csvput_filename.c_str()) );
              }
              treg_csvput = reg_csv_file.get();
            }
          }

          if( !treg_output || !*treg_output )
          {
            cout << endl;
            errors << priority(error)
                   << "Cannot open output file: '" << treg_output_filename
                   << "'.  Skipping analysis..." << endl;
            continue;
          }

          if( !treg_detput || !*treg_detput )
          {
            cout << endl;
            errors << priority(error)
                   << "Cannot open detailed output file: '" << treg_detput_filename
                   << "'.  Skipping analysis..." << endl;
            continue;
          }

          if(    local_trait_parser.get_output_options().export_out
              && (!treg_csvput || !*treg_csvput) )
          {
            cout << endl;
            errors << priority(error)
                   << "Cannot open export output file: '" << treg_csvput_filename
                   << "'.  Skipping analysis..." << endl;
            continue;
          }

          cout << "Computing trait regression................" << flush;

          sibpal_analysis s_analysis(errors);

          s_analysis.run_regression_analysis(rp, local_trait_parser, *treg_output, *treg_detput, *treg_csvput);

          cout << "done." << endl;

          ++inside_tregs;

          treg_outfile.close();
          treg_detfile.close();

          if( local_trait_parser.get_output_options().export_out )
            treg_csvfile.close();
        }
        else if(    name == "MEAN_TEST"         || name == "MEANS_TEST"
                 || name == "MARKER_REGRESSION" || name == "MEAN_REGRESSION" )
        {
          cout << endl << "SIBPAL analysis......."
               << setw(3) << ++sibpals  << endl << endl << flush;

          meantest_parser local_mean_parser = global_mean_parser;

          local_mean_parser.parse_test_parameter_section(param);

          if( output_name_root.size() )
          {
            mean_output_filename = output_name_root + ".mean";

            if( !inside_means )
            {
              mean_outfile.open( mean_output_filename.c_str() );

              if( mean_outfile )
                print_title(mean_outfile);
            }
            else
              mean_outfile.open( mean_output_filename.c_str(), ofstream::out | ofstream::app );

            mean_output = &mean_outfile;

            if( local_mean_parser.csv_output() )
            {
              mean_csvput_filename = output_name_root + ".mean_export";

              if( !inside_means )
                mean_csvfile.open( mean_csvput_filename.c_str() );
              else
                mean_outfile.open( mean_output_filename.c_str(), ofstream::out | ofstream::app );

              mean_csvput = &mean_csvfile;
            }
          }
          else
          {
            mean_output_filename = "means.out";
            if( !mean_file )
            {
              mean_file.reset( new ofstream(mean_output_filename.c_str()) );
              if( *mean_file )
                print_title(*mean_file);
            }
            mean_output = mean_file.get();

            if( local_mean_parser.csv_output() )
            {
              mean_csvput_filename = "means.export";
              if( !mean_csv_file )
              {
                mean_csv_file.reset( new ofstream(mean_csvput_filename.c_str()) );
              }
              mean_csvput = mean_csv_file.get();
            }
          }

          if( !mean_output || !*mean_output )
          {
            cout << endl;
            errors << priority(error)
                   << "Cannot open output file: '" << mean_output_filename
                   << "'.  Skipping analysis..." << endl;
            continue;
          }

          if( local_mean_parser.csv_output() && (!mean_csvput || !*mean_csvput) )
          {
            cout << endl;
            errors << priority(error)
                   << "Cannot open export output file: '" << mean_csvput_filename
                   << "'.  Skipping analysis..." << endl;
            continue;
          }

          cout << "Computing mean test......................." << flush;

          sibpal_analysis s_analysis(errors);

          s_analysis.run_mean_test(rp, local_mean_parser, *mean_output, *mean_csvput);

          cout << "done." << endl;

          ++inside_means;

          mean_outfile.close();

          if( local_mean_parser.csv_output() )
            mean_csvfile.close();
        }
      }
    }
  }

  for( size_t anal = 0; anal < sdata.treg_analysis().size(); ++anal )
  {
    LSF_ptr<LSFBase> params = sdata.treg_analysis()[anal].first;
    string output_name_root = sdata.treg_analysis()[anal].second;

    if( !params ) continue;

    regression_parser local_trait_parser = global_trait_parser;

    regression_type local_reg_type = default_reg_type;
    
    if( has_attr(params, "single") || has_attr(params, "single_marker") )
      local_reg_type = SINGLE_MARKER;
    else if( has_attr(params, "multiple") || has_attr(params, "multiple_marker") )
      local_reg_type = MULTIPLE_MARKER;
    else if( has_attr(params, "zero") || has_attr(params, "zero_marker") )
      local_reg_type = ZERO_MARKER;

    cout << endl << "SIBPAL analysis......."
         << setw(3) << ++sibpals  << endl << endl << flush;

    local_trait_parser.set_regression_type(local_reg_type);
    local_trait_parser.parse_test_parameter_section(params);

    // Check if valid pair_info_file exist, then read the file.
    //
    if( local_trait_parser.get_pair_info_file().size() )
      PALBASE::read_pair_info_file(local_trait_parser.get_pair_info_file(), rp, errors, info_file);

    //local_trait_parser.dump_parser(cout);

    ostream* output = NULL;
    ostream* csvput = NULL;
    ostream* detput = NULL;

    ofstream outfile;
    ofstream csvfile;
    ofstream detfile;

    string   output_filename;
    string   csvput_filename;
    string   detput_filename;

    if( output_name_root.size() )
    {
      output_filename = output_name_root + ".treg";
      detput_filename = output_name_root + ".treg_det";

      outfile.open( output_filename.c_str() );
      detfile.open( detput_filename.c_str() );

      if( outfile )
        print_title(outfile);
      output = &outfile;

      if( detfile )
        print_title(detfile);
      detput = &detfile;

      if( local_trait_parser.get_output_options().export_out )
      {
        csvput_filename = output_name_root + ".treg_export";
        csvfile.open( csvput_filename.c_str() );

        csvput = &csvfile;
      }
    }
    else
    {
      output_filename = "traits.out";
      if( !reg_file )
      {
        reg_file.reset( new ofstream(output_filename.c_str()) );
        if( *reg_file )
          print_title(*reg_file);
      }
      output = reg_file.get();

      detput_filename = "traits.det";
      if( !reg_det_file )
      {
        reg_det_file.reset( new ofstream(detput_filename.c_str()) );
        if( *reg_det_file )
          print_title(*reg_det_file);
      }
      detput = reg_det_file.get();

      if( local_trait_parser.get_output_options().export_out )
      {
        csvput_filename = "traits.export";
        if( !reg_csv_file )
        {
          reg_csv_file.reset( new ofstream(csvput_filename.c_str()) );
        }
        csvput = reg_csv_file.get();
      }
    }

    if( !output || !*output )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << output_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( !detput || !*detput )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open detailed output file: '" << detput_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if(    local_trait_parser.get_output_options().export_out
        && (!csvput || !*csvput) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open export output file: '" << csvput_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    cout << "Computing trait regression................" << flush;

    sibpal_analysis s_analysis(errors);

    s_analysis.run_regression_analysis(rp, local_trait_parser, *output, *detput, *csvput);

    cout << "done." << endl;
  }

  for( size_t anal = 0; anal < sdata.mean_analysis().size(); ++anal )
  {
    LSF_ptr<LSFBase> params = sdata.mean_analysis()[anal].first;
    string output_name_root = sdata.mean_analysis()[anal].second;

    if( !params ) continue;

    cout << endl << "SIBPAL analysis......."
         << setw(3) << ++sibpals  << endl << endl << flush;

    meantest_parser local_mean_parser = global_mean_parser;

    local_mean_parser.parse_test_parameter_section(params);

    ostream* output = NULL;
    ostream* csvput = NULL;

    ofstream outfile;
    ofstream csvfile;

    string   output_filename;
    string   csvput_filename;

    if( output_name_root.size() )
    {
      output_filename = output_name_root + ".mean";
      outfile.open( output_filename.c_str() );

      if( outfile )
        print_title(outfile);
      output = &outfile;

      if( local_mean_parser.csv_output() )
      {
        csvput_filename = output_name_root + ".mean_export";
        csvfile.open( csvput_filename.c_str() );

        csvput = &csvfile;
      }
    }
    else
    {
      output_filename = "means.out";
      if( !mean_file )
      {
        mean_file.reset( new ofstream(output_filename.c_str()) );
        if( *mean_file )
          print_title(*mean_file);
      }
      output = mean_file.get();

      if( local_mean_parser.csv_output() )
      {
        csvput_filename = "means.export";
        if( !mean_csv_file )
        {
          mean_csv_file.reset( new ofstream(csvput_filename.c_str()) );
        }
        csvput = mean_csv_file.get();
      }
    }

    if( !output || !*output )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open output file: '" << output_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    if( local_mean_parser.csv_output() && (!csvput || !*csvput) )
    {
      cout << endl;
      errors << priority(error)
             << "Cannot open export output file: '" << csvput_filename
             << "'.  Skipping analysis..." << endl;
      continue;
    }

    cout << "Computing mean test......................." << flush;

    sibpal_analysis s_analysis(errors);

    s_analysis.run_mean_test(rp, local_mean_parser, *output, *csvput);

    cout << "done." << endl;
  }

  if( !sibpals )
  {
    ostream* treg_output = NULL;
    ostream* treg_detput = NULL;
    ostream* treg_csvput = NULL;

    ostream* mean_output = NULL;
    ostream* mean_csvput = NULL;

    if( !reg_file )
    {
      reg_file.reset( new ofstream("traits.out") );
      if( *reg_file )
        print_title(*reg_file);
    }
    treg_output = reg_file.get();

    if( !reg_det_file )
    {
      reg_det_file.reset( new ofstream("traits.det") );
      if( *reg_det_file )
        print_title(*reg_det_file);
    }
    treg_detput = reg_det_file.get();

    if( !mean_file )
    {
      mean_file.reset( new ofstream("means.out") );
      if(*mean_file)
        print_title(*mean_file);
    }
    mean_output = mean_file.get();

    cout << endl << "  No analyses specified."                            
         << endl << "Performing default analysis..."
         << endl << endl << flush;

    sibpal_analysis s_analysis(errors);

    cout << "Performing mean tests....................." << flush;

    s_analysis.run_mean_test(rp, global_mean_parser, *mean_output, *mean_csvput);

    cout << "done." << endl;

    cout << "Performing trait regressions.............." << flush;

    s_analysis.run_regression_analysis(rp, global_trait_parser, *treg_output, *treg_detput, *treg_csvput);

    cout << "done." << endl;
  }

  cout << endl << "Analysis complete!" << endl << endl;

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

} // end of namespace SIBPAL
} // end of namespace SAGE

int main(int argc, char* argv[])
{
  free(malloc(1));

  boost::shared_ptr<SAGE::SIBPAL::Sibpal> sibpal(new SAGE::SIBPAL::Sibpal(argc, argv));

  if( !sibpal.get() )
    exit(EXIT_FAILURE);

  sibpal->main();

  return EXIT_SUCCESS;
}
