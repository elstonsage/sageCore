//============================================================================
// File:      lodlink.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 9/20/2
//                                                                          
// Notes:     lodlink class and program.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/lodlink.h"

using namespace std;

namespace SAGE
{

namespace LODLINK
{

void
build_headers()
{
  cleves_elston_result::build_headers();
  non_ss_lod_ratio_result::build_headers();
  ss_lod_ratio_result::build_headers();
  non_ss_smiths_result::build_headers();
  ss_smiths_result::build_headers();
  non_ss_faraways_result::build_headers();
  ss_faraways_result::build_headers();  
  non_ss_mortons_result::build_headers();
  ss_mortons_result::build_headers();    
  non_ss_genotype_result::build_headers();
  ss_genotype_result::build_headers();
}

//============================================================================
// IMPLEMENTATION:  lodlink
//============================================================================
//
lodlink::lodlink(int argc, char** argv)
       : APP::SAGEapp(APP::APP_LODLINK, true, argc, argv)
{
  LSFInit();
}

void lodlink::print_title(ostream& o)
{
  SAGEapp::print_title(o);
  o << "\n" << endl;
}

void
lodlink::check_for_loops(const RPED::RefMultiPedigree& mped, cerrorstream& errors)
{
  typedef SAGE::RPED::RefMultiPedigree::pedigree_const_iterator  pedigree_const_iterator;
  typedef SAGE::RPED::RefMultiPedigree::subpedigree_const_iterator  subpedigree_const_iterator;

  set<string>  loop_subpeds;
  
  pedigree_const_iterator  ped_iter = mped.pedigree_begin();
  for(; ped_iter != mped.pedigree_end(); ++ped_iter)
  {
    subpedigree_const_iterator  subped_iter = ped_iter->subpedigree_begin();
    for(; subped_iter != ped_iter->subpedigree_end(); ++subped_iter)
    {
      if(MPED::mp_utilities::has_loops(*subped_iter))
      {
        string  ped_name                  = subped_iter->pedigree()->name();
        string  member_name               = subped_iter->member_begin()->name();
        string  constituent_pedigree_name = "  constituent pedigree in pedigree " + ped_name +
                                            " containing member " + member_name;
        loop_subpeds.insert(constituent_pedigree_name);
      }
    }
  }
  
  if(! loop_subpeds.empty())
  {
    string  loop_names;
    set<string>::const_iterator  iter     = loop_subpeds.begin();
    set<string>::const_iterator  end_iter = loop_subpeds.end();
    for(; iter != end_iter; ++iter)
    {
      loop_names += *iter + "\n"; 
    }
    
    errors << priority(critical) << "The following constituent pedigree(s) contain loops and/or marriage rings:\n\n" 
           << loop_names << "\n\nPlease remove loops and/or marriage rings (Run Pedinfo for help "
           << "in locating loops and marriage rings.), and rerun the program.  Exiting ..." 
           << endl;
    exit(EXIT_FAILURE);
  }
}

// - Check for allele frequencies of 0 and issue a warning for any that are
//   found.
//
void
lodlink::check_allele_freqs(const RPED::RefMultiPedigree& mped, cerrorstream& errors)
{
  const MLOCUS::inheritance_model_map&  loci = mped.info().markers();
  
  MLOCUS::inheritance_model_map::index_const_iterator  l_iter = loci.index_begin();
  MLOCUS::inheritance_model_map::index_const_iterator  l_end_iter = loci.index_end();
  for(; l_iter != l_end_iter; l_iter++)
  {
    MLOCUS::allele_iterator  a_iter = l_iter->second.allele_begin();
    MLOCUS::allele_iterator  a_end_iter = l_iter->second.allele_end();
    for(; a_iter != a_end_iter; a_iter++)
    {
      if(a_iter->frequency() == 0.0)
      {
        errors << priority(warning) << "Frequency of allele '" << a_iter->name() 
                                    << "' at locus '" << l_iter->first << "' is 0."
                                    << endl;
      }
    }
  }
}

int lodlink::main()
{
  // - Read data.
  //
  lodlink_data  my_data(name, debug());

  print_title(my_data.info());

  my_data.input(argc, argv);
  
  const RPED::RefMultiPedigree*  my_mped = &(my_data.pedigrees());
  assert(my_mped != 0);

  check_for_loops(*my_mped, my_data.errors());
  check_allele_freqs(*my_mped, my_data.errors());

  cout << "Generating statistics....................." << flush;

  build_headers();
  parser  my_parser(my_mped, my_data, cout, my_data.errors());
 
  // - Do the analysis specified in each parameter file analysis block.
  //
  bool  analysis_specified = false;
  for(size_t a = 0; a < my_data.analyses().size(); ++a)
  {
    LSF_ptr<LSFBase> analysis_ptr = my_data.analyses()[a];
    assert(analysis_ptr != 0);

    analysis_specified = true;

    my_parser.parse(analysis_ptr);
    
    ofstream  summary_file; 
    string    summary_file_name = my_parser.user_instructions().file_name_root + ".sum"; 
    summary_file.open(summary_file_name.c_str());
    if(! summary_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file '" << summary_file_name << "'.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }
    
    ofstream  detail_file; 
    string    detail_file_name = my_parser.user_instructions().file_name_root + ".det"; 
    detail_file.open(detail_file_name.c_str());
    if(! detail_file)
    {
      sage_cerr << priority(fatal) 
                << "Cannot open output file '" << detail_file_name << "'.  Exiting..." << endl;
      exit(EXIT_FAILURE);
    }    
    
    analysis  my_analysis(my_data.errors(), *my_mped, my_parser.user_instructions());
    my_analysis.build();
    my_analysis.analyze();
    my_analysis.write(summary_file, detail_file, *this);
    analysis::clear();    
    
    ge_models::clear_models();
  }

  if(! analysis_specified)
  {
    my_data.errors() << priority(error) << "Parameter file contains no analysis." << endl;
  }
  else
  {
    cout << "done." << endl;
    cout << endl << "Analysis complete!" << endl;
  }

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

}
}


//============================================================================
// LODLINK, the program
//============================================================================
//
int main(int argc, char* argv[])
{
  free(malloc(1));

  SAGE::LODLINK::lodlink lodlink_inst(argc, argv);

  lodlink_inst.main();
  
  exit(EXIT_SUCCESS);
}


