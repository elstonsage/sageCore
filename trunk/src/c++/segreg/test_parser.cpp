//============================================================================
// File:      test_parser.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 5/16/2001.                                                   
//                                                                          
// Notes:     Tests the segreg parser class.
//    
//            Program options (argv[3]).
//              0 -> output all sub-models.
//              1 -> output non-sub-model options and mean, transformation and 
//                   frequency sub-models.
//              2 -> output transmission, resid, variance and fpmm sub-models.
//              3 -> call misc sub-model member functions which do calculations.
//              4 -> output sub-model dump() functions.
//              5 -> output onset sub-model.
//              6 -> output ascertainment sub-model.
//              7 -> call set_two_to_dom() and set_two_to_rec() for mean and 
//                   variance sub-models.
//              8 -> output prevalence constraint sub-model. 
//              9 -> output prevalence estimate sub-model. 
//             10 -> output type susceptibility sub-model. 
//             11 -> output primary trait type.
//             12 -> output covariate sub-models.
//             13 -> output covariates from prevalence sub-blocks.
//                                                                      
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================
#define protected public

#include <string>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "segreg/parser.h"
#include "segreg/RegPenetranceCommonAdjustments.h"
#include "util/StringUtils.h"
#include "output/Output.h"

using namespace std;
using namespace SAGE;
using namespace SEGREG;

inline ostream& operator<< (ostream& out, const MAXFUN::TransformationSubmodel& sm)
{
  int  old_precision = out.precision();
  out.precision(12);

  out << "\n" << sm.name() << " values: \n";
  out << "Option: " << sm.option_description() << std::endl;
  vector<MAXFUN::ParameterInput>::const_iterator  iter;
  for(iter = sm.my_parameters.begin(); iter != sm.my_parameters.end(); ++iter)
  {
    out << (OUTPUT::Table()
        << (OUTPUT::TableRow() << "maxfun parameter" << iter->param_name)
        << (OUTPUT::TableRow() << "value"            << iter->initial_estimate)
        << (OUTPUT::TableRow() << "status"           << MAXFUN::ParamTypeEnum2str(iter->initial_type))
        << (OUTPUT::TableRow() << "lower bound"      << iter->lower_bound)
        << (OUTPUT::TableRow() << "upper bound"      << iter->upper_bound));
  }

  out.precision(old_precision);
  
  return out;
}


enum print_option { ALL, G1, G2, CALC, DUMP, ONSET, ASCER, TWO, PREV_CONSTR, PREV_EST, 
                    TYPE_SUSC, TRAIT_TYPE, COV, PREV, INVALID };

print_option  
arg_2_print_option(const string& arg)
{
  print_option  option;

  if(arg == "0")
  {
    option = ALL;
  }
  else if(arg == "1")
  {
    option = G1;
  }
  else if(arg == "2")
  {
    option = G2;
  }
  else if(arg == "3")
  {
    option = CALC;
  }
  else if(arg == "4")
  {
    option = DUMP;
  }
  else if(arg == "5")
  {
    option = ONSET;
  }
  else if(arg == "6")
  {
    option = ASCER;
  }
  else if(arg == "7")
  {
    option = TWO;
  }
  else if(arg == "8")
  {
    option = PREV_CONSTR;
  }
  else if(arg == "9")
  {
    option = PREV_EST;
  }
  else if(arg == "10")
  {
    option = TYPE_SUSC;
  }
  else if(arg == "11")
  {
    option = TRAIT_TYPE;
  }
  else if(arg == "12")
  {
    option = COV;
  }
  else if(arg == "13")
  {
    option = PREV;
  }
  else
  {
    option = INVALID;
  }
  
  return option;
}

int main(int argc, char* argv[])
{
  //print_title(cout);
  
  //SAGEapp::expire();

  if (! (argc == 4 || argc == 5))
  {
    cerr << "usage: " << argv[0] << " <parameters> <pedigree> < test number> [output file]"
         << endl << endl
         << "Command line parameters:"             << endl
         << "  parameters   - Parameter File"      << endl
         << "  pedigree     - Pedigree Data File"  << endl
         << "  test output  - Integer Indicating The Output" << endl
         << "  output file  - File to append out to" << endl
         << endl << endl;
    exit(EXIT_FAILURE);
  }
  
  print_option  option = arg_2_print_option(string(argv[3]));
  if(option == INVALID)
  {
    cerr << "Option '" << argv[3] << "' not understood." << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTSEGREGPARSER-%P: ");

  // - Read parameter file.
  //
  cout << "Loading parameters..............";
  LSFBase *params = loadLSFfile(argv[1], "TESTSEGREGPARSER Parameter file", sage_cerr, false);
  cout << "done." << endl;

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }
  
  ofstream info_file;
  info_file.open("testparser.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testparser.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  if(argc == 5)
  {
    out_file.open(argv[4], ios::out | ios::app);
  }
  else
  {
    out_file.open("testparser.out");
  }

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testparser.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  //print_title(info_file);
  
  if(option != DUMP)
  {
    //print_title(out_file);
  }

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);
  error_file.prefix("%%TESTSEGREGPARSER-%P:");
  errors.insert(sage_cerr);
  errors.restrict(r_ge, error);
  errors.insert(error_file);
  errors.restrict(r_ge, information);
  
  // - Read the pedigree data.
  //
  RPED::RefMultiPedigree p;
  bool pedigree_loaded = false;
  LSFList::const_iterator i;
  AttrVal a;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !*i ) continue;

    if( toUpper( (*i)->name() ) == "PEDIGREE" )
    {
      if( (*i)->attrs() && (*i)->attrs()->has_attr("column") )
      {
        RPED::RefLSFFortranPedigreeFile ped_reader(errors);

        ped_reader.set_force_skip_markers(false);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(true);

        ped_reader.process_parameters(p.info(), *i);

        cout << "Reading pedigrees..............." << flush;
        if( !ped_reader.input(p, argv[2], info_file) )
        {
          errors << priority(fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, argv[2], info_file, true, false);
        
        break;
      }
      else
      {
        RPED::RefLSFDelimitedPedigreeFile ped_reader(errors);

        ped_reader.set_force_skip_markers(false);
        ped_reader.set_force_skip_traits(false);
        ped_reader.set_force_dynamic_markers(true);

        ped_reader.process_parameters(p.info(), *i);

        cout << "Reading pedigrees..............." << flush;
        if( !ped_reader.input(p,argv[2], info_file) )
        {
          errors << priority(fatal) << "Error reading pedigree data" << endl;
          exit(EXIT_FAILURE);
        }
        pedigree_loaded = true;
        
        // Need to be added to print 10 individuals. - yjs 
        //
        ped_reader.print_mped(p, argv[2], info_file, true, false);
        
        break;
      }
    }
  }  

  if( !pedigree_loaded )
  {
    errors << priority(fatal) << "Fatal Error: No pedigree specified!  Terminating..." << endl;
    exit(EXIT_FAILURE);
  }

  cout << "done." << endl;

  if( !p.pedigree_count() )
  {
    errors << priority(critical) << "No pedigrees to analyze... aborting." << endl;
    exit(EXIT_FAILURE);
  }

  // - Check for errors in pedigree data.
  //
  size_t error_count = 0;
  RPED::RefMultiPedigree::pedigree_const_iterator j;
  for( j = p.pedigree_begin(); j != p.pedigree_end(); ++j)
    error_count += j->error_count();

  if( error_count )
    errors << priority(error) << "Errors appear in pedigree data.  " 
           << "Results may be incomplete." << endl;
           

  RPED::RefMultiPedigree::pedigree_iterator jj;
  for( jj = p.pedigree_begin(); jj != p.pedigree_end(); ++jj)
  {
    PedigreeSort( *jj );
  }
  
  
  // - Parse segreg_analysis blocks.
  //
  parser  test_parser(&p, cout, errors);
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !*i ) continue;

    string  block_name = toUpper((*i)->name());
    
    if(block_name == "SEGREG_ANALYSIS" ||
       block_name == "SEGREG"            )
    {
      size_t  id = test_parser.analysis_number() + 1;
    
      errors << priority(error) << endl;
      errors << priority(error) << "---------- " << id << " New Analysis Block " << id << " ---------- " << endl;
      errors << priority(error) << endl;
      
      test_parser.parse_test_parameter_section(*i);
      
      if(option != DUMP)
      {
        out_file << endl << "---------- " << id << " New Analysis Block " << id << " ---------- " << endl << endl;
      }
      
      if(option == ALL)
      {
        out_file << "mean missing: " << boolalpha << test_parser.get_model().get_type_missing() << "\n"
                 << "trans missing: " << test_parser.get_model().get_trans_missing() << noboolalpha
                 << endl;
      }
      
      if(option == G1 || option == ALL)
      {      
        out_file << "File name root: " << test_parser.get_model().get_file_name_root() << endl;
      }
      
      if(option != DUMP)
      {
        out_file << "Model class: " << model_class_2_string(test_parser.get_model().get_model_class()) << endl;
      }
      
      if(option == G1 || option == ALL)
      {
        out_file << "Each pedigree: "               << test_parser.get_model().get_each_pedigree()                      << endl
                 << "Penetration function output: " << test_parser.get_model().get_pen_func_output()                    << endl
                 << "Type probability: "            << test_parser.get_model().get_type_prob()                          << endl
                 << "Primary trait: "               << test_parser.get_model().get_primary_trait()              << "\n" << endl
                 <<                                    test_parser.get_model().mean_sub_model
                 << "mean default:  "               << boolalpha << test_parser.get_model().mean_sub_model.is_default() << endl
                 <<                                    test_parser.get_model().transf_sub_model
                 <<                                    test_parser.get_model().freq_sub_model
                 << "freq default:  "               << boolalpha << test_parser.get_model().freq_sub_model.default_()   << endl
                 <<                                                                                                        endl;
      }
      
      // - Covariate sub-models.
      //
      if(option == COV || option == ALL)
      {
        out_file << test_parser.get_model().mean_cov_sub_model;

        vector<CovariateSubmodel::covariate> covs = test_parser.get_model().mean_cov_sub_model.covariates();
        
        for(int i = 0; i < (int)covs.size(); ++i)
        {
          out_file << covs[i];
        }
        
        out_file << test_parser.get_model().var_cov_sub_model;

        covs = test_parser.get_model().var_cov_sub_model.covariates();

        for(int i = 0; i < (int)covs.size(); ++i)
        {
          out_file << covs[i];
        }
        
        out_file << test_parser.get_model().susc_cov_sub_model;

        covs = test_parser.get_model().susc_cov_sub_model.covariates();

        for(int i = 0; i < (int)covs.size(); ++i)
        {
          out_file << covs[i];
        }
        
        out_file << test_parser.get_model().comp_trait_sub_model;

        covs = test_parser.get_model().comp_trait_sub_model.covariates();

        for(int i = 0; i < (int)covs.size(); ++i)
        {
          out_file << covs[i];
        }
      }
      
      if(option == G2 || option == ALL)
      {
        out_file <<                                     test_parser.get_model().transm_sub_model
                 << "transm default:  " << boolalpha << test_parser.get_model().transm_sub_model.default_() << endl
                 <<                                     test_parser.get_model().resid_sub_model
                 << "resid default:  "  << boolalpha << test_parser.get_model().resid_sub_model.is_in_default_mode()  << endl
                 <<                                     test_parser.get_model().var_sub_model
                 << "var default:  "    << boolalpha << test_parser.get_model().var_sub_model.is_default()  << endl
                 <<                                     test_parser.get_model().fpmm_sub_model
                 << "fpmm frequency: "               << test_parser.get_model().fpmm_sub_model.frequency()  << endl
                 << "fpmm loci: "                    << test_parser.get_model().fpmm_sub_model.loci()       << endl
                 <<                                                                                            endl;
      }
      
      if(option == CALC)
      {
        out_file << test_parser.get_model().resid_sub_model
                 << "alpha mother: father's trait known, mother's trait known    "
                 << test_parser.get_model().resid_sub_model.alpha_mother(2.0,2.0) << endl
                 << "alpha mother: father's trait missing, mother's trait known    "
                 << test_parser.get_model().resid_sub_model.alpha_mother(QNAN,2.0) << endl
                 << "alpha mother: father's trait missing, mother's trait missing    "
                 << test_parser.get_model().resid_sub_model.alpha_mother(QNAN,QNAN) << endl
                 << "alpha father: father's trait known, mother's trait known    "
                 << test_parser.get_model().resid_sub_model.alpha_father(2.0,2.0) << endl
                 << "alpha father: father's trait known, mother's trait missing    "
                 << test_parser.get_model().resid_sub_model.alpha_father(2.0,QNAN) << endl
                 << "alpha father: father's trait missing, mother's trait missing    "
                 << test_parser.get_model().resid_sub_model.alpha_father(QNAN,QNAN) << endl
                 << test_parser.get_model().freq_sub_model
                 << test_parser.get_model().transm_sub_model
                 << (OUTPUT::Table()
                 << (OUTPUT::TableRow() << "genotype probability: indiv AA, mother BB, father AA" << test_parser.get_model().transm_sub_model.prob(index_AA, index_BB, index_AA)))
                 << test_parser.get_model().transf_sub_model;
  
        std::vector<double>  traits;
        traits.push_back(1);
        traits.push_back(QNAN);
        traits.push_back(3);
        test_parser.get_model().transf_sub_model.calculate_geom_mean(traits);
        test_parser.get_model().transf_sub_model.transform(traits);  
        
        OUTPUT::Table t;
        
        for(std::vector<double>::const_iterator iter = traits.begin(); iter != traits.end(); ++iter)
        {
          t << (OUTPUT::TableRow() << *iter);
        }
        
        out_file << t;
      }
      
      
      // - Conditions below meant to insure that only non default 
      //   sub-models are dumped.
      //
      if(option == DUMP)
      {
        if(! SAGE::isnan(test_parser.get_model().fpmm_sub_model.variance()))
        {
          out_file << "\n";
          test_parser.get_model().fpmm_sub_model.dump(out_file);
        }
        
        // - Violates original intention of only dumping if user
        //   actually specifies something, but does not cause a problem.
        //
        if(true)
        {
          out_file << "\n";
          test_parser.get_model().transf_sub_model.dump(out_file);
        }
        
        if(! test_parser.get_model().mean_sub_model.is_default())
        {
          out_file << "\n";
          test_parser.get_model().mean_sub_model.dump(out_file);
        }
        
        if(! test_parser.get_model().var_sub_model.is_default())
        {
          out_file << "\n";
          test_parser.get_model().var_sub_model.dump(out_file);
        }
        
        if(! test_parser.get_model().resid_sub_model.is_in_default_mode())
        {
          out_file << "\n";
          test_parser.get_model().resid_sub_model.dump(out_file);

          test_common_reg_corr_adjustments(out_file,
                                           test_parser.get_model().get_model_class(),
                                           test_parser.get_model().resid_sub_model);
        }
        
        if(! test_parser.get_model().freq_sub_model.default_())
        {
          out_file << "\n";
          test_parser.get_model().freq_sub_model.dump(out_file);
        }
        
        if(test_parser.get_model().comp_trait_sub_model.covariates().size())
        {
          out_file << "\n";
          test_parser.get_model().comp_trait_sub_model.dump(out_file);
        }
        
        if(test_parser.get_model().mean_cov_sub_model.covariates().size())
        {
          out_file << "\n";
          test_parser.get_model().mean_cov_sub_model.dump(out_file);
        }
        
        if(test_parser.get_model().var_cov_sub_model.covariates().size())
        {
          out_file << "\n";
          test_parser.get_model().var_cov_sub_model.dump(out_file);
        }
        
        if(test_parser.get_model().susc_cov_sub_model.covariates().size())
        {
          out_file << "\n";
          test_parser.get_model().susc_cov_sub_model.dump(out_file);
        }
        
        out_file << "\n";
        test_parser.get_model().transm_sub_model.dump(out_file);
        
        if(ONSET_AVAILABLE)
        {
          out_file << "\n";
          test_parser.get_model().ons_sub_model.dump(out_file);
        }
        
        out_file << "\n";
        test_parser.get_model().ascer_sub_model.dump(out_file);
      }
      
      if(option == ONSET || option == ALL)
      {
        out_file << test_parser.get_model().ons_sub_model
                 << "type option: " << test_parser.get_model().ons_sub_model.t_option_description() << "\n"
                 << "multi option: " << test_parser.get_model().ons_sub_model.m_option_description() << "\n"
                 << "affection status: " << test_parser.get_model().ons_sub_model.affection_status() << "\n"
                 << "age of onset: " << test_parser.get_model().ons_sub_model.age_of_onset() << "\n"
                 << "age at exam: " << test_parser.get_model().ons_sub_model.age_at_exam() << "\n"
                 << endl;

        if(option == ONSET)
        {
          out_file << "mean missing: " << boolalpha << test_parser.get_model().get_type_missing() << "\n"
                   << "trans missing: " << test_parser.get_model().get_trans_missing() << noboolalpha
                   << endl << endl
                   << test_parser.get_model().ascer_sub_model;
        }
      }
      
      if(option == ASCER || option == ALL)
      {
        out_file << test_parser.get_model().ascer_sub_model
                 << "psf indicator: " << test_parser.get_model().ascer_sub_model.psf_indicator() << "\n"
                 << (OUTPUT::Table()
                 << (OUTPUT::TableRow() << "thresh"           << test_parser.get_model().ascer_sub_model.thresh())
                 << (OUTPUT::TableRow() << "thresh indicator" << test_parser.get_model().ascer_sub_model.thresh_indicator())
                 << (OUTPUT::TableRow() << "thresh high"      << test_parser.get_model().ascer_sub_model.thresh_high())
                 << (OUTPUT::TableRow() << "thresh low"       << test_parser.get_model().ascer_sub_model.thresh_low()))
                 << "includes: ";
        
        const vector<double>&  includes = test_parser.get_model().ascer_sub_model.includes();
        vector<double>::const_iterator  iter;
        for(iter = includes.begin(); iter != includes.end(); iter++)
        {
          out_file << *iter << "  ";
        }

        out_file << endl << "indic_thresh: " << test_parser.get_model().ascer_sub_model.indic_thresh() << "\n" << endl;
      }
      
      // - Type susceptibility.
      //
      if(option == TYPE_SUSC || option == ALL)
      {
        out_file << test_parser.get_model().susc_sub_model << endl;
      }
      
      if(option == TWO)
      {
        // - mean sub-model.
        //
        genotype_specific_mean_sub_model  mean_sm = model(test_parser.get_model()).mean_sub_model;

        out_file << (mean_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << mean_sm;
        
        mean_sm.my_two_is_dom = false;

        out_file << (mean_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << mean_sm;
        
        mean_sm.my_two_is_dom = true;

        out_file << (mean_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << mean_sm << endl;
        
        
        // - variance sub-model.
        //
        genotype_specific_variance_sub_model  var_sm = model(test_parser.get_model()).var_sub_model;

        out_file << (var_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << var_sm;
        
        var_sm.my_two_is_dom = false;

        out_file << (var_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << var_sm;
        
        var_sm.my_two_is_dom = true;

        out_file << (var_sm.is_two_dom() ? "dominant" : "recessive") << "\n" << var_sm << endl;
      }
      if(option == TRAIT_TYPE || option == ALL)
      {
        out_file << "Primary trait type:  "  << primary_type_2_string(test_parser.get_model().get_primary_trait_type()) << endl;
      }
      
      if(option == PREV || option == ALL)
      {
        const prevalence_sub_model& psm = test_parser.get_model().prev_sub_model;
        out_file << "\nPrevalence - \n";

        for(size_t j = 0; j < psm.get_estimate_count(); ++j)
        {
          for(int i = 0; i < (int)psm.get_estimate_susc_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(S): " << psm.get_estimate_susc_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_estimate_susc_covariate_value(j, i) << "\n"
                     << std::endl;
          }
          for(int i = 0; i < (int)psm.get_estimate_mean_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(M): " << psm.get_estimate_mean_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_estimate_mean_covariate_value(j, i) << "\n"
                     << std::endl;
          }
          for(int i = 0; i < (int)psm.get_estimate_var_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(V): " << psm.get_estimate_var_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_estimate_var_covariate_value(j, i) << "\n"
                     << std::endl;
          }
          if(test_parser.get_model().get_primary_trait_type() == pt_ONSET)
            out_file << (OUTPUT::Table() << (OUTPUT::TableRow() << "AGE" << psm.get_estimate_age(j)));
        }

        for(size_t j = 0; j < psm.get_constraint_count(); ++j)
        {
          out_file << "P Constraint: " << j << endl;
          for(int i = 0; i < (int)psm.get_constraint_susc_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(S): " << psm.get_constraint_susc_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_constraint_susc_covariate_value(j, i) << "\n"
                     << std::endl;
          }
          for(int i = 0; i < (int)psm.get_constraint_mean_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(M): " << psm.get_constraint_mean_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_constraint_mean_covariate_value(j, i) << "\n"
                     << std::endl;
          }
          for(int i = 0; i < (int)psm.get_constraint_var_covariate_count(j); ++i)
          {
            out_file << "\n" << "Name(V): " << psm.get_constraint_var_covariate_name(j, i) << "\n"
                     << "Coefficient: " << psm.get_constraint_var_covariate_value(j, i) << "\n"
                     << std::endl;
          }
        
          if(test_parser.get_model().get_primary_trait_type() == pt_ONSET)
            out_file << (OUTPUT::Table() << (OUTPUT::TableRow() << "AGE" << psm.get_constraint_age(j)));

          out_file << "R:  " << psm.get_constraint_number_affected(j) << endl;
          out_file << "N:  " << psm.get_constraint_sample_size(j) << endl;
        }
        
        out_file << endl;
      }
    }
  }
}

