//****************************************************************************
//* File:      analysis.cpp                                                  *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                  Feb 08 01 *
//*                                                                          *
//* Notes:     This source file has fcor_analysis run.                       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/analysis.h"

namespace SAGE {
namespace FCOR {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     FcorAnalysis                                                ~
// ~                                                                        ~
// ~ Purpose:   Run fcor analysis.                                          ~
// ~                                                                        ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//--------------------------------------------------------------------------- 
// Out-of-Line Implementation of FcorAnalysis
//--------------------------------------------------------------------------- 

FcorAnalysis::FcorAnalysis(const FcorParser& f_parser, cerrorstream& err)
{
  my_parser = &f_parser;
  errors    = err;
}

void
FcorAnalysis::run_analysis(ostream& sum_out, ostream& det_out, ostream& pair_out,
                           ostream& cov_out, ostream& xls_out)
{
  my_pairsetdata.build_pairset(*my_parser);

  bool sub_cor  = false;
  bool main_cor = false;

  // subtype.
  //
  if( my_parser->get_analysis_options().pairset != MAINTYPES )
  {
    cout << endl << "Computing subtype correlations............" << flush;

    compute_correlations(my_pairsetdata.get_subtype_pairset(),
                         my_pairsetdata.get_subtype_info(),
                         my_subtype_corinfo,
                         my_subtype_result);

    cout << "done." << endl;

    if( my_parser->get_analysis_options().standard_error )
    {
      cout << "Computing standard errors................." << flush;

      compute_standard_errors(my_subtype_corinfo, my_subtype_result);

      cout << "done." << endl;  
    }

    sub_cor = true;
    //view_results(my_pairsetdata.get_subtype_info(), my_subtype_result, cout);
  }

  // maintype.
  //
  if( my_parser->get_analysis_options().pairset != SUBTYPES )
  {
    cout << endl << "Computing main type correlations.........." << flush;

    compute_correlations(my_pairsetdata.get_maintype_pairset(),
                         my_pairsetdata.get_maintype_info(),
                         my_maintype_corinfo,
                         my_maintype_result);

    cout << "done." << endl;

    if( my_parser->get_analysis_options().standard_error )
    {
      cout << "Computing standard errors................." << flush;

      compute_standard_errors(my_maintype_corinfo, my_maintype_result);

      cout << "done." << endl;  
    }

    main_cor = true;
    //view_results(my_pairsetdata.get_maintype_info(), my_maintype_result, cout);

    if( my_parser->get_analysis_options().KE_with_optimal )
    {
      cout << "Computing KE standard errors.............." << flush;

      do_KE_analysis(my_pairsetdata.get_maintype_pairset(),
                     my_pairsetdata.get_maintype_info(),
                     my_maintype_result);

      cout << "done." << endl;  
    }
    //view_results(my_pairsetdata.get_maintype_info(), my_maintype_result, cout);
  }

  // homogeneity test.
  //
  if( my_parser->get_analysis_options().homogeneity_test )
  {
    cout << endl << "Peforming homogeneity test................" << flush;

    if( !sub_cor )
    {
      compute_correlations(my_pairsetdata.get_subtype_pairset(),
                           my_pairsetdata.get_subtype_info(),
                           my_subtype_corinfo,
                           my_subtype_result);

      if( my_parser->get_analysis_options().standard_error )
        compute_standard_errors(my_subtype_corinfo, my_subtype_result);
    }

    Htest  h(my_pairsetdata, my_subtype_corinfo);
    h.compute_homogeneity_test(my_subtype_result, my_Htest_result);

    cout << "done." << endl;
  }

  // print output files.
  //
  cout << "Writing output files......................" << flush;

  FcorView fv(my_pairsetdata);

  fv.print_analysis_results(sum_out, SUMMARY, my_subtype_result, my_maintype_result, my_Htest_result);

  if( my_parser->get_analysis_options().detailed_output )
    fv.print_analysis_results(det_out, DETAILED, my_subtype_result, my_maintype_result, my_Htest_result);

  if( my_parser->get_analysis_options().pair_output )
    fv.print_analysis_results(pair_out, PAIR_COUNT, my_subtype_result, my_maintype_result, my_Htest_result);

  if( my_parser->get_analysis_options().xls_output )
    fv.print_analysis_results(xls_out, XLS_FORMAT, my_subtype_result, my_maintype_result, my_Htest_result);

  cout << "done." << endl;

  // variance-covariance matrix.
  //
  if( my_parser->get_analysis_options().var_covs.size() )
  {
    cout << endl << "Computing variance-covariance matrices...." << flush;

    if( !sub_cor )
    {
      compute_correlations(my_pairsetdata.get_subtype_pairset(),
                           my_pairsetdata.get_subtype_info(),
                           my_subtype_corinfo,
                           my_subtype_result);

      if( my_parser->get_analysis_options().standard_error )
        compute_standard_errors(my_subtype_corinfo, my_subtype_result);
    }

    if( !main_cor )
    {
      compute_correlations(my_pairsetdata.get_maintype_pairset(),
                           my_pairsetdata.get_maintype_info(),
                           my_maintype_corinfo,
                           my_maintype_result);

      if( my_parser->get_analysis_options().standard_error )
        compute_standard_errors(my_maintype_corinfo, my_maintype_result);
    }

    VarCovCal  vcc(my_subtype_corinfo, my_maintype_corinfo);
    vcc.compute_variance_covariances(my_subtype_result, my_maintype_result, my_var_cov_result);

    cout << "done." << endl;

    cout << "Writing output files......................" << flush;

    fv.print_var_cov_results(cov_out, my_var_cov_result);

    cout << "done." << endl;
  }
}

void
FcorAnalysis::compute_correlations(const pairset_vector&        pairset,
                                   const pairset_info_vector&   pinfos,
                                         CorrelationCal&        corinfo,
                                         pairset_result_vector& results)
{
  corinfo.set_parser(*my_parser);
  corinfo.compute_correlations(pairset, pinfos, results);
  //corinfo.view_corinfo(cout);
}

void
FcorAnalysis::compute_standard_errors(const CorrelationCal&        corinfo,
                                            pairset_result_vector& results)
{
  StdErrCal  se(corinfo);
  se.compute_standard_errors(results);
}

void
FcorAnalysis::do_KE_analysis(const pairset_vector&        pairsets,
                             const pairset_info_vector&   pinfos,
                                   pairset_result_vector& results)
{
  // Compute pedigree-specific optimal weights.
  //
  optimal_weight_finder owf;
  owf.find_optimal_weights(pairsets, results);

  const weight_matrix_vector& oweight = owf.get_weight_matrix();

  // Compute KE correlations with optimal weights.
  //
  CorrelationCal KE_correlation;
  KE_correlation.set_parser(*my_parser);
  KE_correlation.compute_correlations(pairsets, pinfos, oweight, results);
  //KE_correlation.view_corinfo(cout);

  // Compute KE variances with optimal weights.
  //
  StdErrCal  KE_se(KE_correlation);
  KE_se.compute_standard_errors(oweight, results);
}

void
FcorAnalysis::view_results(const pairset_info_vector&   pinfos,
                           const pairset_result_vector& results,
                           ostream& out) const
{
  for( size_t r = 0; r < results.size(); ++r )
  {
    out << pinfos[r].gname << endl;

    const Matrix2D< correlation_result >&    cor = results[r].corr;
    const Matrix2D< standard_error_result >& std = results[r].std_err;

    for( size_t t1 = 0; t1 < cor.rows(); ++t1 )
    {
      for( size_t t2 = 0; t2 < cor.cols(); ++t2 )
      {
        out << setw(5) << cor(t1, t2).count << " ";

        for( size_t w = 0; w < cor(t1, t2).correlation.size(); ++w )
          out << setw(5) << cor(t1, t2).correlation[w] << " ";
      }
      out << endl;
    }
    out << endl;

    for( size_t t1 = 0; t1 < std.rows(); ++t1 )
    {
      for( size_t t2 = 0; t2 < std.cols(); ++t2 )
      {
        out << setw(5) << std(t1, t2).least_pair_count << " "
            << setw(3) << std(t1, t2).replaced << " "
            << setw(5) << std(t1, t2).w1 << " ";

        for( size_t w = 0; w < std(t1, t2).standard_error.size(); ++w )
          out << setw(5) << std(t1, t2).standard_error[w] << " ";
      }
      out << endl;
    }
    out << endl;
  }
  out << endl;
}

// end of FcorAnalysis Implementation

} // end of namespace FCOR
} // end of namespace SAGE
