#ifndef FCOR_ANALYSIS_H
#define FCOR_ANALYSIS_H

//****************************************************************************
//* File:      analysis.h                                                    *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                  Jan 31 01 *
//*                                                                          *
//* Notes:     This header file defines fcor_analysis.                       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/output.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     FcorAnaysis                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines functions for fcor_analysis.                         ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class FcorAnalysis
{
  public:

    FcorAnalysis(const FcorParser& f_parser, cerrorstream& err = sage_cerr);
   
    void  run_analysis(ostream& sum_out, ostream& det_out, ostream& pair_out,
                       ostream& cov_out, ostream& xls_out);

  protected:

    void  compute_correlations(const pairset_vector&        pairset,
                               const pairset_info_vector&   pinfos,
                                     CorrelationCal&        corinfo,
                                     pairset_result_vector& results);

    void  compute_standard_errors(const CorrelationCal&        corinfo,
                                        pairset_result_vector& results);

    void  do_KE_analysis(const pairset_vector&        pairset,
                         const pairset_info_vector&   pinfos,
                               pairset_result_vector& results);

    void  view_results(const pairset_info_vector&   pinfos,
                       const pairset_result_vector& results,
                       ostream& out) const;

    const FcorParser*               my_parser;

    PairSetData                     my_pairsetdata;

    CorrelationCal                  my_subtype_corinfo;
    CorrelationCal                  my_maintype_corinfo;

    pairset_result_vector           my_subtype_result;
    pairset_result_vector           my_maintype_result;

    Htest_result_vector             my_Htest_result;
    var_cov_result_vector           my_var_cov_result;

    cerrorstream                    errors;
};

} // end of namespace FCOR
} // end of namespace SAGE

#endif
