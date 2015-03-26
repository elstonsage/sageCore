#ifndef SEGREG_OUTPUT_H
#define SEGREG_OUTPUT_H
//****************************************************************************
//* File:      analysis_out.h                                                *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                        yjs May. 01 *
//*                                                                          *
//* Notes:     This header file defines analysis_out class for SEGREG.       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************


#include <iostream>
#include "segreg/analysis.h"
#include "error/errorstream.h"
#include "error/errormanip.h"

namespace SAGE {
namespace SEGREG {

/** The analysis_output class collects several static functions for use elsewhere.
 *  These functions output various tables of SEGREG results, and/or multiply
 *  used sections of code.
 */
class analysis_output
{
public:

  static void output_model            (ostream&, const model&);
  static void output_initial_estimates(ostream&, const primary_analysis_results&);
  static void output_final_estimates  (ostream&, const primary_analysis_results&);
  static void output_likelihoods      (ostream&, const primary_analysis_results&);
  static void output_final_error      (ostream&, const primary_analysis_results&);
  static void output_vc_matrix        (ostream&, const primary_analysis_results&);

  static void output_likelihood_table_header  (ostream&, model& test_model);
  static void output_likelihood_table_results (ostream&, const primary_analysis_results&);
  static void output_likelihood_table_footer  (ostream&);

  static void output_mlm_resid_corr_results   (ostream&, const primary_analysis_results&);

  static bool skip_poly_locus; // due to JA
  static double deriv_sum_sq; // due to JA
  
  /// Output the asymptotic_p_value table of segregation results
  ///
  /// \param out     The output stream
  /// \param results The results objects to include in the table
  static void output_asymptotic_p_value_results
      (ostream&                                out,
       const vector<primary_analysis_results>& results);  

protected:

  static void output_model_continuous (ostream&, const model&);
  static void output_model_binary     (ostream&, const model&);

};

//----------------------------------------------------------------------------
//  Class:     analysis_viewer
//
//  Purpose:   Base class for analysis output classes.
//----------------------------------------------------------------------------

class analysis_viewer
{
  public:

    // Required to make compiler happy.
    virtual ~analysis_viewer() { }

    analysis_viewer(ostream& output) : my_output(output) { }

    virtual void print_header(const primary_analysis_results&);
    virtual void print_results(const primary_analysis_results&) = 0;
    virtual void print_footer(const primary_analysis_results&)  = 0;

    ostream&          output_stream()        { return my_output; }

  protected:

    ostream& my_output;
};

//----------------------------------------------------------------------------
//  Class:     analysis_result_file
//
//  Purpose:   Class for analysis result.
//----------------------------------------------------------------------------

class analysis_result_file : public analysis_viewer
{
  public:

    analysis_result_file(ostream& output)
      : analysis_viewer(output)
    { }

    virtual void print_results(const primary_analysis_results&);
    virtual void print_footer(const primary_analysis_results&);
};

//----------------------------------------------------------------------------
//  Class:     analysis_detailed_file
//
//  Purpose:   Class for analysis detailed result.
//
//----------------------------------------------------------------------------

class analysis_detailed_file : public analysis_viewer
{
  public:

    analysis_detailed_file(ostream& output, bool debug = false)
      : analysis_viewer(output), my_debug(debug)
    { }

    virtual void print_results(const primary_analysis_results&);
    virtual void print_footer(const primary_analysis_results&);

  protected:
  
    bool my_debug;
};

//----------------------------------------------------------------------------
//  Class:     analysis_intermediate_file
//
//  Purpose:   Class for analysis intermediate result.
//
//----------------------------------------------------------------------------

class analysis_intermediate_file : public analysis_viewer
{
  public:

    analysis_intermediate_file(ostream& output)
      : analysis_viewer(output)
    { }

    virtual void print_results(const primary_analysis_results&);
    virtual void print_footer(const primary_analysis_results&);
};

//#include "segreg/analysis_out.ipp"

} // end of namespace SEGREG
} // end of namespace SAGE

#endif
