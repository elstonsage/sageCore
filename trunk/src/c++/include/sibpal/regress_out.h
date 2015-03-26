#ifndef SIBPAL_REGRESSION_OUTPUT_H
#define SIBPAL_REGRESSION_OUTPUT_H

//=============================================================================
// File:    regress_out.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//
//            class trait_regression_viewer
//            class trait_regression_textfile : public trait_regression_viewer
//            class trait_regression_csvfile  : public trait_regression_viewer
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress_variants.h"

namespace SAGE   {
namespace SIBPAL {

class regression_outfile
{
  public:

    regression_outfile(ostream& sum, ostream& det, ostream& exp, cerrorstream& err = sage_cerr);

    ~regression_outfile() {}

    ostream&  sum_output_stream()        { return my_sum_output; }
    ostream&  det_output_stream()        { return my_det_output; }
    ostream&  exp_output_stream()        { return my_exp_output; }

    const TraitRegression*  get_analysis() const   { return my_analysis; }

    void print_header(const TraitRegression& reg);
    void print_results(size_t test_num);
    void print_footer();

    void print_instability_note();

  protected:

    // Header
    //
    void print_common_header(ostream& o, bool det = false);
    void print_exp_header();

    // Common
    void print_analysis_title   (ostream& o);
    void print_trait_name       (ostream& o);
    void print_dependent_name   (ostream& o);
    void print_subset_info      (ostream& o);
    void print_other_option_info(ostream& o);
    void print_legend           (ostream& o);
    void print_double_line      (ostream& o, bool det = false);

    // Results
    //
    void print_sum_results(size_t test_num);
    void print_det_results(size_t test_num);
    void print_exp_results(size_t test_num);

    // Sum
    void print_sum_table_heading();
    void print_sum_single_line();
    void print_sum_result(size_t test_number, size_t i);
    void print_joint_F_test();

    // Det
    void print_model_info();
    void print_sample_info();
    void print_trait_info();
    void print_dependent_info(size_t test_num);
    void print_regression_info();

    // Exp
    void print_exp_result(size_t test_number, size_t i);

    const TraitRegression*    my_analysis;

    bool                      my_detailed_out;
    bool                      my_export_out;

    ostream&                  my_sum_output;
    ostream&                  my_det_output;
    ostream&                  my_exp_output;

    size_t                    my_pre_width;
    size_t                    my_param_name_max_size;
    size_t                    my_param_effect_name_max_size;
    size_t                    my_x_type_count;

    cerrorstream&             errors;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
