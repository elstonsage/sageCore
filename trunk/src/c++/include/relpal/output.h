#ifndef RELPAL_OUTPUT_H
#define RELPAL_OUTPUT_H

//=============================================================================
// File:    output.h
//
// Author:  Yeunjoo Song
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structures.
//
//            class relpal_viewer
//            class relpal_textfile : public relpal_viewer
//            class relpal_csvfile  : public relpal_viewer
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/analysis.h"

namespace SAGE   {
namespace RELPAL {

 
class relpal_outfile
{
  public:

    relpal_outfile(ostream& sum, ostream& det, ostream& exp, const relpal_analysis& ra);

    virtual ~relpal_outfile() {}

    ostream&                sum_output_stream()        { return my_sum_output; }
    ostream&                det_output_stream()        { return my_det_output; }
    ostream&                exp_output_stream()        { return my_exp_output; }

    const relpal_analysis&  get_analysis() const   { return my_analysis;   }

    void print_header(bool first_level);
    void print_results(size_t test_num, bool first_level, bool valid = true, bool reliable = true);
    void print_footer(bool first_level);

  protected:

    // Header
    //
    void print_sum_header(bool first_level);
    void print_det_header(bool first_level);
    void print_exp_header();

    // Sum & Det common heading

    void print_analysis_title(ostream& o, bool first_level, string s);

    void print_trait_summary(ostream& o);
    void print_subset_summary(ostream& o);
    void print_pval_option(ostream& o);
    void print_legend(ostream& o);

    void print_table_heading_line1(ostream& o, bool print_cM = true, bool print_count = true);
    void print_table_heading_line2(ostream& o, bool sum, bool print_cM = true, bool print_count = true);

    void print_double_line(ostream& o, bool sum = true);
    void print_single_line(ostream& o, bool print_cM = true, bool print_count = true);

    void print_reg_table_heading  (ostream& o, bool sum, bool print_cM = true, bool print_count = true);
    void print_score_table_heading(ostream& o, bool sum, bool print_cM = true, bool print_count = true);

    void print_reg_single_line  (ostream& o, bool sum, bool print_cM = true, bool print_count = true);
    void print_score_single_line(ostream& o, bool print_cM = true, bool print_count = true);

    // Results
    //
    void print_first_level_results (size_t test_num, bool valid = true);
    void print_second_level_results(size_t test_num, bool valid = true, bool reliable = true);

    // Sum & Det common result

    void print_reg_result_line(ostream& o,
                               string name,  size_t max_name_size,
                               string cM,    string count,
                               double T_val, double p_val);

    void print_reg_result_block(ostream& o,
                                string name,        size_t max_name_size,
                                string cM,          string count,
                                const matrix& beta, const matrix& var,
                                double T_val, double p_val);

    void print_score_result_line(ostream& o,
                                 string name, size_t max_name_size,
                                 string cM,   string count,
                                 double Tun,  double correct,
                                 double empp, size_t rep);

    void print_fail_line(ostream& o,
                         string name, size_t max_name_size,
                         string cM,   string count);

    void get_reg_param_results(size_t i, size_t t,
                               const matrix& beta, const matrix& var,
                               matrix& param_beta, matrix& param_var);

    pair<double, double> compute_multivariate_T_value(const matrix& beta, const matrix& var);
    pair<double, double> compute_one_sided_T_value(const matrix& beta, const matrix& var);

    // Det specific result

    void print_det_first_level_results (size_t test_num);
    void print_det_second_level_results(size_t test_num, bool valid, bool reliable);

    void print_model_heading();
    void print_first_model_info();
    void print_second_model_info();

    void print_sample_info();

    void print_reg_result_heading(string h);
    void print_score_result_heading();

    void print_reg_estimates(const regression_model& model, const regression_result& reg_result);

    string get_trait_names() const;
    string get_independent_variable_names(const regression_model& model) const;

    // Members
    //
    const relpal_analysis&    my_analysis;

    bool                      detailed_out;
    bool                      export_out;

    vector<string>            my_var_opt;

    ostream&                  my_sum_output;
    ostream&                  my_det_output;
    ostream&                  my_exp_output;
};

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
