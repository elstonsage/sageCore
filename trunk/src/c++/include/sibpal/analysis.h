#ifndef SIBPAL_ANALYSIS_H
#define SIBPAL_ANALYSIS_H

//=============================================================================
// File:    sibpal_analysis.h
//
// Author:  Yeunjoo Song
//
// History: Version 0.0 Initial implementation                          Jun. 05
//
// Notes:   This file contains definition for following data structures.
//            class SibpalAnalysis
//
// Copyright (c) 2005   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/meantest_out.h"
#include "sibpal/regress_out.h"

namespace SAGE   {
namespace SIBPAL {

class sibpal_analysis
{
  public:

    sibpal_analysis(cerrorstream& err = sage_cerr);

    void run_mean_test          (relative_pairs& pairs, meantest_parser&   mparser, ostream &out, ostream &csv);
    void run_regression_analysis(relative_pairs& pairs, regression_parser& rparser, ostream &out, ostream &det, ostream &csv);

  protected:

    // Build models
    //
    void build_models(const regression_parser& rparser);

    void add_markers(regression_model& a_model, const regression_parser& r_parser);
    void add_covariates(regression_model& a_model, const regression_parser& r_parser);
    void add_interactions(regression_model& a_model, const regression_parser& r_parser);
    void add_batch_interactions(regression_model& a_model, const regression_parser& r_parser, const marker_type& mar, bool d);
    void add_other_options(regression_model& a_model, const regression_parser& r_parser);

    void add_a_trait(regression_model& a_model, const trait_type& trt);
    void add_a_marker(regression_model& a_model, const marker_type& mar);
    void add_a_covariate(regression_model& a_model, const covariate_type& cov);

    bool check_main_effects(const regression_model&      reg,  
                            const independent_variable&  interaction) const;

    bool check_main_marker_effect(const regression_model& reg, const marker_type& mar) const;
    bool check_main_covariate_effect(const regression_model& reg, const covariate_type& cov) const;

    void set_output_options(const regression_parser& r_parser);

    bool is_valid_model(const regression_model& reg) const;

    const relative_pairs*     my_pairs;

    vector<regression_model>  my_analysis_models;
    regression_model          my_current_model;

    cerrorstream&             errors;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
