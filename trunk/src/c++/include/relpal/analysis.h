#ifndef RELPAL_ANALYSIS_H
#define RELPAL_ANALYSIS_H

//=============================================================================
// File:    relpal_analysis.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                      Mar. 07
//
// Notes:   This file contains definition for following data structures.
//            class relpal_analysis
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/two_level_test.h"

namespace SAGE   {
namespace RELPAL {

class relpal_analysis
{
  public:

    relpal_analysis(cerrorstream& err = sage_cerr);

    bool run_analysis(relative_pairs& pairs, relpal_parser& parser, bool multipoint_ibd,
                      ostream& sco, ostream& csv, ostream& reg);

    const relative_pairs*   get_pairs()              const;
    const regression_type&  get_regression_type()    const;
    const regression_model& get_current_model()      const;
    const regression_model& get_current_null_model() const;
    const analysis_result&  get_current_result()     const;

    size_t                  get_member_count()       const;
    size_t                  get_pair_count()         const;
    size_t                  get_fsib_pair_count()    const;
    size_t                  get_hsib_pair_count()    const;

    const data_options&     get_data_options()       const;
    const output_options&   get_output_options()     const;

  protected:

    void build_by_pedigree(relative_pairs& r_pairs, const regression_model& mo, ostream& out);

    pair_filter make_filter(const regression_model& r_model);

    bool do_two_level_regression(ostream* out);

    void dump_by_pedigree(ostream &out) const;
    void dump_residuals(const vector<matrix>& res, ostream &out) const;

    // Build models
    //
    void build_models(const relpal_parser& parser);

    void add_traits(regression_model& a_model, const relpal_parser& parser);
    void add_ind_level_effects(regression_model& a_model, const relpal_parser& parser);
    void add_ped_level_null_effects(regression_model& a_model, const relpal_parser& parser, bool addm);
    void add_other_options(regression_model& a_model, const relpal_parser& parser);

    void add_ind_covariates  (regression_model& mo, const relpal_parser& parser);
    void add_a_ind_covariate (regression_model& mo, const covariate_type& cov);
    void add_ind_interactions(regression_model& mo, const relpal_parser& parser);

    void add_a_marker         (regression_model& mo, const marker_type& mar, bool test_variable);
    void add_a_ped_covariate  (regression_model& mo, const covariate_type& cov, bool test_variable);
    void add_a_ped_interaction(regression_model& mo, const interaction_type& inter, bool test_variable);

    void set_data_options(const relpal_parser& parser);
    void set_output_options(const relpal_parser& parser);

    // Members
    //
    const relative_pairs*     my_pairs;

    regression_type           my_reg_type;

    vector<regression_model>  my_analysis_models;
    vector<regression_model>  my_null_models;

    regression_model          my_filter_model;
    regression_model          my_current_model;
    regression_model          my_current_null_model;

    analysis_data_type        my_current_data;

    analysis_result           my_current_result;

    size_t                    my_member_count;
    size_t                    my_pair_count;
    size_t                    my_fsib_pair_count;
    size_t                    my_hsib_pair_count;

    data_options              my_data_opt;
    output_options            my_output_opt;

    cerrorstream&             errors;
};

#include "relpal/analysis.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif

