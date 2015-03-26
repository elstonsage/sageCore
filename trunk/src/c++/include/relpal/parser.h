#ifndef RELPAL_PARSER_H
#define RELPAL_PARSER_H

//=============================================================================
// File:    relpal_parser.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                  yjs Mar. 07
//
// Notes:   This file contains definition for following data structures.
//
//            class relpal_parser : public APP::BasicParser
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "app/aparser.h"
#include "relpal/model.h"

namespace SAGE   {
namespace RELPAL {

class relpal_parser : public APP::BasicParser
{
  public:

    relpal_parser(const relative_pairs& relpairs, cerrorstream& err = SAGE::sage_cerr);

    void parse_symbols(const SymbolTable* syms);
    void parse_parameter(const LSFBase* param);
    void parse_test_parameter_section(const LSFBase* params);
    void parse_test_parameter(const LSFBase* param);

    void clear();

    void set_default_use_members(bool mibd);

    bool do_ind_batch_test()                                         const;
    bool do_first_level_test()                                       const;

    const regression_type&               get_regression_type()       const;

    const vector<dependent_variable>&    get_traits()                const;

    const vector<covariate_type>&        get_ind_covariates()        const;
    const vector<covariate_type>&        get_ind_batch_covariates()  const;
    const vector<covariate_type>&        get_ped_null_covariates()   const;
    const vector<covariate_type>&        get_ped_test_covariates()   const;

    const vector<marker_type>&           get_null_markers()          const;
    const vector<marker_type>&           get_test_markers()          const;

    const vector<interaction_type>&      get_ind_interactions()      const;
    const vector<interaction_type>&      get_ped_null_interactions() const;
    const vector<interaction_type>&      get_ped_test_interactions() const;

    const data_options&                  get_data_options()          const;
    const analysis_options&              get_analysis_options()      const;
    const pvalue_options&                get_pvalue_options()        const;
    const output_options&                get_output_options()        const;

    void dump_parser(ostream &out)                                   const;

  private:

    void parse_trait(const LSFBase* param);
    void parse_model_type(const LSFBase* param);

    bool parse_effect_parameters(const LSFBase* params, bool first);

    void parse_marker(const LSFBase* param);
    void parse_covariate(const LSFBase* param, bool first);
    void parse_interaction(const LSFBase* param, bool first);

    void parse_marker_type(     const LSFBase* param, marker_type& mar);
    void parse_covariate_type(  const LSFBase* param, covariate_type& cov);
    void parse_interaction_type(const LSFBase* param, interaction_type& inter, bool first);

    void parse_data_options(const LSFBase* params);

    void parse_subset(const LSFBase* param);
    void parse_use_pairs(const LSFBase* param);
    void parse_use_members(const LSFBase* param);

    void parse_output_options(const LSFBase* params);

    void parse_pvalue_options(const LSFBase* params);

    void make_MM_interaction_type(interaction_type& inter, const marker_type&    mar1,
                                                           const marker_type&    mar2) const;

    void make_MC_interaction_type(interaction_type& inter, const marker_type&    cov,
                                                           const covariate_type& mar)  const;

    void process_batch_MM_interaction(vector<interaction_type>& ped_test_inters)      const;
    void process_batch_MC_interaction(vector<interaction_type>& ped_test_inters)      const;

    void update_null_markers();
    void update_null_covariates();

    void add_test_markers();
    void switch_to_single_marker_model();

    bool is_in_first_level_model(size_t batch_cov) const;

    // Data members

    const relative_pairs&         my_pairs;

    regression_type               my_reg_type;

    vector<dependent_variable>    my_traits;

    vector<covariate_type>        my_ind_covariates;
    vector<covariate_type>        my_ind_batch_covariates;
    vector<covariate_type>        my_ped_null_covariates;
    vector<covariate_type>        my_ped_test_covariates;

    vector<marker_type>           my_null_markers;
    vector<marker_type>           my_test_markers;

    vector<interaction_type>      my_ind_interactions;
    vector<interaction_type>      my_ped_null_interactions;
    vector<interaction_type>      my_ped_test_interactions;

    data_options                  my_data_opt;

    analysis_options              my_analysis_opt;

    pvalue_options                my_pvalue_opt;

    output_options                my_output_opt;

    bool                          my_ind_batch_test;
};

#include "relpal/parser.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
