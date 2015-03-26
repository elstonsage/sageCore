#ifndef SIBPAL_REGRESSION_PARAMS_H
#define SIBPAL_REGRESSION_PARAMS_H

//=============================================================================
// File:    regress_params.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Took it out from sibregress.h                           yjs  Jun.05
//
// Notes:   This file contains definition for following data structures.
//            struct dependent_variable
//            struct covariate_type
//            struct marker_type
//            struct independent_variable
//            class regression_parameters
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/definitions.h"

namespace SAGE {
namespace SIBPAL {

struct trait_type
{
  explicit trait_type(size_t t = (size_t)-1, bool bin = false,
                      double m = QNAN, double p = QNAN);

  string name(const relative_pairs& pairs) const;

  bool operator==(const trait_type& t) const;
  bool operator!=(const trait_type& t) const;

  bool           binary;
  bool           valid;

  size_t         trait_index;

  double         fixed_mean;
  double         fixed_male_mean;
  double         fixed_female_mean;

  double         fixed_correlation;       // User entered trait correlation

  SampleInfo     trait_used_sibs_info;    // member-wise trait info for used sibs

  SampleInfo     trait_all_sibs_info;     // member-wise trait info for all sibs
  SampleInfo     trait_male_sibs_info;
  SampleInfo     trait_female_sibs_info;

  double         var_b;
  double         var_r;

  // Correlations of trait
  //
  double         correlation;             // intra sibship trait correlation
  double         fisher_correlation;      // pair-wise trait correlation
  double         fisher_fsib_correlation;
  double         fisher_hsib_correlation;
  double         sum_diff_correlation;    // correlation between squared sum & diff
};

struct dependent_variable : public trait_type
{
  explicit dependent_variable(size_t t = (size_t)-1, bool bin = false,
                              double m = QNAN, double p = QNAN);

  double get_mean() const;

  bool operator==(const dependent_variable& t) const;
  bool operator!=(const dependent_variable& t) const;

  // Correlations of dependent_variable
  //
  vector<double> p_correlation;
  vector<double> p_fsib_correlation;
  vector<double> p_hsib_correlation;

  vector<double> p_empirical_correlation;
  vector<double> p_fsib_empirical_correlation;
  vector<double> p_hsib_empirical_correlation;
  vector<double> p_fh_empirical_correlation;

  // Parameter correlation for dependent_variable W4
  //  correlation between residuals from sum & diff regression
  vector<double> p_sum_diff_correlation;
  vector<double> p_fsib_sum_diff_correlation;
  vector<double> p_hsib_sum_diff_correlation;
  vector<double> p_fh_sum_diff_correlation;

  // Dependent variable info
  //
  SampleInfo     info;                 // pair-wise dependent variable info
  SampleInfo     mm_pair_info;         // mm pair dependent variable info
  SampleInfo     mf_pair_info;         // mf pair dependent variable info
  SampleInfo     ff_pair_info;         // ff pair dependent variable info
};

struct covariate_type
{
  enum cov_op { none, sum, diff, prod, both, avg, single, pair };

  covariate_type();

  explicit covariate_type(size_t c, cov_op op, double powr = 1.0, double fixed_mean = QNAN);

  string name(const relative_pairs& pairs) const;
  string effect_name() const;
  string short_effect_name() const;

  bool operator==(const covariate_type& c) const;
  bool operator!=(const covariate_type& c) const;
  bool operator< (const covariate_type& c) const;

  size_t         covariate_index;
  cov_op         operation;

  double         power;
  double         fixed_mean;

  SampleInfo     info;
};

struct marker_type
{
  enum effect_type { TOTAL, ADDITIVE, DOMINANCE, OTHER };

  explicit marker_type(size_t m = (size_t)-1, effect_type e = TOTAL, bool x = false);

  string name(const relative_pairs& pairs) const;
  string effect_name() const;
  string short_effect_name() const;

  bool operator==(const marker_type& m) const;
  bool operator!=(const marker_type& m) const;
  bool operator< (const marker_type& m) const;

  size_t         marker_index;
  effect_type    effect;

  bool           x_linked;
};

struct independent_variable
{
  enum indep_var_type { MARKER=0, COVARIATE, INTERACTION, INVALID };

  independent_variable();

  string    name(const relative_pairs& pairs)             const;
  string    effect_name()                                 const;
  string    short_effect_name()                           const;

  string    covariate_name(const relative_pairs& pairs)   const;
  string    covariate_effect_name(bool shortname=false)   const;
  string    marker_name(const relative_pairs& pairs)      const;
  string    marker_effect_name(bool shortname=false)      const;

  bool operator==(const independent_variable& c) const;
  bool operator!=(const independent_variable& c) const;
  bool operator< (const independent_variable& c) const;

  void clear();

  bool                    valid;

  vector<marker_type>     markers;
  vector<covariate_type>  covariates;
  indep_var_type          type;

  SampleInfo              info;
  SampleInfo              mm_pair_info;
  SampleInfo              mf_pair_info;
  SampleInfo              ff_pair_info;
};

struct data_options
{
  data_options();

  string dump(string pre_space) const;

  use_pair_type  use_pairs;

  bool           skip_uninformative_pairs;
  bool           use_mm_pair;
  bool           use_mf_pair;
  bool           use_ff_pair;
};

struct analysis_options
{
  analysis_options();

  string dump(string pre_space) const;

  bool           standardize_parameter;
  bool           robust_variance;
  bool           leverage_adjustment;
  bool           identity_weight;
  bool           pool_correlation;
  bool           blup_mean;
  bool           sibship_mean;
  size_t         sibship_mean_threshold;
  size_t         max_iterations;
  double         w1;
};

struct empirical_pvalue_options
{
  empirical_pvalue_options();

  string dump(string pre_space) const;

  bool           is_on;

  int            seed;

  double         width;
  double         confidence;
  double         threshold;

  size_t         replicates;
  size_t         min_replicates;
  size_t         max_replicates;

};

typedef empirical_pvalue_options pvalue_options;

struct output_options
{
  output_options();

  string dump(string pre_space) const;

  bool           wide_out;
  bool           detailed_out;
  bool           export_out;

  bool           print_design_matrix;
  bool           print_correl_matrix;
  bool           pvalues_scientific_notation;

  size_t         design_matrix_rows;
  string         design_matrix_file;
  string         correl_matrix_file;

  size_t         param_name_max_size;
  size_t         param_effect_name_max_size;

  bool           dump_data;
  bool           dump_debug;
  string         dump_data_file;
};

class regression_model
{
  public:

    friend class sibpal_analysis;
    friend class TraitRegression;

    typedef vector<independent_variable>       parameter_vector;
    typedef parameter_vector::iterator         parameter_iterator;
    typedef parameter_vector::const_iterator   parameter_const_iterator;

    typedef vector<trait_type>                 trait_vector;
    typedef trait_vector::iterator             trait_iterator;
    typedef trait_vector::const_iterator       trait_const_iterator;

    typedef pair<trait_iterator, bool>         tib_value;
    typedef pair<parameter_iterator, bool>     pib_value;

    regression_model();

    void    clear();
    void    invalidate();

    bool    valid()                                             const;
    bool    is_x_linked()                                       const;

    size_t                         get_parameter_count()        const;
    size_t                         get_subset_count()           const;

    dependent_variable&            get_trait();
    independent_variable&          get_parameter(size_t i);
    trait_type&                    get_subset(size_t t);

    const dependent_variable&      get_trait()                  const;
    const independent_variable&    get_parameter(size_t i)      const;
    const trait_type&              get_subset(size_t t)         const;

    regression_type                get_regression_type()        const;
    string                         get_regression_method_name() const;

    const data_options&            get_data_options()           const;
    const analysis_options&        get_analysis_options()       const;
    const pvalue_options&          get_pvalue_options()         const;
    const output_options&          get_output_options()         const;

    parameter_iterator             parameter_begin();
    parameter_iterator             parameter_end();

    parameter_const_iterator       parameter_begin()            const;
    parameter_const_iterator       parameter_end()              const;

    void      set_trait(size_t t, bool b);
    void      set_trait(const dependent_variable& t);

    void      set_regression_type(regression_type r);
    void      set_regression_method_name(string n);

    void      set_data_options(const data_options& op);
    void      set_analysis_options(const analysis_options& op);
    void      set_pvalue_options(const pvalue_options& op);
    void      set_output_options(const output_options& op);

    void      clear_subsets();

    tib_value add_subset(const trait_type& t);
    tib_value add_subset(size_t t);
    tib_value set_subset(size_t t);

    void      clear_parameters();
    void      clear_markers();
    void      clear_invalid_parameters();

    pib_value add_parameter(const independent_variable& p);
    pib_value add_marker(size_t m, marker_type::effect_type e = marker_type::TOTAL, bool x_linked = false);
    pib_value add_covariate(size_t c, covariate_type::cov_op op, double power = 1.0, double fixed_mean = QNAN);

    void      dump_model(const relative_pairs& relpairs, ostream &out) const;

  protected:

    void      validate();

  private:

    bool                       my_valid;
    bool                       my_x_linked;

    dependent_variable         my_trait;
    trait_vector               my_subsets;
    parameter_vector           my_parameters;

    regression_type            my_reg_type;
    string                     my_reg_method_name;

    data_options               my_data_opt;
    analysis_options           my_analysis_opt;
    pvalue_options             my_pvalue_opt;
    output_options             my_output_opt;
};

//typedef regression_model regression_parameters;

#include "sibpal/regress_params.ipp"

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
