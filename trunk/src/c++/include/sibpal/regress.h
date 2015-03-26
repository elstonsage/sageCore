#ifndef SIBPAL_REGRESSION_H
#define SIBPAL_REGRESSION_H

//=============================================================================
// File:    regress.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Half sib added.                                       yjs Mar. 2006
//          X-linkage added.                                      yjs Jun. 2006
//
// Notes:   This file contains definition for following data structures.
//
//            class RegressionVariant
//            class TraitRegression
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/parser.h"
#include "sibpal/regress_result.h"

namespace SAGE   {
namespace SIBPAL {

class RegressionVariant
{
  public:

    typedef dependent_variable             trait_parameter;

    RegressionVariant(string, cerrorstream& err = sage_cerr);

    virtual ~RegressionVariant();

    string          name()        const;

    virtual string  description() const = 0;

    virtual void    regress() = 0;
    virtual void    build() = 0;
    virtual void    update() = 0;
    virtual void    simulate() = 0;

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& t,
                                 bool                   use_empirical_correlations = false,
                                 bool                   center = false) const = 0;

    virtual void    weight_matrix(const sib_cluster&     ship,
                                  matrix&                W,
                                  const trait_parameter& t,
                                  bool                   use_empirical_correlations,
                                  weight_status_type&    status,
                                  double                 fsib_var,
                                  double                 hsib_var) const = 0;

    virtual bool    apply_weighting() const;

    virtual const w_map&  get_w() const;

    virtual void    clear_w();

  protected:

    void            set_name(const string& n);

    cerrorstream&         errors;
    mutable sib_weights   weights;
    mutable w_map w;

  private:

    string                my_name;
};

class TraitRegression
{
  public:

    friend class WeightedVariance;
    friend class WeightedCorrelatedVariance;
    friend class WeightedCorrelatedVarianceAndTraits;

    typedef dependent_variable                     trait_parameter;
    typedef boost::shared_ptr<RegressionVariant>   regression_variant_pointer;

    TraitRegression(relative_pairs& p, cerrorstream& err = sage_cerr);

    virtual ~TraitRegression();

    void  invalidate();
    void  validate();
    void  reset();

    void  regress();
    void  do_regress();

    void  set_regression_method(RegressionVariant* r);
    void  set_regression_method(regression_variant_pointer r);
    void  disown_regression_method();

    void  set_model(const regression_model& p);
    void  set_sib_clusters(const vector<sib_cluster>& sc);
    void  set_filter(const pair_filter& pf);
    void  set_use_pairs(const pair<bool, bool>& up);
    void  set_pair_counts(const pair<size_t, size_t>& pc);
    void  set_use_x_pairs(const vector<bool>& x);
    void  set_fix_x_pairs(const vector<bool>& x);
    void  set_x_pair_counts(const vector<size_t>& pc);
    void  set_x_types(const vector<pair_type>& pt);
    void  set_reg_results(const regression_results& r);

          relative_pairs&       get_pairs()         const;

          RegressionVariant*    regression_method();
    const RegressionVariant*    regression_method() const;

          regression_model&     get_model();
    const regression_model&     get_model()         const;

    const vector<sib_cluster>&  get_sib_clusters()  const;
    const pair_filter&          get_filter()        const;
          pair<bool, bool>      get_use_pairs()     const;
          pair<size_t, size_t>  get_pair_counts()   const;
    const vector<bool>&         get_use_x_pairs()   const;
    const vector<bool>&         get_fix_x_pairs()   const;
    const vector<size_t>&       get_x_pair_counts() const;
    const vector<pair_type>&    get_x_types()       const;

          regression_results&   get_reg_results();
    const regression_results&   get_reg_results()   const;

          regression_results&   get_reg_intercept_results();
    const regression_results&   get_reg_intercept_results()   const;

    bool   built()                     const;
    bool   valid()                     const;

    size_t pair_count()                const;
    size_t fsib_pair_count()           const;
    size_t hsib_pair_count()           const;
    size_t valid_pair_type_count()     const;
    size_t valid_intercept_count()     const;
    size_t mm_pair_count()             const;
    size_t mf_pair_count()             const;
    size_t ff_pair_count()             const;

    size_t valid_trait_count()         const;
    size_t valid_parameter_count()     const;
    size_t parameter_count()           const;

    size_t get_svd_return_code()       const;

    cerrorstream&  error_stream();
    void           set_error_stream(cerrorstream& e);

    // Get the trait/covariate/ for a sib member.
    double member_trait(const MPED::member_base* m, size_t t) const;

    // Get the trait means for the pair (sex-specific for X-linkage).
    pair<double, double> get_trait_means(  const sib_pair& pair) const;
    pair<double, double> get_trait_means_x(const sib_pair& pair) const;

    // Get the sibship specific means for the pair (sex-specific for X-linkage).
    pair<double, double> get_sibship_means(  const sib_pair& pair)  const;
    pair<double, double> get_sibship_means_x(const sib_pair& pair)  const;

    // Get the BLUP means for the pair (sex-specific for X-linkage).
    pair<double, double> get_blup_means(  const sib_pair& pair)  const;
    pair<double, double> get_blup_means_x(const sib_pair& pair)  const;

    // Compute the correlation between pairs of dependent variables
    void estimate_dependent_variable_correlation();

    // Simulation
    //
    void  simulate();
    void  do_simulate();

  protected:

    //----------------------
    // do_regress functions
    //----------------------

    void make_filter();
    void estimate_trait_sample_mean();

    void regress_univariate(bool do_sim = false, bool check_intercept = false);

    void build_x_types();

    void do_regress_univariate();
    void do_regress_univariate_sub(map<const sib_cluster*, matrix>& As,
                                   map<const sib_cluster*, matrix>& Ws,
                                   map<const sib_cluster*, matrix>& Ys,
                                   map<const sib_cluster*, matrix>& Rs,
                                   size_t m, double fv, double hv);
    void check_intercepts();
    void do_regress_F_test();

    void make_sib_cluster();
    void build();
    void copy_build(const TraitRegression*);

    //-----------------
    // Build functions
    //-----------------

    // Fixup effect types (called from build)
    void normalize_effects();

    bool params_dom_equal(const independent_variable& total_param,
                          const independent_variable& param) const;

    // Compute moments of independent variates
    void estimate_independent_variable_parameters();

    // Compute the correlation between SUM & DIFF
    void estimate_trait_sum_diff_correlation();

    void estimate_trait_correlation();
    void estimate_sibship_blup_mean();

    //---------------------------------
    // do_regress_univariate functions
    //---------------------------------

    void regress_markers(GLS3& gls, bool use_empirical_correlations,
                         map<const sib_cluster*, matrix>& As,
                         map<const sib_cluster*, matrix>& Ws,
                         map<const sib_cluster*, matrix>& Ys,
                         double fv, double hv);

    bool build_residuals(GLS3& gls, map<const sib_cluster*, matrix>& Rs, double* fv, double* hv);

    void weight_matrix(const sib_cluster&     ship,
                       matrix&                W,
                       const trait_parameter& t,
                       bool                   use_empirical_correlations,
                       weight_status_type&    status,
                       double                 fsib_var,
                       double                 hsib_var) const;

    // Compute a vector traits for a given sib_cluster for the given trait parameter
    void trait_vector(const sib_cluster&     ship,
                      matrix&                y,
                      const trait_parameter& t,
                      bool                   use_empirical_correlations = false,
                      bool                   center = false) const;

    // Compute the marker sharing matrix based on the marker parameters
    void design_matrix(const sib_cluster& ship,
                       matrix&            A,
                       bool               standardize) const;

    void design_matrix_x(const sib_cluster& ship,
                         matrix&            A,
                         bool               standardize) const;

    // Decode results
    void compute_robust_variance(GLS3& gls, bool use_empirical_correlations);
    void compute_parameter_result(size_t j, double beta, double ss, size_t df, pair_type pt = MIXED);
    void compute_intercept_result(size_t r, double beta, double ss, pair_type pt = MIXED);

    //-------------------------
    // Design_matrix functions
    //-------------------------

    // Compute the value of a marker sharing parameter for a sib pair
    double marker_sharing(const sib_pair& pair, const marker_type& param, bool permute = true, bool normalize = false) const;
    double marker_sharing(const sib_pair& pair, const independent_variable& param, bool permute = true, bool normalize = false) const;

    // Compute the value of a covariate for a sib pair
    double covariate_value(const sib_pair& pair, const covariate_type& param, bool mean_correct = true) const;
    double covariate_value(const sib_pair& pair, const independent_variable& param, bool mean_correct = true) const;

    double parameter_value(const sib_pair& pair, const independent_variable& param,
                           bool mean_correct_covariates = true,
                           bool permute_markers = true,
                           bool normalize_markers = false) const;

    // Get the traits/covariates/ for each pair in a sib pair -- these values
    // are used throughout and can be randomized.
    pair<double, double> pair_traits(const sib_pair&        pair,
                                     const trait_parameter& param) const;

    //--------------------------
    // Optional print functions
    //--------------------------

    void print_optional_output();

    double compute_b(const matrix& A, const matrix& W, const matrix& Y);
    void dump_data(ostream &out);

    void print_design_matrix_header(ostream& out) const;
    void print_correlation_matrix_header(ostream& out) const;
    void print_sibship_data_matrix(const sib_cluster& ship,
                                   const matrix &w,
                                   ostream& out,
                                   size_t r = 0,
                                   bool rows_restricted = false) const;
    void print_correlation_matrix(ostream& out) const;

    //----------------------
    // Simulation functions
    //----------------------

    bool simulating() const;

    void simulate_univariate();
    void do_simulation_univariate();
    void do_simulation_univariate_replicates(TraitRegression& sim, 
                                             double precision, size_t min_replicates,
                                                               size_t max_replicates);

    // Build simulation data structures (not needed for non-simulation analyses)
    void build_simulation();

    void randomize_marker_data();
    void randomize_marker_data_x();

    // Build during a simulation when we've permuted the design matrix
    void simulating_build();

    // Get a pair based on the current permutation vector
    const sib_pair permutation_pair(const sib_pair& pr) const;

    void copy_simulation_vectors(const TraitRegression& reg);

    //---------
    // Members
    //---------

    relative_pairs&                 pairs;

    regression_variant_pointer      my_regression_method;

    regression_model                my_model;
    pair_filter                     my_filter;

    vector<sib_cluster>             my_sib_clusters;

    size_t                          my_valid_trait_count;
    size_t                          my_valid_parameter_count;

    size_t                          my_fsib_pair_count;
    size_t                          my_hsib_pair_count;

    bool                            my_use_fsib;
    bool                            my_use_hsib;

    bool                            my_built;
    bool                            my_built_trait_all_sibs_info;
    bool                            my_valid_regression;
    bool                            my_simulating;

    vector<bool>                    my_use_x_pair;
    vector<bool>                    my_fix_x_pair;
    vector<size_t>                  my_x_pair_count;
    vector<pair_type>               my_x_types;

    vector<double>                  my_fsib_intercept;
    vector<double>                  my_hsib_intercept;
    vector<vector<double> >         my_x_intercept;
    vector<pair_type>               my_x_intercept_types;

    sib_mean_map                    my_sib_mean_map;
    sib_mean_map                    my_male_sib_mean_map;
    sib_mean_map                    my_female_sib_mean_map;

    regression_results              my_reg_results;
    regression_results              my_reg_results_intercept;

    GLS3                            my_gls;

    simulation_map_type             my_simulation_map;
    group_permutation_vector_type   my_group_permutation_vector;

    MTRandomizer                    randomizer;

    correlation_matrix_map          my_correlation_matrix_map;

    map<const sib_cluster*, matrix> my_final_As;
    map<const sib_cluster*, matrix> my_final_Ws;
    map<const sib_cluster*, matrix> my_final_Ys;
    map<const sib_cluster*, matrix> my_final_Rs;

    cerrorstream&                   errors;

};

boost::shared_ptr<RegressionVariant>
getRegressionVariant(TraitRegression& reg, string name = "");

#include "sibpal/regress.ipp"

} //end of namespace SIBPAL
} //end of namespace SAGE

#endif
