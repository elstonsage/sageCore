#ifndef RELPAL_TWO_LEVEL_TEST_H
#define RELPAL_TWO_LEVEL_TEST_H

//=============================================================================
// File:    two_level_test.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                      Mar. 07
//
// Notes:   This file contains definition for following data structures.
//            class two_level_base
//            class two_level_regression
//            class two_level_score_test
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/two_level_calculator.h"
#include "relpal/semi_definite.h"
#include "relpal/covariance_calculator.h"

namespace SAGE   {
namespace RELPAL {

class two_level_base
{
  public:

    two_level_base(cerrorstream& err = sage_cerr);

    void set_pairs(const relative_pairs* r);
    void set_regression_type(const regression_type& r);
    void set_model(const regression_model& m);
    void set_data(const analysis_data_type& d);

    void set_data_out(ostream* o);
    void set_debug_out(ostream* o);

    const relative_pairs*      get_pairs()           const;
    const regression_type&     get_regression_type() const;
    const regression_model&    get_model()           const;
    const relpal_least_square& get_gls_1()           const;

  protected:

    void    compute_covariate_means();

    void    get_ind_x_y(const vector< mem_pointer >& sp, matrix& x, matrix& y) const;

    bool    build_ind_residuals(const vector<matrix>& ys,
                                const vector<matrix>& xs,
                                      vector<matrix>& resids);

    double  get_trait_value(size_t t, mem_pointer mem) const;

    double  get_ind_parameter_value(size_t t, size_t p, mem_pointer mem) const;

    double  get_ped_parameter_value(size_t t1, size_t t2,
                                    size_t p, size_t pair_index,
                                    mem_pointer mem1,
                                    mem_pointer mem2) const;

    double  get_marker_sharing_value(size_t p, size_t pair_index,
                                     mem_pointer mem1,
                                     mem_pointer mem2) const;

    double  get_covariate_value(size_t p, size_t pair_index,
                                mem_pointer mem1,
                                mem_pointer mem2) const;

    double  get_interaction_value(size_t p, size_t pair_index,
                                  mem_pointer mem1,
                                  mem_pointer mem2) const;

    void    get_V(const matrix& w, matrix& v) const;

    void    print_variances(ostream &out) const;

    ostream& data_out();
    ostream& debug_out();

    // Members
    //
    const relative_pairs*     my_pairs;

    regression_type           my_reg_type;

    regression_model          my_model;

    analysis_data_type        my_data;

    relpal_least_square       my_gls_1;

    vector<double>            my_variances;

    ostream*                  my_data_out;
    ostream*                  my_debug_out;

    cerrorstream&             errors;
};

class two_level_regression : public two_level_base
{
  public:

    two_level_regression(cerrorstream& err = sage_cerr);

    bool    do_two_level_regression();
    bool    get_null_weight(vector< pair<trimatrix, trimatrix> >& null_ws,
                            vector<matrix>& resids);

    const relpal_least_square& get_gls_2() const;

  protected:

    void    do_individual_level(vector<matrix>& ys,
                                vector<matrix>& xs,
                                vector< pair<trimatrix, trimatrix> >& vs);

    void    do_pedigree_level(const vector< pair<trimatrix, trimatrix> >& vs,
                              const vector<matrix>&    resids);

    bool    check_positive_semidefinite();

    double  find_max_diff(const vector<double>& prev_var, const vector<double>& curr_var) const;

    // Members
    //
    relpal_least_square     my_gls_2;

};

class two_level_score_test : public two_level_base
{
  public:

    two_level_score_test(cerrorstream& err = sage_cerr);

    void    compute_null_weight(const regression_model& null_model,
                                regression_result&      re);

    bool    do_two_level_score_test(score_test_result& re);

    const relpal_score&   get_score_2() const;
    const vector<matrix>& get_residuals() const;

    double  get_correction(size_t type = 1) const;
    double  get_emp_pvalue(size_t type = 1) const;
    size_t  get_rep_count(size_t type = 1)  const;

    bool    has_valid_null()     const;
    bool    has_reliable_score() const;

  protected:

    bool    do_pedigree_level();

    double  compute_correction(const matrix& U_star, const matrix& sigma_star_i);

    double  maximize_correction_1(semi_definite& correction__fn, double lb = NE_INF);
    double  maximize_correction_2(semi_definite& correction__fn, double lb = NE_INF);

    double  compute_empirical_p_value(size_t& rep, size_t test_df, double T, double correction,
                                      const matrix& sigma_star_i,  const matrix& sigma_star_sqrt);
    double  get_random_normal();

    void    correction_tests();

    void    get_B_matrix(size_t n, const trimatrix& b, const matrix& S, matrix& B);

    void    compute_covariance_matrix();
    void    compute_covariance_matrix_given_m();

    // Members
    //
    vector< pair<trimatrix, trimatrix> >  my_null_weights;

    vector<matrix>     my_residuals;

    vector<trimatrix>  my_null_IBD_covariances;
    vector<trimatrix>  my_IBD_covariances;

    relpal_score       my_score_2;

    double             my_correction_na;
    double             my_correction_sw;
    double             my_correction_al;
    double             my_correction_ib;

    double             my_emp_pvalue_na;
    double             my_emp_pvalue_sw;
    double             my_emp_pvalue_al;
    double             my_emp_pvalue_ib;

    size_t             my_rep_count_na;
    size_t             my_rep_count_sw;
    size_t             my_rep_count_al;
    size_t             my_rep_count_ib;

    bool               my_valid_null;
    bool               my_reliable_score;

    MersenneTwister    my_random_generator;
};

#include "relpal/two_level_test.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
