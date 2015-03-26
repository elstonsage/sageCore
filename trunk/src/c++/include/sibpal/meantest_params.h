#ifndef SIBPAL_MEAN_TEST_PARAMS_H
#define SIBPAL_MEAN_TEST_PARAMS_H

//=============================================================================
// File:    mean_test_params.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Took it out from sibmeantest.h                          yjs  Jun.05
//
// Notes:   This file contains definition for following data structures.
//            class mean_estimate
//            struct marker_parameter
//            struct trait_parameter
//            class test_parameters
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/definitions.h"

namespace SAGE   {
namespace SIBPAL {

class mean_estimate
{
  public:

    mean_estimate();

    double        f(size_t n)                   const;
    double        residual_variance()           const;
    double        standard_error()              const;
    double        residual_variance(size_t n)   const;
    double        standard_error(size_t n)      const;

    void          set_pi(double p);
    void          set_f1(double f);
    void          set_cov(const matrix& c);

    void          set_w(double w);

    double        pi()   const;
    double        f1()   const;
    const matrix& cov()  const;

    double        w()    const;

    void          clear();

  private:

    double my_pi;
    double my_f1;
    matrix my_ss;

    double my_w;
};

struct marker_parameter
{
  public:

    marker_parameter();
    marker_parameter(size_t m);

    void clear();

    string name(const relative_pairs& pairs) const;

    bool operator==(const marker_parameter& p) const;
    bool operator!=(const marker_parameter& p) const;
    bool operator< (const marker_parameter& p) const;

    double t_value(double mean)                    const;
    double t_value(size_t n, double mean)          const;
    double p_value(long df, double mean)           const;
    double p_value(size_t n, long df, double mean) const;

    int                 affecteds;
    size_t              pair_count;
    size_t              marker;
    mean_estimate       estimate;

};

struct trait_parameter
{
  public:

    trait_parameter();
    trait_parameter(size_t t, int a = -1);

    string name(const relative_pairs& pairs) const;

    size_t affected_count() const;
    size_t affected_types() const;

    void set_affection(int a);
    void clear();

    bool is_set() const;

    bool operator==(const trait_parameter& p) const;
    bool operator!=(const trait_parameter& p) const;
    bool operator< (const trait_parameter& p) const;

    size_t trait;
    bool   affection[3];
};

class meantest_parameters
{
  public:

    friend class SibMeanTest;

    typedef std::vector<marker_parameter>             marker_parameter_vector;
    typedef marker_parameter_vector::iterator         marker_parameter_iterator;
    typedef marker_parameter_vector::const_iterator   marker_parameter_const_iterator;

    typedef std::vector<trait_parameter>              trait_parameter_vector;
    typedef trait_parameter_vector::iterator          trait_parameter_iterator;
    typedef trait_parameter_vector::const_iterator    trait_parameter_const_iterator;

    meantest_parameters();

    const marker_parameter& beta(size_t i)          const;
    const marker_parameter& operator[](size_t i)    const;
    const trait_parameter&  subsets(size_t i) const;

    void clear_marker_parameters();

    std::pair<marker_parameter_const_iterator, bool> add_marker_parameter(const marker_parameter& p);
    std::pair<marker_parameter_const_iterator, bool> set_marker_parameter(const marker_parameter& p);
    std::pair<marker_parameter_const_iterator, bool> add_marker(size_t m);
    std::pair<marker_parameter_const_iterator, bool> set_marker(size_t m);

    std::pair<trait_parameter_const_iterator, bool>  add_subset(const trait_parameter& t);
    std::pair<trait_parameter_const_iterator, bool>  add_subset(size_t t);
    std::pair<trait_parameter_const_iterator, bool>  add_subset(size_t t, int a);

    void                    set_affection_status(int a);
    void                    set_trait(size_t t);
    void                    set_trait(size_t t, int a);
    void                    validate();
    void                    invalidate();

    bool                    valid()               const;
    int                     affection_status()    const;
    const trait_parameter&  trait()               const;
    bool                    trait_set()           const;
    size_t                  marker_count()        const;
    size_t                  subset_count()        const;

    void                    set_pvalues_scientific_notation(bool a);
    bool                    get_pvalues_scientific_notation() const;

    void                    set_export_output(bool a);
    bool                    get_export_output() const;

    void                    set_use_full_sibs(bool s);
    bool                    get_use_full_sibs() const;

    void                    set_use_half_sibs(bool s);
    bool                    get_use_half_sibs() const;

    void                    set_w(double w);
    double                  get_w() const;

  protected:

    marker_parameter& beta(size_t i);
    marker_parameter& operator[](size_t i);
    trait_parameter&  subsets(size_t i);

  private:

    bool                      my_valid_regression;
    trait_parameter           my_trait;
    marker_parameter_vector   my_betas;
    trait_parameter_vector    my_subsets;

    bool                      my_pval_s_notation;
    bool                      my_export_output;

    bool                      my_use_full_sibs;
    bool                      my_use_half_sibs;

    double                    my_w;
};

#include "sibpal/meantest_params.ipp"

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
