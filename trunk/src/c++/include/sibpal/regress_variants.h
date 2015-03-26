#ifndef SIBPAL_REGRESSION_VARIANTS_H
#define SIBPAL_REGRESSION_VARIANTS_H

//=============================================================================
// File:    regress_variants.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Took it out from sibregress.h & regressvariants.cpp.    yjs  Jun.05
//
// Notes:   This file contains definition for following data structures.
//
//            class VariantImplementation : public RegressionVariant
//
//            class SumSquared                 : public VariantImplementation
//            class DifferenceSquared          : public VariantImplementation
//            class CrossProduct               : public VariantImplementation
//            class WeightedVariance           : public VariantImplementation
//
//            class WeightedCorrelatedVariance : public WeightedVariance
//            class WeightedCorrelatedVarianceAndTraits : public WeightedCorrelatedVariance
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress.h"

namespace SAGE   {
namespace SIBPAL {

class VariantImplementation : public RegressionVariant
{
  public:

    typedef RegressionVariant        base_type;
    typedef dependent_variable       trait_parameter;

    VariantImplementation(string           name,
                          TraitRegression& reg,
                          cerrorstream&    err = sage_cerr);

    virtual void regress();
    virtual void simulate();
    virtual void build();
    virtual void update();

    virtual std::string description() const;

    virtual void   weight_matrix(const sib_cluster&     ship,
                                  matrix&                W,
                                  const trait_parameter& trait,
                                  bool                   use_empirical_correlations,
                                  weight_status_type&    status,
                                  double                 fsib_var,
                                  double                 hsib_var) const;

    TraitRegression*  regression;
};

class SumSquared : public VariantImplementation
{
  public:

    typedef VariantImplementation  base_type;

    SumSquared(TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& trait,
                                 bool                   use_empirical_correlations, 
                                 bool                   center) const;
};

class DifferenceSquared : public VariantImplementation
{
  public:

    typedef VariantImplementation  base_type;

    DifferenceSquared(TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& trait,
                                 bool                   use_empirical_correlations, 
                                 bool                   center) const;
};

class CrossProduct : public VariantImplementation
{
  public:

    typedef VariantImplementation  base_type;

    CrossProduct(TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& trait,
                                 bool                   use_empirical_correlations, 
                                 bool                   center) const;
};

class WeightedVariance : public VariantImplementation
{
  public:

    typedef VariantImplementation  base_type;

    WeightedVariance(TraitRegression& reg, cerrorstream& err = sage_cerr);

    WeightedVariance(std::string name, TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual void    build();
    virtual void    update();

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& trait,
                                 bool                   use_empirical_correlations, 
                                 bool                   center) const;

    mutable TraitRegression sum_regression;
    mutable TraitRegression diff_regression;

    double  ss_d;
    double  ss_s;
    double  tss_d;
    double  tss_s;
};

class WeightedCorrelatedVariance : public WeightedVariance
{
  public:

    typedef WeightedVariance  base_type;

    WeightedCorrelatedVariance(TraitRegression& reg, cerrorstream& err = sage_cerr);

    WeightedCorrelatedVariance(std::string name, TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual bool    apply_weighting() const;

    virtual void    trait_vector(const sib_cluster&     ship,
                                 matrix&                y,
                                 const trait_parameter& trait,
                                 bool                   use_empirical_correlations, 
                                 bool                   center) const;

    virtual void    weight_matrix(const sib_cluster&     ship,
                                  matrix&                W,
                                  const trait_parameter& trait,
                                  bool                   use_empirical_correlations,
                                  weight_status_type&    status,
                                  double                 fsib_var,
                                  double                 hsib_var) const;

    void load_weights(const sib_cluster& ship, bool use_empirical_correlations) const;

  protected:

    // Various temporaries that I'd rather keep around than re-allocate every time
    mutable matrix Wy_d;
    mutable matrix Wy_s;
    mutable matrix W_s;
    mutable matrix W_d;
    mutable matrix y_s;
    mutable matrix y_d;
    mutable matrix W_sum;
    mutable matrix y_sum;
    mutable matrix temp;
};

class WeightedCorrelatedVarianceAndTraits : public WeightedCorrelatedVariance
{
  public:

    typedef WeightedCorrelatedVariance base_type;

    WeightedCorrelatedVarianceAndTraits(TraitRegression& reg, cerrorstream& err = sage_cerr);

    virtual string  description() const;

    virtual void    update();

    virtual void    weight_matrix(const sib_cluster&     ship,
                                  matrix&                W,
                                  const trait_parameter& trait,
                                  bool                   use_empirical_correlations,
                                  weight_status_type&    status,
                                  double                 fsib_var,
                                  double                 hsib_var) const;

  protected:

    // Various temporaries that I'd rather keep around than re-allocate every time
    mutable matrix W_sd;
    mutable matrix sigma;
    mutable matrix sigma_i;
    mutable matrix sigma_adj;
};

#include "regress_variants.ipp"

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
