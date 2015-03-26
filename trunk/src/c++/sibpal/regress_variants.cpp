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

#include "sibpal/regress_variants.h"

#define DEBUG_W3 0
#define DEBUG_W4 0

namespace SAGE   {
namespace SIBPAL {

//
//------------------------------------------------------------------------
//

VariantImplementation::VariantImplementation(string name, TraitRegression& reg, cerrorstream& err)
                     : base_type(name, err)
{
  regression = &reg;
}

void
VariantImplementation::regress()
{
  regression->do_regress(); 
}

void
VariantImplementation::simulate()
{
  regression->do_simulate();
}

void
VariantImplementation::build()
{
  update();
}

void
VariantImplementation::update()
{
  // Compute the sib dependent variable correlations
  regression->estimate_dependent_variable_correlation();
}

void
VariantImplementation::weight_matrix(const sib_cluster&     ship,
                                     matrix&                W,
                                     const trait_parameter& trait,
                                     bool                   use_empirical_correlations,
                                     weight_status_type&    status,
                                     double                 fsib_var,
                                     double                 hsib_var) const
{
#if 0
  cout << "VariantImplementation::weight_matrix()..." << endl;
#endif

  W.clear();
  size_t n = ship.valid_pair_count();

  if(!n)
  {
    W.setstate(matrix::badbit);
    return;
  }

  W = eye<double>(n);

  if( regression->get_model().get_analysis_options().identity_weight )
  {
    if( status == BESTW )
      status = INVERSEW;
    return;
  }

  double c = hsib_var / fsib_var;
  
  for( size_t i = 0; i < n; ++i )
  {
    if( ship[i].is_hsib_pair() )
    {
      W(i, i) = c;
    }
  }

  if( !use_empirical_correlations )
  {
    status = NORMALW;
    return;
  }

  double p2_ff = 1.0, p2_hh = 1.0, p2_fh = 1.0;
  double p1_ff = trait.p_fsib_empirical_correlation[1];
  double p1_hh = trait.p_hsib_empirical_correlation[1];
  double p1_fh = trait.p_fh_empirical_correlation[1];
  double p0_ff = trait.p_fsib_empirical_correlation[0];
  double p0_hh = trait.p_hsib_empirical_correlation[0];
  double p0_fh = trait.p_fh_empirical_correlation[0];

  sib_matrix_pattern p_ff(p2_ff, p1_ff, p0_ff, 0.);
  sib_matrix_pattern p_hh(p2_hh, p1_hh, p0_hh, 0.);
  sib_matrix_pattern p_fh(p2_fh, p1_fh, p0_fh, 0.);

  weights.weight_matrix_combined(ship, W, c, p_ff, p_hh, p_fh, status);

#if 0
  cout << "W =" << endl;
  print_matrix(W, cout);
  cout << "end of VariantImplementation::weight_matrix()..." << endl;
#endif

  return;
}

//
//------------------------------------------------------------------------
//

SumSquared::SumSquared(TraitRegression& reg, cerrorstream& err)
          : base_type("SUM", reg, err)
{}

void
SumSquared::trait_vector(const sib_cluster&     ship,
                         matrix&                y,
                         const trait_parameter& trait,
                         bool                   use_empirical_correlations,
                         bool                   center) const
{
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  y.resize_nofill(pair_count,1);

  pair<double, double> means;
  double tmean = 0;

  if(center)
    tmean = trait.info.mean();
#if 0
  cout << "tmean = " << tmean << endl;
#endif
  for( size_t j = 0; j < pair_count; ++j )
  {
    double t1 = regression->member_trait(ship[j].rels().pair.first,  trait.trait_index);
    double t2 = regression->member_trait(ship[j].rels().pair.second, trait.trait_index);

    if( regression->get_model().get_analysis_options().sibship_mean )
    {
      means = regression->get_sibship_means(ship[j]);
    }
    else if( regression->get_model().get_analysis_options().blup_mean )
    {
      means = regression->get_blup_means(ship[j]);
    }
    else
      means = regression->get_trait_means(ship[j]);
#if 0
  cout << "means = (" << means.first << ", " << means.second << ")"  << endl;
#endif
    double t1m  = t1 - means.first;
    double t2m  = t2 - means.second;

    y(j,0) = 0.50*(t1m + t2m)*(t1m + t2m)-tmean;
  }
}

//
//------------------------------------------------------------------------
//

DifferenceSquared::DifferenceSquared(TraitRegression& reg, cerrorstream& err)
                 : base_type("DIFFERENCE", reg, err)
{}

void
DifferenceSquared::trait_vector(const sib_cluster&     ship,
                                matrix&                y,
                                const trait_parameter& trait,
                                bool                   use_empirical_correlations,
                                bool                   center) const
{
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  double tmean = 0;
  if(center)
    tmean = trait.info.mean();

  y.resize_nofill(pair_count,1);

  for( size_t j = 0; j < pair_count; ++j )
  {
    double t1 = regression->member_trait(ship[j].rels().pair.first,  trait.trait_index);
    double t2 = regression->member_trait(ship[j].rels().pair.second, trait.trait_index);

    y(j,0) = -0.50 * (t1 - t2)*(t1 - t2) - tmean;
  }
}

//
//------------------------------------------------------------------------
//

CrossProduct::CrossProduct(TraitRegression& reg, cerrorstream& err)
            : base_type("PRODUCT",reg,err)
{}

void
CrossProduct::trait_vector(const sib_cluster&     ship,
                           matrix&                y,
                           const trait_parameter& trait,
                           bool                   use_empirical_correlations,
                           bool                   center) const
{
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  y.resize_nofill(pair_count,1);

  pair<double, double> means;
  double tmean = 0;

  if(center)
    tmean = trait.info.mean();

  for( size_t j = 0; j < pair_count; ++j )
  {
    double t1 = regression->member_trait(ship[j].rels().pair.first,  trait.trait_index);
    double t2 = regression->member_trait(ship[j].rels().pair.second, trait.trait_index);

    if( regression->get_model().get_analysis_options().sibship_mean )
    {
      means = regression->get_sibship_means(ship[j]);
    }
    else if( regression->get_model().get_analysis_options().blup_mean )
    {
      means = regression->get_blup_means(ship[j]);
    }
    else
      means = regression->get_trait_means(ship[j]);

    double t1m  = t1 - means.first;
    double t2m  = t2 - means.second;

    y(j,0) = t1m * t2m - tmean;
  }
}

//
//------------------------------------------------------------------------
//

WeightedVariance::WeightedVariance(TraitRegression& reg, cerrorstream& err)
                : base_type("WEIGHTED2", reg, err),
                  sum_regression(reg.get_pairs(), reg.error_stream()),
                  diff_regression(reg.get_pairs(), reg.error_stream())
{}

WeightedVariance::WeightedVariance(std::string name, TraitRegression& reg, cerrorstream& err)
                : base_type(name, reg, err),
                  sum_regression(reg.get_pairs(), reg.error_stream()),
                  diff_regression(reg.get_pairs(), reg.error_stream())
{}

void
WeightedVariance::build()
{
#if 0
  cout << "WeightedVariance::build()..." << endl;
  cout << "SUM reg ===" << endl;
#endif
  sum_regression.invalidate();
  sum_regression.set_model( regression->get_model() );
  sum_regression.get_model().invalidate();
  sum_regression.get_model().set_regression_method_name("SUM");
  sum_regression.copy_build(regression);
#if 0
  cout << "DIFF reg ===" << endl;
#endif
  diff_regression.invalidate();
  diff_regression.set_model( regression->get_model() );
  diff_regression.get_model().invalidate();
  diff_regression.get_model().set_regression_method_name("DIFF");
  diff_regression.copy_build(regression);

  update();

#if 0
  cout << "end of WeightedVariance::build()..." << endl;
#endif
}

void
WeightedVariance::update()
{
#if 0
  cout << "WeightedVariance::update()..." << endl;
#endif

  double fsv = regression->get_reg_results().get_full_sum_residual_variance();
  double hsv = regression->get_reg_results().get_half_sum_residual_variance();
  double fdv = regression->get_reg_results().get_full_diff_residual_variance();
  double hdv = regression->get_reg_results().get_half_diff_residual_variance();

#if 0
  cout << "SUM reg ===" << endl;
  cout << "fsv = " << fsv << ", hsv = " << hsv << endl;
#endif
  sum_regression.get_reg_results().set_full_sum_residual_variance(fsv);
  sum_regression.get_reg_results().set_half_sum_residual_variance(hsv);
  sum_regression.do_regress_univariate();

#if 0
  cout << "DIFF reg ===" << endl;
  cout << "fdv = " << fdv << ", hdv = " << hdv << endl;
#endif
  diff_regression.get_reg_results().set_full_diff_residual_variance(fdv);
  diff_regression.get_reg_results().set_half_diff_residual_variance(hdv);
  diff_regression.do_regress_univariate();

  ss_s = ss_d = 1;

  if(!sum_regression.valid() || !diff_regression.valid())
    return;

  ss_s = sum_regression.get_reg_results().get_residual_variance();
  ss_d = diff_regression.get_reg_results().get_residual_variance();
  tss_s = sum_regression.get_reg_results().get_total_variance();
  tss_d = diff_regression.get_reg_results().get_total_variance();

#if 0 //DEBUG_WEIGHTS
  cout << "W2:  ss_d=" <<  ss_d << ",  ss_s=" <<  ss_s << endl;
  cout << "W2: tss_d=" << tss_d << ", tss_s=" << tss_s << endl;
#endif

  // Compute the sib dependent variable correlations
  regression->estimate_dependent_variable_correlation();

  regression->get_reg_results().set_sum_residual_variance(ss_s);
  regression->get_reg_results().set_diff_residual_variance(ss_d);

  if( regression->get_use_pairs().first && !regression->get_use_pairs().second )
  {
    regression->get_reg_results().set_full_sum_residual_variance(ss_s);
    regression->get_reg_results().set_full_diff_residual_variance(ss_d);
  }
  else if( !regression->get_use_pairs().first && regression->get_use_pairs().second )
  {
    regression->get_reg_results().set_half_sum_residual_variance(ss_s);
    regression->get_reg_results().set_half_diff_residual_variance(ss_d);
  }

#if 0
  cout << "end of WeightedVariance::update()..." << endl;
#endif
}

void
WeightedVariance::trait_vector(const sib_cluster&     ship,
                               matrix&                y,
                               const trait_parameter& trait,
                               bool                   use_empirical_correlations,
                               bool                   center) const
{
#if 0
  cout << "WeightedVariance::trait_vector()..." << endl;
#endif
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  y.resize_nofill(pair_count,1);

  double w1 = ss_d;
  double w2 = ss_s;

  w[ship.valid_sib_count()] = make_pair(w1, w2);

  pair<double, double> means;
  double tmean = 0;

  if(center)
    tmean = trait.info.mean();

  for( size_t j = 0; j < pair_count; ++j )
  {
    double t1 = regression->member_trait(ship[j].rels().pair.first,  trait.trait_index);
    double t2 = regression->member_trait(ship[j].rels().pair.second, trait.trait_index);

    if( regression->get_model().get_analysis_options().sibship_mean )
    {
      means = regression->get_sibship_means(ship[j]);
    }
    else if( regression->get_model().get_analysis_options().blup_mean )
    {
      means = regression->get_blup_means(ship[j]);
    }
    else
      means = regression->get_trait_means(ship[j]);

    double t1m  = t1 - means.first;
    double t2m  = t2 - means.second;

    y(j,0) = 0.50 * (w2 * (t1m + t2m)*(t1m + t2m)
                   - w1 * (t1m - t2m)*(t1m - t2m))
           / (w1 + w2) - tmean;
  }

#if 0
  cout << "tmean = " << tmean << endl;
  cout << "y =" << endl;
  print_matrix(y, cout);
  cout << "end of WeightedVariance::trait_vector()..." << endl;
#endif
}

//
//------------------------------------------------------------------------
//

WeightedCorrelatedVariance::WeightedCorrelatedVariance(TraitRegression& reg, cerrorstream& err)
                          : base_type("WEIGHTED3", reg, err)
{}

WeightedCorrelatedVariance::WeightedCorrelatedVariance(string name, TraitRegression& reg, cerrorstream& err)
                          : base_type(name, reg, err)
{}

//"Weighted square trait sum and difference (W3)"
void
WeightedCorrelatedVariance::trait_vector(const sib_cluster&     ship,
                                         matrix&                y,
                                         const trait_parameter& trait,
                                         bool                   use_empirical_correlations,
                                         bool                   center) const
{

#if 0
  cout << "WeightedCorrelatedVariance::trait_vector()..." << endl;
#endif
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  // Load y_d and y_s
  sum_regression.trait_vector(ship, y_s, trait, use_empirical_correlations, true);
  diff_regression.trait_vector(ship, y_d, trait, use_empirical_correlations, true);

  // Load W_d and W_s
  load_weights(ship, use_empirical_correlations);

  w[ship.valid_sib_count()] = make_pair(W_d(0,0), W_s(0,0));

  if( ss_s == 0.0 )
    y = y_d;

  else if( ss_d == 0.0 )
    y = y_s;

  else
  {
    multiply(W_d, y_d, Wy_d);
    multiply(W_s, y_s, Wy_s);

    y_sum  = Wy_s;
    y_sum += Wy_d;

    W_sum  = W_s;
    W_sum += W_d;

#if 0
    cout << "y_s =" << endl;
    print_matrix(y_s, cout);
    cout << "y_d =" << endl;
    print_matrix(y_d, cout);
    cout << "W_s =" << endl;
    print_matrix(W_s, cout);
    cout << "W_d =" << endl;
    print_matrix(W_d, cout);
    cout << "Wy_s =" << endl;
    print_matrix(Wy_s, cout);
    cout << "Wy_d =" << endl;
    print_matrix(Wy_d, cout);
    cout << "y_sum =" << endl;
    print_matrix(y_sum, cout);
    cout << "W_sum =" << endl;
    print_matrix(W_sum, cout);
#endif

    SVDinverse(W_sum,temp);

    multiply(temp, y_sum, y);

#if 0
    cout << "temp =" << endl;
    print_matrix(temp, cout);
    cout << "y =" << endl;
    print_matrix(y, cout);
#endif

    if(center)
      y -= trait.info.mean();
  }

#if 0
  cout << "mean = " << trait.info.mean() << endl;
  cout << "y =" << endl;
  print_matrix(y, cout);
#endif

#if 0
  cout << "end of WeightedCorrelatedVariance::trait_vector()..." << endl;
#endif
}

void
WeightedCorrelatedVariance::weight_matrix(const sib_cluster&     ship,
                                          matrix&                W,
                                          const trait_parameter& trait,
                                          bool                   use_empirical_correlations,
                                          weight_status_type&    status,
                                          double                 fsib_var,
                                          double                 hsib_var) const
{
#if 0
  cout << "WeightedCorrelatedVariance::weight_matrix()..." << endl;
#endif

  int n = ship.valid_pair_count();

  if(!n)
  {
    W.setstate(matrix::badbit);
    return;
  }

  W.clear();

  // Load W_d and W_s
  load_weights(ship, use_empirical_correlations);

#if 0
  cout << "W_s = " << endl;
  print_matrix(W_s, cout);
  cout << "W_d = " << endl;
  print_matrix(W_d, cout);
#endif

  W  = W_s;
  W += W_d;

  if( status == NORMALW )
  {
    temp = W;
    SVDinverse(temp,W);
  }
  else
    status = INVERSEW;

#if 0
  cout << "W =" << endl;
  print_matrix(W, cout);
  cout << "end of WeightedCorrelatedVariance::weight_matrix()..." << endl;
#endif
}

void
WeightedCorrelatedVariance::load_weights(const sib_cluster& ship, bool use_empirical_correlations) const
{
#if 0
  cout << "WeightedCorrelatedVariance::load_weight()..." << endl;
#endif

  // trait_vector in W3
  trait_parameter& sum_trait  = sum_regression.get_model().get_trait();
  trait_parameter& diff_trait = diff_regression.get_model().get_trait();

  weight_status_type xstatus = INVERSEW;

  double fsv = 1.0;
  double hsv = 1.0;
  double fdv = 1.0;
  double hdv = 1.0;

  if( regression->get_use_pairs().first && regression->get_use_pairs().second )
  {
    fsv = regression->get_reg_results().get_full_sum_residual_variance();
    hsv = regression->get_reg_results().get_half_sum_residual_variance();
    fdv = regression->get_reg_results().get_full_diff_residual_variance();
    hdv = regression->get_reg_results().get_half_diff_residual_variance();
  }
#if 0
  cout << "fsv = " << fsv << ", hsv = " << hsv << endl;
#endif

  sum_regression.weight_matrix(ship, W_s, sum_trait, use_empirical_correlations, xstatus, fsv, hsv);

#if 0
  cout << "fdv = " << fdv << ", hdv = " << hdv << endl;
#endif

  diff_regression.weight_matrix(ship, W_d, diff_trait, use_empirical_correlations, xstatus, fdv, hdv);

#if 0
  cout << "W_s = " << endl;
  print_matrix(W_s, cout);
  cout << "W_d = " << endl;
  print_matrix(W_d, cout);
  cout << "ss_s  = " << ss_s << endl;
  cout << "ss_d  = " << ss_d << endl;
#endif

  if( ss_s != 0.0 )
    W_s /= ss_s;
  else
    W_s = zero<double>(ship.valid_pair_count());

  if( ss_d != 0.0 )
    W_d /= ss_d;
  else
    W_d = zero<double>(ship.valid_pair_count());

#if 0
  cout << "W_s = " << endl;
  print_matrix(W_s, cout);
  cout << "W_d = " << endl;
  print_matrix(W_d, cout);
#endif
}

//
//------------------------------------------------------------------------
//

WeightedCorrelatedVarianceAndTraits::WeightedCorrelatedVarianceAndTraits(TraitRegression& reg, cerrorstream& err)
                                   : base_type("WEIGHTED4", reg, err)
{}

void
WeightedCorrelatedVarianceAndTraits::update()
{
#if 0
  cout << "WeightedCorrelatedVarianceAndTraits::update()..." << endl;
#endif

  WeightedVariance::update();

  trait_parameter& trait = regression->get_model().get_trait();

  const trait_parameter& sum_trait  = sum_regression.get_model().get_trait();
  const trait_parameter& diff_trait = diff_regression.get_model().get_trait();

  SimpleCorrelationInfo          pool_correlation;
  SimpleCorrelationInfo          pool_fsib_correlation;
  SimpleCorrelationInfo          pool_hsib_correlation;
  SimpleCorrelationInfo          pool_fh_correlation;

  vector<SimpleCorrelationInfo>  p_sum_diff_correlations(3);
  vector<SimpleCorrelationInfo>  p_fsib_sum_diff_correlations(3);
  vector<SimpleCorrelationInfo>  p_hsib_sum_diff_correlations(3);
  vector<SimpleCorrelationInfo>  p_fh_sum_diff_correlations(3);

  matrix beta_s, beta_d;
  matrix r_s, r_d;
  matrix A_s, A_d;

  // Iterate over all sib_cluster
  for( size_t c = 0; c < sum_regression.get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = sum_regression.get_sib_clusters()[c];

    if( regression->get_model().is_x_linked() )
    {
      sum_regression.design_matrix_x(sc, A_s, false);
      diff_regression.design_matrix_x(sc, A_d, false);
    }
    else
    {
      sum_regression.design_matrix(sc, A_s, false);
      diff_regression.design_matrix(sc, A_d, false);
    }

    sum_regression.trait_vector(sc,  y_s, sum_trait, true, true);
    diff_regression.trait_vector(sc, y_d, diff_trait, true, true);

    sum_regression.my_gls.build_residuals(y_s,  A_s, r_s);
    diff_regression.my_gls.build_residuals(y_d, A_d, r_d);

#if 0
    cout << "r_s = " << endl;
    print_matrix(r_s, cout);
    cout << "r_d = " << endl;
    print_matrix(r_d, cout);
#endif

    for( size_t i = 0; i < sc.valid_pair_count(); ++i )
    {
      for( size_t j = i; j < sc.valid_pair_count(); ++j )
      {
        switch( sibs_shared( sc[i].rels().pair, sc[j].rels().pair ) )
        {
          case 2:

            p_sum_diff_correlations[2].add(r_s(i,0), r_d(j,0));
            pool_correlation.add(r_s(i,0), r_d(j,0));

            if( sc[i].is_fsib_pair() && sc[j].is_fsib_pair() )
            {
              p_fsib_sum_diff_correlations[2].add(r_s(i,0), r_d(j,0));
              pool_fsib_correlation.add(r_s(i,0), r_d(j,0));
            }
            else if( sc[i].is_hsib_pair() && sc[j].is_hsib_pair() )
            {
              p_hsib_sum_diff_correlations[2].add(r_s(i,0), r_d(j,0));
              pool_hsib_correlation.add(r_s(i,0), r_d(j,0));
            }
            else
            {
              p_fh_sum_diff_correlations[2].add(r_s(i,0), r_d(j,0));
              pool_fh_correlation.add(r_s(i,0), r_d(j,0));
            }

            break;

          case 1:  // One sib in common

            p_sum_diff_correlations[1].add(r_s(i,0), r_d(j,0));
            p_sum_diff_correlations[1].add(r_s(j,0), r_d(i,0));

            pool_correlation.add(r_s(i,0), r_d(j,0));
            pool_correlation.add(r_s(j,0), r_d(i,0));

            if( sc[i].is_fsib_pair() && sc[j].is_fsib_pair() )
            {
              p_fsib_sum_diff_correlations[1].add(r_s(i,0), r_d(j,0));
              p_fsib_sum_diff_correlations[1].add(r_s(j,0), r_d(i,0));

              pool_fsib_correlation.add(r_s(i,0), r_d(j,0));
              pool_fsib_correlation.add(r_s(j,0), r_d(i,0));
            }
            else if( sc[i].is_hsib_pair() && sc[j].is_hsib_pair() )
            {
              p_hsib_sum_diff_correlations[1].add(r_s(i,0), r_d(j,0));
              p_hsib_sum_diff_correlations[1].add(r_s(j,0), r_d(i,0));

              pool_hsib_correlation.add(r_s(i,0), r_d(j,0));
              pool_hsib_correlation.add(r_s(j,0), r_d(i,0));
            }
            else
            {
              p_fh_sum_diff_correlations[1].add(r_s(i,0), r_d(j,0));
              p_fh_sum_diff_correlations[1].add(r_s(j,0), r_d(i,0));

              pool_fh_correlation.add(r_s(i,0), r_d(j,0));
              pool_fh_correlation.add(r_s(j,0), r_d(i,0));
            }

            break;

          case 0:  // No sibs in common

            p_sum_diff_correlations[0].add(r_s(i,0), r_d(j,0));
            p_sum_diff_correlations[0].add(r_s(j,0), r_d(i,0));

            pool_correlation.add(r_s(i,0), r_d(j,0));
            pool_correlation.add(r_s(j,0), r_d(i,0));

            if( sc[i].is_fsib_pair() && sc[j].is_fsib_pair() )
            {
              p_fsib_sum_diff_correlations[0].add(r_s(i,0), r_d(j,0));
              p_fsib_sum_diff_correlations[0].add(r_s(j,0), r_d(i,0));

              pool_fsib_correlation.add(r_s(i,0), r_d(j,0));
              pool_fsib_correlation.add(r_s(j,0), r_d(i,0));
            }
            else if( sc[i].is_hsib_pair() && sc[j].is_hsib_pair() )
            {
              p_hsib_sum_diff_correlations[0].add(r_s(i,0), r_d(j,0));
              p_hsib_sum_diff_correlations[0].add(r_s(j,0), r_d(i,0));

              pool_hsib_correlation.add(r_s(i,0), r_d(j,0));
              pool_hsib_correlation.add(r_s(j,0), r_d(i,0));
            }
            else
            {
              p_fh_sum_diff_correlations[0].add(r_s(i,0), r_d(j,0));
              p_fh_sum_diff_correlations[0].add(r_s(j,0), r_d(i,0));

              pool_fh_correlation.add(r_s(i,0), r_d(j,0));
              pool_fh_correlation.add(r_s(j,0), r_d(i,0));
            }

            break;
          default:
            break;
        }
      }
    }
  }

  double p0 = p_sum_diff_correlations[0].correlation();
  double p1 = p_sum_diff_correlations[1].correlation();
  double p2 = p_sum_diff_correlations[2].correlation();

  double p0_ff = p_fsib_sum_diff_correlations[0].correlation();
  double p1_ff = p_fsib_sum_diff_correlations[1].correlation();
  double p2_ff = p_fsib_sum_diff_correlations[2].correlation();

  double p0_hh = p_hsib_sum_diff_correlations[0].correlation();
  double p1_hh = p_hsib_sum_diff_correlations[1].correlation();
  double p2_hh = p_hsib_sum_diff_correlations[2].correlation();

  double p0_fh = p_fh_sum_diff_correlations[0].correlation();
  double p1_fh = p_fh_sum_diff_correlations[1].correlation();
  double p2_fh = p_fh_sum_diff_correlations[2].correlation();

#if 0
  cout << "p0_icor = " << p0 << endl;
  cout << "p1_icor = " << p1 << endl;
  cout << "p2_icor = " << p2 << endl;
  cout << "p0_ff_icor = " << p0_ff << endl;
  cout << "p1_ff_icor = " << p1_ff << endl;
  cout << "p2_ff_icor = " << p2_ff << endl;
  cout << "p0_hh_icor = " << p0_hh << endl;
  cout << "p1_hh_icor = " << p1_hh << endl;
  cout << "p2_hh_icor = " << p2_hh << endl;
  cout << "p0_fh_icor = " << p0_fh << endl;
  cout << "p1_fh_icor = " << p1_fh << endl;
  cout << "p2_fh_icor = " << p2_fh << endl;
#endif

  if(!finite(p0))
    p0 = 0;
  if(!finite(p1))
    p1 = 0;
  if(!finite(p2))
    p2 = 0;

  if(!finite(p0_ff))
    p0_ff = 0;
  if(!finite(p1_ff))
    p1_ff = 0;
  if(!finite(p2_ff))
    p2_ff = 0;

  if(!finite(p0_hh))
    p0_hh = 0;
  if(!finite(p1_hh))
    p1_hh = 0;
  if(!finite(p2_hh))
    p2_hh = 0;

  if(!finite(p0_fh))
    p0_fh = 0;
  if(!finite(p1_fh))
    p1_fh = 0;
  if(!finite(p2_fh))
    p2_fh = 0;

#if 0
  double restricted_p2 = max(0.,max(p1,p2));
  double restricted_p1 = max(0.,p1);
  double restricted_p0 = max(0.,min(p1,p0));

  double pooled_cor = pool_correlation.correlation();

  if(    (p0 > p1 || p1 > p2)
      && regression->get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = restricted_p2 = max(0.,pooled_cor);

  trait.p_sum_diff_correlation[2] = restricted_p2;
  trait.p_sum_diff_correlation[1] = restricted_p1;
  trait.p_sum_diff_correlation[0] = restricted_p0;

  restricted_p2 = max(0.,max(p1_ff,p2_ff));
  restricted_p1 = max(0.,p1_ff);
  restricted_p0 = max(0.,min(p1_ff,p0_ff));

  pooled_cor = pool_fsib_correlation.correlation();

  if(    (p0_ff > p1_ff || p1_ff > p2_ff)
      && regression->get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = restricted_p2 = max(0.,pooled_cor);

  trait.p_fsib_sum_diff_correlation[2] = restricted_p2;
  trait.p_fsib_sum_diff_correlation[1] = restricted_p1;
  trait.p_fsib_sum_diff_correlation[0] = restricted_p0;

  restricted_p2 = max(0.,max(p1_hh,p2_hh));
  restricted_p1 = max(0.,p1_hh);
  restricted_p0 = max(0.,min(p1_hh,p0_hh));

  pooled_cor = pool_hsib_correlation.correlation();

  if(    (p0_hh > p1_hh || p1_hh > p2_hh)
      && regression->get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = restricted_p2 = max(0.,pooled_cor);

  trait.p_hsib_sum_diff_correlation[2] = restricted_p2;
  trait.p_hsib_sum_diff_correlation[1] = restricted_p1;
  trait.p_hsib_sum_diff_correlation[0] = restricted_p0;

  restricted_p2 = max(0.,max(p1_fh,p2_fh));
  restricted_p1 = max(0.,p1_fh);
  restricted_p0 = max(0.,min(p1_fh,p0_fh));

  pooled_cor = pool_fh_correlation.correlation();

  if(    (p0_fh > p1_fh || p1_fh > p2_fh)
      && regression->get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = restricted_p2 = max(0.,pooled_cor);

  trait.p_fh_sum_diff_correlation[2] = restricted_p2;
  trait.p_fh_sum_diff_correlation[1] = restricted_p1;
  trait.p_fh_sum_diff_correlation[0] = restricted_p0;
#else
  trait.p_sum_diff_correlation[2] = max(0.,p2);
  trait.p_sum_diff_correlation[1] = max(0.,p1);
  trait.p_sum_diff_correlation[0] = max(0.,p0);

  trait.p_fsib_sum_diff_correlation[2] = max(0.,p2_ff);
  trait.p_fsib_sum_diff_correlation[1] = max(0.,p1_ff);
  trait.p_fsib_sum_diff_correlation[0] = max(0.,p0_ff);

  trait.p_hsib_sum_diff_correlation[2] = max(0.,p2_hh);
  trait.p_hsib_sum_diff_correlation[1] = max(0.,p1_hh);
  trait.p_hsib_sum_diff_correlation[0] = max(0.,p0_hh);

  trait.p_fh_sum_diff_correlation[2] = max(0.,p2_fh);
  trait.p_fh_sum_diff_correlation[1] = max(0.,p1_fh);
  trait.p_fh_sum_diff_correlation[0] = max(0.,p0_fh);
#endif

#if 0
  cout << "end of WeightedCorrelatedVarianceAndTraits::update()..." << endl;
#endif
}

void
WeightedCorrelatedVarianceAndTraits::weight_matrix(const sib_cluster&     ship,
                                                   matrix&                W,
                                                   const trait_parameter& trait,
                                                   bool                   use_empirical_correlations,
                                                   weight_status_type&    status,
                                                   double                 fsib_var,
                                                   double                 hsib_var) const
{
#if 0
  cout << "WeightedCorrelatedVarianceAndTraits::weight_matrix()..." << endl;
#endif

  W.clear();
  size_t n = ship.valid_pair_count();

  if(!n)
  {
    W.setstate(matrix::badbit);
    return;
  }

  double p0_ff = trait.p_fsib_sum_diff_correlation[0];
  double p1_ff = trait.p_fsib_sum_diff_correlation[1];
  double p2_ff = trait.p_fsib_sum_diff_correlation[2];
  double p0_hh = trait.p_hsib_sum_diff_correlation[0];
  double p1_hh = trait.p_hsib_sum_diff_correlation[1];
  double p2_hh = trait.p_hsib_sum_diff_correlation[2];
  double p0_fh = trait.p_fh_sum_diff_correlation[0];
  double p1_fh = trait.p_fh_sum_diff_correlation[1];
  double p2_fh = trait.p_fh_sum_diff_correlation[2];

  W_sd.clear();
  W_sd.resize_fill(n, n, 0.0);

  double c = hsib_var / fsib_var;
  
  for( size_t i = 0; i < n; ++i )
  {
    if( ship[i].is_fsib_pair() )
    {
      W_sd(i, i) = p2_ff;
    }
    else if( ship[i].is_hsib_pair() )
    {
      W_sd(i, i) = c * p2_hh;
    }
  }

  if( use_empirical_correlations )
  {

    weight_status_type xstatus = NORMALW;

    sib_matrix_pattern p_ff(p2_ff, p1_ff, p0_ff, 0.);
    sib_matrix_pattern p_hh(p2_hh, p1_hh, p0_hh, 0.);
    sib_matrix_pattern p_fh(p2_fh, p1_fh, p0_fh, 0.);

    weights.weight_matrix_combined(ship, W_sd, c, p_ff, p_hh, p_fh, xstatus);
  }

#if 0
  cout << "p0_ff = " << p0_ff << ", p1_ff = " << p1_ff << ", p2_ff = " << p2_ff << endl;
  cout << "p0_hh = " << p0_hh << ", p1_hh = " << p1_hh << ", p2_hh = " << p2_hh << endl;
  cout << "p0_fh = " << p0_fh << ", p1_fh = " << p1_fh << ", p2_fh = " << p2_fh << endl;

  cout << "W_sd = " << endl;
  print_matrix(W_sd, cout);
#endif

  W_sd *= sqrt(ss_s*ss_d);

#if 0
  cout << "W_sd = " << endl;
  print_matrix(W_sd, cout);
#endif

  // Load W_d and W_s
  load_weights(ship, use_empirical_correlations);

#if 0
  cout << "W_s = " << endl;
  print_matrix(W_s, cout);
  cout << "W_d = " << endl;
  print_matrix(W_d, cout);
#endif

  W_sum  = W_s;
  W_sum += W_d;

  SVDinverse(W_sum,sigma);

  sigma_adj  = W_s * W_sd * W_d;

#if 0
  cout << "W_s * W_sd * W_d = " << endl;
  print_matrix(sigma_adj, cout);
  cout << "sigma_adj_transpose = " << endl;
  print_matrix(transpose(sigma_adj), cout);
#endif

  sigma_adj += transpose(sigma_adj);

  sigma_adj += W_sum;

  W = sigma * sigma_adj * sigma;

  if(status == INVERSEW)
  {
    temp = W;
    SVDinverse(temp,W);
  }
  else
    status = NORMALW;
    //status = INVERSEW;

#if 0
  cout << "W =" << endl;
  print_matrix(W, cout);
  cout << "end of WeightedCorrelatedVarianceAndTraits::weight_matrix()..." << endl;
#endif

  return;
}

//
//------------------------------------------------------------------------
//

boost::shared_ptr<RegressionVariant>
getRegressionVariant(TraitRegression& reg, string name)
{
  typedef boost::shared_ptr<RegressionVariant>  variant_pointer;

  if( !name.size() )
    name = reg.get_model().get_regression_method_name();

  if( !name.size() )
    return variant_pointer();

  name = toUpper(name);

  if(name == "DIFF" || name == "DIFFERENCE" || name == "HE" || name == "HE1" )
    return variant_pointer(new DifferenceSquared(reg, reg.error_stream()));

  if(name == "SUM" || name == "SUMM")
    return variant_pointer(new SumSquared(reg, reg.error_stream()));

  if(name == "PROD" || name == "PRODUCT" || name == "HE2")
    return variant_pointer(new CrossProduct(reg, reg.error_stream()));

  if(name == "W2" || name == "WEIGHTED2" || name == "WEIGHT2")
    return variant_pointer(new WeightedVariance(reg, reg.error_stream()));

  if(name == "W3" || name == "WEIGHTED3" || name == "WEIGHT3")
    return variant_pointer(new WeightedCorrelatedVariance(reg, reg.error_stream()));

  if(name == "W4" || name == "WEIGHTED4" || name == "WEIGHT4")
    return variant_pointer(new WeightedCorrelatedVarianceAndTraits(reg, reg.error_stream()));

  return variant_pointer();
};

} // end of namespace SIBPAL
} // end of namespace SAGE
