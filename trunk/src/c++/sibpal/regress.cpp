//=============================================================================
// File:    regress.cpp
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains implementation for following data structures.
//
//            class RegressionVariant
//            class TraitRegression
//
// Copyright (c) 2001 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/regress.h"

#undef DEBUG_TRAIT_VECTOR
#undef DEBUG_MARKER_VECTOR
#define DEBUGCOR(x)

#if __KCC
#  define BROKEN_STL_RANDOM_SHUFFLE
#endif

namespace SAGE   {
namespace SIBPAL {

//
//------------------------------------------------------------------------
//

RegressionVariant::RegressionVariant(std::string name, cerrorstream& err)
                 : errors(err)
{
  my_name = name;
}

RegressionVariant::~RegressionVariant()
{}
                   
//
//------------------------------------------------------------------------
//

//------------------------------------------------------------
// public member functions
//------------------------------------------------------------

TraitRegression::TraitRegression(relative_pairs& p, cerrorstream& err)
               : pairs(p), errors(err)
{
  invalidate();
  reset();
}

TraitRegression::~TraitRegression()
{
  disown_regression_method();
}

void
TraitRegression::invalidate()
{

  my_valid_parameter_count = my_valid_trait_count = 0;
  my_built = false;
  my_simulating = false;

  my_reg_results.clear();
  my_reg_results_intercept.clear();

  if( my_use_fsib && !my_use_hsib )
  {
    my_reg_results.clear_full();
  }
  else if( !my_use_fsib && my_use_hsib )
  {
    my_reg_results.clear_half();
  }

  my_valid_regression = false;

  return;
}

void
TraitRegression::reset()
{
  my_built_trait_all_sibs_info = false;
  my_use_fsib = my_use_hsib = false;
  my_fsib_pair_count = my_hsib_pair_count = 0;

  my_use_x_pair.resize(3, true);
  my_fix_x_pair.resize(3, false);
  my_x_pair_count.resize(3, 0);
  my_x_types.resize(0);
  my_x_types.push_back(MIXED);

  my_fsib_intercept.resize(0);
  my_hsib_intercept.resize(0);
  my_x_intercept.resize(3);
  my_x_intercept_types.resize(0);

  return;
}

void
TraitRegression::regress()
{
  if( !regression_method() )
    set_regression_method( getRegressionVariant(*this) );

  if( regression_method() )
  {
    regression_method()->regress();
  }
  else
    errors << priority(error) << "Invalid or unsupported regression type: "
           << get_model().get_regression_method_name() << endl;

  return;
}

void
TraitRegression::do_regress()
{
  make_filter();
  estimate_trait_sample_mean();

  if( get_model().get_data_options().use_pairs == SIB )
  {
    if( pairs.fsib_pair_count() && pairs.hsib_pair_count() )
    {
#if 0
      cout << "using both full & half sibs" << endl;
      cout << "*** first using full sibs only" << endl;
#endif
      my_use_fsib = true;
      my_use_hsib = false;

      my_fsib_intercept.resize(1, 1.0);
      my_hsib_intercept.resize(0);
      
      regress_univariate(false, false);

      if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
      {
#if 0
        const w_map& ws = regression_method()->get_w();

        w_map::const_iterator wi = ws.begin();
        for( ; wi != ws.end(); ++wi )
          cout << "sibship size " << wi->first
               << " : w1 = " << wi->second.first
               << ", w2 = " << wi->second.second << endl;
#endif
        my_reg_results.set_full_w(regression_method()->get_w());
        regression_method()->clear_w();
      }
#if 0
      my_reg_results.dump(cout);
      cout << "*** second using half sibs only" << endl;
#endif
      my_use_fsib = false;
      my_use_hsib = true;

      my_fsib_intercept.resize(0);
      my_hsib_intercept.resize(1, 1.0);

      regress_univariate(false, false);

      if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
      {
#if 0
        const w_map& ws = regression_method()->get_w();

        w_map::const_iterator wi = ws.begin();
        for( ; wi != ws.end(); ++wi )
          cout << "sibship size " << wi->first
               << " : w1 = " << wi->second.first
               << ", w2 = " << wi->second.second << endl;
#endif
        my_reg_results.set_half_w(regression_method()->get_w());
        regression_method()->clear_w();
      }
#if 0
      my_reg_results.dump(cout);
      cout << "*** third using full & half sibs" << endl;
#endif
      if(    my_reg_results.get_full_residual_variance() < EPS
          || my_reg_results.get_half_residual_variance() < EPS )
      {
        errors << priority(error)
               << "Can't combine full & half sib pairs..." << endl;

        return;
      }

      // do regress_combined
      my_use_fsib = true;
      my_use_hsib = true;

      my_fsib_intercept.resize(2, 1.0);
      my_fsib_intercept[1] = 0.0;
      my_hsib_intercept.resize(2, 1.0);
      my_hsib_intercept[0] = 0.0;

      regress_univariate(false, true);
    }
    else if( pairs.fsib_pair_count() )
    {
      my_use_fsib = true;
      my_use_hsib = false;

      my_fsib_intercept.resize(1, 1.0);
      my_hsib_intercept.resize(0);

      //cout << "using full sibs" << endl;

      regress_univariate(true, true);

      if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
      {
        my_reg_results.set_full_w(regression_method()->get_w());
        regression_method()->clear_w();
      }
    }
    else if( pairs.hsib_pair_count() )
    {
      my_use_fsib = false;
      my_use_hsib = true;

      my_fsib_intercept.resize(0);
      my_hsib_intercept.resize(1, 1.0);

      //cout << "using half sibs" << endl;

      regress_univariate(false, true);

      if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
      {
        my_reg_results.set_half_w(regression_method()->get_w());
        regression_method()->clear_w();
      }
    }
    else
      errors << priority(error)
             << "No sib pairs exist.  Exiting..." << endl;
  }
  else if( get_model().get_data_options().use_pairs == FSIB )
  {
    if( get_model().is_x_linked() )
    {
      if( pairs.fsib_pair_count() )
      {
        my_use_fsib = true;
        my_use_hsib = false;

        my_fsib_intercept.resize(1, 1.0);
        my_hsib_intercept.resize(0);

        if( !get_model().get_data_options().use_mm_pair || !pairs.mm_pair_count() )
          my_use_x_pair[MM] = false;

        if( !get_model().get_data_options().use_mf_pair || !pairs.mf_pair_count() )
          my_use_x_pair[MF] = false;

        if( !get_model().get_data_options().use_ff_pair || !pairs.ff_pair_count() )
          my_use_x_pair[FF] = false;

        build_x_types();

        regress_univariate(true, true);

        if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
        {
          my_reg_results.set_full_w(regression_method()->get_w());
          regression_method()->clear_w();
        }
      }
      else
        errors << priority(error)
               << "No full sib pairs exist.  Exiting..." << endl;

    }
    else
    {  
      if( pairs.fsib_pair_count() )
      {
        //cout << "using full sibs" << endl;

        my_use_fsib = true;
        my_use_hsib = false;

        my_fsib_intercept.resize(1, 1.0);
        my_hsib_intercept.resize(0);

        regress_univariate(true, true);

        if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
        {
          my_reg_results.set_full_w(regression_method()->get_w());
          regression_method()->clear_w();
        }
      }
      else
        errors << priority(error)
               << "No full sib pairs exist.  Exiting..." << endl;
    }
  }
  else if( get_model().get_data_options().use_pairs == HSIB )
  {
    if( pairs.hsib_pair_count() )
    {
      //cout << "using half sibs" << endl;

      my_use_fsib = false;
      my_use_hsib = true;

      my_fsib_intercept.resize(0);
      my_hsib_intercept.resize(1, 1.0);

      regress_univariate(false, true);

      if( toUpper(get_model().get_regression_method_name().substr(0,1)) == "W" )
      {
        my_reg_results.set_half_w(regression_method()->get_w());
        regression_method()->clear_w();
      }
    }
    else
      errors << priority(error)
             << "No half sib pairs exist.  Exiting..." << endl;
  }

  return;
}

void
TraitRegression::estimate_dependent_variable_correlation()
{
  SampleInfo             info;
  SampleInfo             mm_pair_info;
  SampleInfo             mf_pair_info;
  SampleInfo             ff_pair_info;

  SimpleCorrelationInfo          pool_correlation;
  vector<SimpleCorrelationInfo>  p_correlation(2);

  matrix y;

  dependent_variable& current_trait = get_model().get_trait();

  // Iterate over all sibships
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    trait_vector(sc, y, current_trait, false, false);

    for(size_t i = 0; i < sc.valid_pair_count(); ++i)
    {
      info += y(i,0);

      if( is_mm_pair(sc[i].rels().pair) )
        mm_pair_info += y(i,0);
      else if( is_mf_pair(sc[i].rels().pair) )
        mf_pair_info += y(i,0);
      else if( is_ff_pair(sc[i].rels().pair) )
        ff_pair_info += y(i,0);

      for(size_t j=i+1; j < sc.valid_pair_count(); ++j)
      {
        switch( sibs_shared( sc[i].rels().pair, sc[j].rels().pair ) )
        {
          case 2:  break;  // Correlation is 1

          case 1:  // One sib in common

            p_correlation[1].add(y(i,0), y(j,0));
            p_correlation[1].add(y(j,0), y(i,0));

            pool_correlation.add(y(i,0), y(j,0));
            pool_correlation.add(y(j,0), y(i,0));

            break;

          case 0:  // No sibs in common

            p_correlation[0].add(y(i,0), y(j,0));
            p_correlation[0].add(y(j,0), y(i,0));

            pool_correlation.add(y(i,0), y(j,0));
            pool_correlation.add(y(j,0), y(i,0));

            break;

          default: // Unrelated pairs
            break;
        }
      }        
    }
  }

  current_trait.info = info;

  current_trait.mm_pair_info = mm_pair_info;
  current_trait.mf_pair_info = mf_pair_info;
  current_trait.ff_pair_info = ff_pair_info;

  double p0 = p_correlation[0].correlation();
  double p1 = p_correlation[1].correlation();

  double restricted_p0 = max(0.,p0);
  double restricted_p1 = max(0.,max(p0,p1));

  double pooled_cor = pool_correlation.correlation();

  if( p0 > p1 && get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = max(0.,pooled_cor);

  current_trait.p_correlation[0] = restricted_p0;
  current_trait.p_correlation[1] = restricted_p1;

  if( my_use_fsib && !my_use_hsib )
  {
    current_trait.p_fsib_correlation[0] = restricted_p0;
    current_trait.p_fsib_correlation[1] = restricted_p1;
  }
  else if( !my_use_fsib && my_use_hsib )
  {
    current_trait.p_hsib_correlation[0] = restricted_p0;
    current_trait.p_hsib_correlation[1] = restricted_p1;
  }
#if 0
  cout << "p0 = " << p0 << ", p1 = " << p1 << ", pooled cor = " << pooled_cor << endl;
  cout << "current_trait.p0_correlation = " << current_trait.p_correlation[0] << endl;
  cout << "current_trait.p1_correlation = " << current_trait.p_correlation[1] << endl;
#endif

  return;
}

void
TraitRegression::simulate()
{
#if 0
  cout << "TraitRegression::simulate()..." << endl;
#endif

  if( simulating() )
    return;

  my_simulating = true;

  if( !regression_method() )
    set_regression_method( getRegressionVariant(*this) );

  if( regression_method() )
    regression_method()->simulate();
  else
    errors << priority(error) << "Invalid or unsupported regression type: "
           << get_model().get_regression_method_name() << endl;

  my_simulating = false;
}

void
TraitRegression::do_simulate()
{
  simulate_univariate();

  return;
}

//------------------------------------------------------------
// protected member functions
//------------------------------------------------------------

void
TraitRegression::make_filter()
{
  my_filter = pair_filter();

  pair_filter& pf = my_filter;

  size_t p = get_model().get_parameter_count();
  size_t s = get_model().get_subset_count();

  double tolerance = INF;

  if( get_model().get_data_options().skip_uninformative_pairs )
    tolerance = 0.001;

  for( size_t i = 0; i < p; ++i )
  {
    const independent_variable& param = get_model().get_parameter(i);

    for(size_t j = 0; j < param.markers.size(); ++j)
    {
      size_t m = param.markers[j].marker_index;
      if( m < pairs.marker_count() )
        pf.add_marker(m, tolerance);
    }

    for(size_t j = 0; j < param.covariates.size(); ++j)
    {
      size_t c = param.covariates[j].covariate_index;

      if( param.covariates[j].operation == covariate_type::pair )
      {
        if( c < pairs.pair_covariate_count() )
          pf.add_pair_covariate(c);
      }
      else if(c < pairs.trait_count())
        pf.add_trait(c);
    }
  }

  size_t tt = get_model().get_trait().trait_index;
  pf.add_trait(tt);

  for(size_t j = 0; j < s; ++j)
  {
    size_t tt = get_model().get_subset(j).trait_index;
    if(tt < pairs.trait_count() )
      pf.add_trait(tt, (size_t)2, EPS);
  }

  return;
}

void
TraitRegression::estimate_trait_sample_mean()
{
  if( my_built_trait_all_sibs_info )
    return;
  
  sib_set my_all_sibs;
  my_all_sibs.clear();

  get_model().get_trait().trait_all_sibs_info.clear();
  get_model().get_trait().trait_all_sibs_info.sample_adjustment(1);

  // X-linkage
  //
  get_model().get_trait().trait_male_sibs_info.clear();
  get_model().get_trait().trait_male_sibs_info.sample_adjustment(1);

  get_model().get_trait().trait_female_sibs_info.clear();
  get_model().get_trait().trait_female_sibs_info.sample_adjustment(1);
  //

  // Iterate over all sib_clusters
  sibship_cluster_const_iterator c_iter = pairs.sibship_cluster_begin();
  sibship_cluster_const_iterator c_end  = pairs.sibship_cluster_end();

  for( ; c_iter != c_end; ++c_iter )
  {
    sib_cluster sc_all(c_iter->fsib_pairs, c_iter->hsib_pairs, pairs, true, true, my_filter);

    if( !sc_all.valid_pair_count() ) continue;

    sib_set::const_iterator si = sc_all.get_valid_sibs().begin();
    for( ; si != sc_all.get_valid_sibs().end(); ++si )
    {
      my_all_sibs.insert(*si);
      get_model().get_trait().trait_all_sibs_info += member_trait( *si, get_model().get_trait().trait_index );

      if( (*si)->is_male() )
        get_model().get_trait().trait_male_sibs_info += member_trait( *si, get_model().get_trait().trait_index );
      else
        get_model().get_trait().trait_female_sibs_info += member_trait( *si, get_model().get_trait().trait_index );
    }
  }

#if 0
  cout << "all_sibs.size = " << my_all_sibs.size() << endl;

//  sib_set::const_iterator pi = my_all_sibs.begin();
//  for( ; pi != my_all_sibs.end(); ++pi )
//  {
//    cout << (*pi)->name() << endl;
//  }
//  cout << endl;

  cout << "trait sample mean = " << get_model().get_trait().trait_all_sibs_info.mean() << endl;
  cout << "trait sample vari = " << get_model().get_trait().trait_all_sibs_info.variance() << endl;

  cout << "trait sample male mean = " << get_model().get_trait().trait_male_sibs_info.mean() << endl;
  cout << "trait sample male vari = " << get_model().get_trait().trait_male_sibs_info.variance() << endl;

  cout << "trait sample female mean = " << get_model().get_trait().trait_female_sibs_info.mean() << endl;
  cout << "trait sample female vari = " << get_model().get_trait().trait_female_sibs_info.variance() << endl;
#endif

  my_built_trait_all_sibs_info = true;

  return;
}

void
TraitRegression::regress_univariate(bool do_sim, bool check_intercept)
{
  make_sib_cluster();

  if( pair_count() == 0 )
  {
    errors << priority(error)
           <<  "No pairs found!  Skipping regression." << endl;

    invalidate();
    return;         // Error no pairs
  }

  // invalidate()
  // normalize_effects()
  // estimate_independent_variable_parameters()
  // estimate_trait_sum_diff_correlation()
  // estimate_trait_correlation()
  // estimate_sibship_blup_mean()
  // regression_method()->build()
  build();

  if( valid_parameter_count() == 0 )
  {
    errors << priority(error) << "No valid parameters.  Skipping regression."
           << endl;

    invalidate();
    return;
  }

  if( pair_count() - valid_parameter_count() < 4 )
  {
    errors << priority(warning)
           << "Too few pairs to perform regression! (n=" << pair_count()
           << ",vp=" << valid_parameter_count()
           << ",vt=" << valid_trait_count() << ")  Skipping regression." << endl;

    invalidate();
    return;
  }

  if( !built() )
  {
    errors << priority(error)
           << "Error building regression data.  Skipping regression."
           << endl;

    invalidate();
    return;
  }

  if( get_model().is_x_linked() )
  {
#if 0
    cout << "mm = " << mm_pair_count() << endl;
    cout << "mf = " << mf_pair_count() << endl;
    cout << "ff = " << ff_pair_count() << endl;
#endif

    if( my_use_x_pair[MM] && (mm_pair_count() - valid_parameter_count() < 4) )
    {
      errors << priority(warning)
             << "Too few brother-brother pairs to perform regression! (n=" << mm_pair_count()
             << ",vp=" << valid_parameter_count()
             << ",vt=" << valid_trait_count() << ")"
             << "  Brother-brother pairs not used in the analysis." << endl;

      my_use_x_pair[MM] = false;
    }

    if( my_use_x_pair[MF] && (mf_pair_count() - valid_parameter_count() < 4) )
    {
      errors << priority(warning)
             << "Too few brother-sister pairs to perform regression! (n=" << mf_pair_count()
             << ",vp=" << valid_parameter_count()
             << ",vt=" << valid_trait_count() << ")"
             << "  Brother-sister pairs not used in the analysis." << endl;

      my_use_x_pair[MF] = false;
    }

    if( my_use_x_pair[FF] && (ff_pair_count() - valid_parameter_count() < 4) )
    {
      errors << priority(warning)
             << "Too few sister-sister pairs to perform regression! (n=" << ff_pair_count()
             << ",vp=" << valid_parameter_count()
             << ",vt=" << valid_trait_count() << ")"
             << "  Sister-sister pairs not used in the analysis." << endl;

      my_use_x_pair[FF] = false;
    }

    if( !my_use_x_pair[MM] && !my_use_x_pair[MF] && !my_use_x_pair[FF] )
    {
      errors << priority(error)
             << "Too few pairs to perform regression on X-linked location!  Skipping regression." << endl;

      invalidate();
      return;
    }

    build_x_types();
  }

  do_regress_univariate();

#if 0
  cout << "done with do_regress_univariate().. "
       << "do_sim = " << do_sim << " " << simulating()
       << " " << get_model().get_pvalue_options().is_on << endl;
#endif

  if( check_intercept && valid() )
  {
    check_intercepts();
    print_optional_output();
  }

  if( get_model().is_x_linked() && valid() )
    do_regress_F_test();

  if( do_sim && !simulating() && valid() && get_model().get_pvalue_options().is_on )
    if( my_use_fsib && !my_use_hsib ) // simulation only allowed for full sib pairs only
      simulate();

  return;
}

void
TraitRegression::build_x_types()
{
  my_x_types.resize(0);
  my_x_intercept_types.resize(0);

  if( my_use_x_pair[MM] )
  {
    my_x_types.push_back(MM);
    my_x_intercept_types.push_back(MM);
  }
  if( my_use_x_pair[MF] )
  {
    my_x_types.push_back(MF);
    my_x_intercept_types.push_back(MF);
  }
  if( my_use_x_pair[FF] )
  {
    my_x_types.push_back(FF);
    my_x_intercept_types.push_back(FF);
  }

  if( my_x_types.size() == 3 )
  {
    my_x_intercept[MM].resize(3, 0.0);
    my_x_intercept[MF].resize(3, 0.0);
    my_x_intercept[FF].resize(3, 0.0);
    my_x_intercept[MM][0] = 1.0;
    my_x_intercept[MF][1] = 1.0;
    my_x_intercept[FF][2] = 1.0;
  }
  else if( my_x_types.size() == 2 )
  {
    my_x_intercept[MM].resize(2, 0.0);
    my_x_intercept[MF].resize(2, 0.0);
    my_x_intercept[FF].resize(2, 0.0);

    if( my_use_x_pair[MM] )
      my_x_intercept[MM][0] = 1.0;
    if( my_use_x_pair[MF] )
      if( !my_use_x_pair[MM] )
        my_x_intercept[MF][0] = 1.0;
      else
        my_x_intercept[MF][1] = 1.0;

    if( my_use_x_pair[FF] )
      my_x_intercept[FF][1] = 1.0;
  }
  else
  {
    my_x_intercept[MM].resize(1, 1.0);
    my_x_intercept[MF].resize(1, 1.0);
    my_x_intercept[FF].resize(1, 1.0);
  }

  return;
}

void
TraitRegression::do_regress_univariate()
{
#if 0
  cout << "do_regress_univariate()..." << endl;
  if( my_fsib_intercept.size() )
  {
    cout << "  my_fsib_intercept: ";
    for( size_t c = 0; c < my_fsib_intercept.size(); ++c )
      cout << my_fsib_intercept[c] << " ";
    cout << endl;
  }
  if( my_hsib_intercept.size() )
  {
    cout << "  my_hsib_intercept: ";
    for( size_t c = 0; c < my_hsib_intercept.size(); ++c )
      cout << my_hsib_intercept[c] << " ";
    cout << endl;
  }
  cout << "  use_x_pair: " << my_use_x_pair[MM] << ", " << my_use_x_pair[MF] << ", " << my_use_x_pair[FF] << endl;
  cout << "  fix_x_pair: " << my_fix_x_pair[MM] << ", " << my_fix_x_pair[MF] << ", " << my_fix_x_pair[FF] << endl;
  cout << "  use x types : " << my_x_types.size();
  cout << " -";
  for( size_t x = 0; x < my_x_types.size(); ++x )
  {
    if( my_x_types[x] == MM )
      cout << " MM";
    if( my_x_types[x] == MF )
      cout << " MF";
    if( my_x_types[x] == FF )
      cout << " FF";
    if( my_x_types[x] == MIXED )
      cout << " MIXED";
  }
  cout << endl;
  if( my_x_intercept_types.size() )
  {
    cout << "  intercept x types : " << my_x_intercept_types.size();
    cout << " -";
    for( size_t x = 0; x < my_x_intercept_types.size(); ++x )
    {
      if( my_x_intercept_types[x] == MM )
        cout << " MM";
      if( my_x_intercept_types[x] == MF )
        cout << " MF";
      if( my_x_intercept_types[x] == FF )
        cout << " FF";
    }
    cout << endl;
  }
  if( my_x_intercept[MM].size() )
  {
    cout << "  my_MM_intercept: ";
    for( size_t c = 0; c < my_x_intercept[MM].size(); ++c )
      cout << my_x_intercept[MM][c] << " ";
    cout << endl;    
  }
  if( my_x_intercept[MF].size() )
  {
    cout << "  my_MF_intercept: ";
    for( size_t c = 0; c < my_x_intercept[MF].size(); ++c )
      cout << my_x_intercept[MF][c] << " ";
    cout << endl;    
  }
  if( my_x_intercept[FF].size() )
  {
    cout << "  my_FF_intercept: ";
    for( size_t c = 0; c < my_x_intercept[FF].size(); ++c )
      cout << my_x_intercept[FF][c] << " ";
    cout << endl;    
  }
#endif

  size_t vi = valid_intercept_count();
  size_t x_type_count = my_x_types.size();

  size_t m = x_type_count * valid_parameter_count();

  m += vi;

  my_gls.reset(m);

  my_gls.set_leverage_adjustment( get_model().get_analysis_options().leverage_adjustment );

  double fv = my_reg_results.get_full_residual_variance();
  double hv = my_reg_results.get_half_residual_variance();

#if 0
  if( my_use_fsib && my_use_hsib )
  {
    cout << "fv = " << fv << ", hv = " << hv << endl;
    cout << "fsv = "   << my_reg_results.get_full_sum_residual_variance()
         << ", hsv = " << my_reg_results.get_half_sum_residual_variance() << endl;
    cout << "fdv = "   << my_reg_results.get_full_diff_residual_variance()
         << ", hdv = " << my_reg_results.get_half_diff_residual_variance() << endl;
  }
#endif

  if( toUpper(get_model().get_regression_method_name()) == "SUM" )
  {
    fv = my_reg_results.get_full_sum_residual_variance();
    hv = my_reg_results.get_half_sum_residual_variance();
  }
  else if( toUpper(get_model().get_regression_method_name()) == "DIFF" )
  {
    fv = my_reg_results.get_full_diff_residual_variance();
    hv = my_reg_results.get_half_diff_residual_variance();
  }

  map<const sib_cluster*, matrix> final_As;
  map<const sib_cluster*, matrix> final_Ws;
  map<const sib_cluster*, matrix> final_Ys;
  map<const sib_cluster*, matrix> final_Rs;

  do_regress_univariate_sub(final_As, final_Ws, final_Ys, final_Rs, m, fv, hv);

  if( get_model().get_analysis_options().robust_variance )
    compute_robust_variance(my_gls, true);

  my_reg_results.set_residual_variance(my_gls.residual_variance(0,0));
  my_reg_results.set_total_variance(my_gls.total_variance(0,0));

  if( toUpper(get_model().get_regression_method_name()) == "SUM" )
  {
    my_reg_results.set_sum_residual_variance(my_gls.residual_variance(0,0));
  }
  else if( toUpper(get_model().get_regression_method_name()) == "DIFF" )
  {
    my_reg_results.set_diff_residual_variance(my_gls.residual_variance(0,0));
  }

  if( my_use_fsib && !my_use_hsib )
  {
    my_reg_results.set_full_residual_variance(my_gls.residual_variance(0,0));
    my_reg_results.set_full_total_variance(my_gls.total_variance(0,0));

    if( toUpper(get_model().get_regression_method_name()) == "SUM" )
    {
      my_reg_results.set_full_sum_residual_variance(my_gls.residual_variance(0,0));
    }
    else if( toUpper(get_model().get_regression_method_name()) == "DIFF" )
    {
      my_reg_results.set_full_diff_residual_variance(my_gls.residual_variance(0,0));
    }
  }
  else if( !my_use_fsib && my_use_hsib )
  {
    my_reg_results.set_half_residual_variance(my_gls.residual_variance(0,0));
    my_reg_results.set_half_total_variance(my_gls.total_variance(0,0));

    if( toUpper(get_model().get_regression_method_name()) == "SUM" )
    {
      my_reg_results.set_half_sum_residual_variance(my_gls.residual_variance(0,0));
    }
    else if( toUpper(get_model().get_regression_method_name()) == "DIFF" )
    {
      my_reg_results.set_half_diff_residual_variance(my_gls.residual_variance(0,0));
    }
  }

#if 0
  if( !get_model().is_x_linked() )
  {
    cout << "my_gls.beta :" << endl;
    print_matrix(my_gls.beta, cout);

    cout << "my_gls.Variance :" << endl;
    print_matrix(my_gls.Variance, cout);

    cout << "my_gls.residual_variance :" << endl;
    print_matrix(my_gls.residual_variance, cout);

    cout << "my_gls.total_variance :" << endl;
    print_matrix(my_gls.total_variance, cout);

    //cout << "SSEc = " << my_residual_square_info.sum() << endl;
    //cout << "SSEr = " << my_residual_square_info_reduced.sum() << endl;
  }
#endif

  const matrix* Variance = &my_gls.Variance;

  if( get_model().get_analysis_options().robust_variance )
    Variance = &my_gls.SandwichVariance;

  size_t df = pair_count() - valid_parameter_count();
  size_t p  = valid_parameter_count();

  if( get_model().is_x_linked() )
  {
    size_t k = 0;
    for( size_t j = 0; j < p; ++j )
    {
      for( size_t x = 0; x < x_type_count; ++x )
      {
        compute_parameter_result(j, my_gls.beta(k,0), (*Variance)(k,k), df, my_x_types[x]);
        ++k;
      }
    }

    for( size_t x = 0; x < vi; ++x, ++k )
    {
      compute_intercept_result(x, my_gls.beta(k,0), (*Variance)(k,k), my_x_intercept_types[x]);
    }
  }
  else
  {
    size_t k = 0;
    for( size_t j = 0; j < p; ++j )
    {
      compute_parameter_result(j, my_gls.beta(k,0), (*Variance)(k,k), df);
      ++k;
    }

    for( ; k < p+vi; ++k )
    {
      compute_intercept_result(k-p, my_gls.beta(k,0), (*Variance)(k,k));
    }
  }

  my_final_As = final_As;
  my_final_Ws = final_Ws;
  my_final_Ys = final_Ys;
  my_final_Rs = final_Rs;

  my_valid_regression = true;

  return;
}

void
TraitRegression::do_regress_univariate_sub(map<const sib_cluster*, matrix>& final_As,
                                           map<const sib_cluster*, matrix>& final_Ws,
                                           map<const sib_cluster*, matrix>& final_Ys,
                                           map<const sib_cluster*, matrix>& final_Rs,
                                           size_t m, double fv, double hv)
{
  bool iterate_1 = true; // using residuals
  bool iterate_2 = true; // using max_diff between previous & current estimates

  size_t max_iterations = max( (size_t)1, get_model().get_analysis_options().max_iterations);

  vector<pair<double, double> > previous_estimates;

  for( size_t i = 0; iterate_2 && i < max_iterations; ++i )
  {
    my_gls.reset();

    if( !(my_use_fsib && my_use_hsib) )
    {
      fv = 1.0;
      hv = 1.0;
    }

#if 0
  if( my_use_fsib && my_use_hsib )
      cout << "iter " << i << " fv = " << fv << " hv = " << hv << " ";
#endif

    regress_markers(my_gls, i != 0, final_As, final_Ws, final_Ys, fv, hv);

    if( !my_gls.beta || !my_gls.beta.rows() )
    {
      errors << priority(error) << "Not enough pairs." << endl;
      return;
    }

    if( my_gls.beta.rows() != m )
    {
#if 0
      errors << priority(error)
             << "Internal error: Parameter count wrong in univeriate regression."
             << "  (b=" << my_gls.beta.rows() <<", p="<<get_model().get_parameter_count()
             << ")" << endl;
#endif
      return;
    }

    if( !my_gls.Variance || !my_gls.Variance.rows() || !my_gls.Variance.cols() )
    {
      errors << priority(error) << "Design matrix singular." << endl;
      return;
    }

    iterate_1 = build_residuals(my_gls, final_Rs, &fv, &hv);

    vector<pair<double, double> > current_estimates(my_gls.beta.rows());

    for( size_t p = 0; p < my_gls.beta.rows(); ++p )
    {
      current_estimates[p].first = my_gls.beta(p,0);

      double variance = my_gls.Variance(p,p);

      if(variance < 0) current_estimates[p].second = QNAN;
      else             current_estimates[p].second = sqrt(variance);

      //cout << p << " : " << my_gls.beta(p,0) << " / " << sqrt(my_gls.Variance(p,p)) << endl;
    }

    if( i )
    {
      double max_diff = 0.;
      for( size_t p = 0; p < my_gls.beta.rows(); ++p )
      {
        double curr = fabs(current_estimates[p].first / current_estimates[p].second);
        double prev = fabs(previous_estimates[p].first / previous_estimates[p].second);

        max_diff = max(max_diff, fabs(prev-curr));

        //cout << "curr = " << curr << ", prev = " << prev << ", max_diff = " << max_diff << endl;
      }

      if( max_diff < 0.1 )
        iterate_2 = false;
    }
 
    previous_estimates = current_estimates;    
  }

  return;
}

void
TraitRegression::check_intercepts()
{
#if 0
  my_reg_results.dump(cout);
#endif
  if( my_use_fsib && my_use_hsib )
  {
    size_t fi = my_reg_results.get_result_count() - 2;
    size_t hi = my_reg_results.get_result_count() - 1;

    double full_intercept = my_reg_results.get_result(fi).estimate.value();
    double half_intercept = my_reg_results.get_result(hi).estimate.value();
#if 0
    cout << "full_intercept = "   << full_intercept
         << ", half_intercept = " << half_intercept << endl;
#endif
    if( full_intercept < 0.0 )
    {
#if 0
      cout << "Now without both intercept..." << endl;
#endif
      // Set both to be 0.
      my_fsib_intercept.resize(0);
      my_hsib_intercept.resize(0);

      my_reg_results_intercept = my_reg_results;
      my_reg_results.clear();

      do_regress_univariate();
#if 0
      cout << "Result with intercept: ";
      my_reg_results_intercept.dump(cout);
      cout << "Result without intercept: ";
      my_reg_results.dump(cout);
#endif
    }
    else
    {
      if( half_intercept < 0.0 )
      {
#if 0
        cout << "Now without half intercept..." << endl;
#endif
        // Set half intercept to be 0.
        my_fsib_intercept.resize(0);
        my_hsib_intercept.resize(0);

        my_fsib_intercept.resize(1, 1.0);
        my_hsib_intercept.resize(1, 0.0);

        my_reg_results_intercept = my_reg_results;
        my_reg_results.clear();

        do_regress_univariate();
#if 0
        cout << "Result with intercept: ";
        my_reg_results_intercept.dump(cout);
        cout << "Result without intercept: ";
        my_reg_results.dump(cout);
#endif
      }
      else if( half_intercept > full_intercept )
      {
#if 0
        cout << "Now with one common intercept..." << endl;
#endif
        // One common intercept.
        my_fsib_intercept.resize(0);
        my_hsib_intercept.resize(0);

        my_fsib_intercept.resize(1, 1.0);
        my_hsib_intercept.resize(1, 1.0);

        my_reg_results_intercept = my_reg_results;
        my_reg_results.clear();

        do_regress_univariate();
#if 0
        cout << "Result with intercept: ";
        my_reg_results_intercept.dump(cout);
        cout << "Result with one common intercept: ";
        my_reg_results.dump(cout);
#endif
        size_t ci = my_reg_results.get_result_count() - 1;
        double intercept = my_reg_results.get_result(ci).estimate.value();
#if 0
        cout << "intercept = " << intercept << endl;
#endif
        if( intercept < 0.0 )
        {
#if 0
          cout << "Now without intercept again..." << endl;
#endif
          my_reg_results.clear();

          my_fsib_intercept.resize(0);
          my_hsib_intercept.resize(0);

          do_regress_univariate();
#if 0
          cout << "Result with intercept: ";
          my_reg_results_intercept.dump(cout);
          cout << "Result without intercept: ";
          my_reg_results.dump(cout);
#endif
        }
      }
    }
  }
  else
  {
    if( get_model().is_x_linked() )
    {
      size_t x_type_count = my_x_types.size();

      vector<double> intercepts;

      for( size_t xt = x_type_count; xt >= 1; --xt )
      {
        size_t xi = my_reg_results.get_result_count() - xt;
        intercepts.push_back(my_reg_results.get_result(xi).estimate.value());
      }

      bool re_do = false;
      for( size_t i = 0; i < intercepts.size(); ++i )
      {
#if 0
        cout << my_x_types[i] << " intercept = " << intercepts[i] << " ";
#endif
        if( intercepts[i] < 0.0 )
        {
          my_x_intercept[my_x_types[i]].resize(0);
          my_x_intercept[my_x_types[i]].resize(x_type_count, 0.0) ;
          re_do = true;
        }
      }
#if 0
      cout << endl;
#endif
      if( re_do )
      {
#if 0
        cout << "Now without intercept..." << endl;
#endif
        vector<bool> keep_cols;
        keep_cols.resize(x_type_count, true);

        for( size_t x = 0; x < x_type_count; ++x )
          if(    my_x_intercept[MM][x] == 0.0
              && my_x_intercept[MF][x] == 0.0
              && my_x_intercept[FF][x] == 0.0 )
            keep_cols[x] = false;

        vector<vector<double> > new_x_intercept;
        new_x_intercept.resize(3);

        for( size_t i = 0; i < 3; ++i )
        {
          vector<pair_type> new_x_intercept_types;

          for( size_t x = 0; x < x_type_count; ++x )
            if( keep_cols[x] )
            {
              new_x_intercept[(pair_type)i].push_back(my_x_intercept[(pair_type)i][x]);
              new_x_intercept_types.push_back(my_x_types[x]);
            }
          my_x_intercept_types = new_x_intercept_types;
        }
            
        my_x_intercept = new_x_intercept;

        my_reg_results_intercept = my_reg_results;
        my_reg_results.clear();

        do_regress_univariate();
#if 0
        cout << "Result with intercept: ";
        my_reg_results_intercept.dump(cout);
        cout << "Result without intercept: ";
        my_reg_results.dump(cout);
#endif
      }
    }
    else
    {
      size_t ci = my_reg_results.get_result_count() - 1;
      double intercept = my_reg_results.get_result(ci).estimate.value();
#if 0
      cout << "intercept = " << intercept << endl;
#endif
      if( intercept < 0.0 )
      {
#if 0
        cout << "Now without intercept..." << endl;
#endif
        my_fsib_intercept.resize(0);
        my_hsib_intercept.resize(0);

        my_reg_results_intercept = my_reg_results;
        my_reg_results.clear();

        do_regress_univariate();
#if 0
        cout << "Result with intercept: ";
        my_reg_results_intercept.dump(cout);
        cout << "Result without intercept: ";
        my_reg_results.dump(cout);
#endif
      }
    }
  }

  if( !my_reg_results_intercept.get_result_count() )
    my_reg_results_intercept = my_reg_results;

  return;
}

void
TraitRegression::do_regress_F_test()
{
  // If # slopes > 1 and value(s) > 0., then re-do the regression without slopes for F-test.
  //
  //cout << my_fix_x_pair[MM] << ", " << my_fix_x_pair[MF] << ", " << my_fix_x_pair[FF] << endl;

  size_t x_type_count = my_x_types.size();

  if( x_type_count < 2 )
    return;

  size_t F_df1 = x_type_count; 

  for( size_t r = 0; r < x_type_count; ++r )
  {
    if( my_gls.beta(r, 0) < 0. )
      F_df1 -= 1;
  }

  //cout << "F_df1 = " << F_df1 << endl;

  if( F_df1 > 1 )
  {
    for( size_t r = 0; r < x_type_count; ++r )
    {
      switch( my_x_types[r] )
      {
        case MM : if( my_gls.beta(r, 0) >= 0. )
                    my_fix_x_pair[MM] = true;
                  break;

        case MF : if( my_gls.beta(r, 0) >= 0. )
                    my_fix_x_pair[MF] = true;
                  break;

        case FF : if( my_gls.beta(r, 0) >= 0. )
                    my_fix_x_pair[FF] = true;
                  break;

        default :
                  break;
      }
    }

    //cout << "re_do_regress.. "
    //     << my_fix_x_pair[MM] << ", " << my_fix_x_pair[MF] << ", " << my_fix_x_pair[FF] << endl;

    size_t vi = valid_intercept_count();
    size_t m  = x_type_count * valid_parameter_count();

    m += vi;

    map<const sib_cluster*, matrix> final_As;
    map<const sib_cluster*, matrix> final_Ws;
    map<const sib_cluster*, matrix> final_Ys;
    map<const sib_cluster*, matrix> final_Rs;

    do_regress_univariate_sub(final_As, final_Ws, final_Ys, final_Rs, m, 1.0, 1.0);

    double SSEc = my_reg_results.get_residual_square_info().sum();
    double SSEr = my_reg_results.get_residual_square_info_reduced().sum();
    double MSEc = my_reg_results.get_residual_variance();

    size_t total_sibs = 0;
    size_t ss = 0;
    for( ; ss < get_sib_clusters().size(); ++ss )
    {
      const sib_cluster& sc = get_sib_clusters()[ss];

      total_sibs += sc.valid_sib_count();
    }

    size_t N = total_sibs - ss;
    size_t p = valid_parameter_count();

    size_t F_df2 = N - p;

    double f_stat = (SSEr - SSEc) / ((double)F_df1 * MSEc);
    double f_pval = 1.0 - F_cdf(F_df1, F_df2, f_stat);

    F_result_type current_F_result;

    current_F_result.set_F_statistic(f_stat);
    current_F_result.set_F_pvalue(f_pval);
    current_F_result.set_df1(F_df1);
    current_F_result.set_df2(F_df2);
    current_F_result.validate();

    my_reg_results.set_F_result(current_F_result);

#if 0
    cout << "my_gls.beta :" << endl;
    print_matrix(my_gls.beta, cout);

    cout << "my_gls.Variance :" << endl;
    print_matrix(my_gls.Variance, cout);

    cout << "my_gls.residual_variance :" << endl;
    print_matrix(my_gls.residual_variance, cout);

    cout << "my_gls.total_variance :" << endl;
    print_matrix(my_gls.total_variance, cout);

    cout << "SSEc = " << SSEc << endl;
    cout << "SSEr = " << SSEr << endl;
    cout << "MSEc = " << MSEc << endl;
    cout << "F_df1 = " << F_df1 << ", F_df2 = " << F_df2 << endl;
    cout << "f stat = " << f_stat << ", f pval = " << f_pval << endl;
#endif
  }

  my_fix_x_pair[MM] = false;
  my_fix_x_pair[MF] = false;
  my_fix_x_pair[FF] = false;

  return;
}

void
TraitRegression::make_sib_cluster()
{
  // Reset the pair counts (see above hack)
  invalidate();
  my_fsib_pair_count = my_hsib_pair_count = 0;

  my_x_pair_count[MM] = 0;
  my_x_pair_count[MF] = 0;
  my_x_pair_count[FF] = 0;
  
  my_sib_clusters.resize(0);

  // Iterate over all sib_clusters
  sibship_cluster_const_iterator c_iter = pairs.sibship_cluster_begin();
  sibship_cluster_const_iterator c_end  = pairs.sibship_cluster_end();

  for( ; c_iter != c_end; ++c_iter )
  {
    if(    my_use_fsib && !my_use_hsib
        && c_iter->full_sibship_map.size() )
    {
      const vector<size_t>& hsib = c_iter->hsib_pairs;

      map< id_pair, vector<size_t> >::const_iterator mi = c_iter->full_sibship_map.begin();
      for( ; mi != c_iter->full_sibship_map.end(); ++mi )
      {
        const vector<size_t>& fsib = mi->second;

        sib_cluster sc(fsib, hsib, pairs, my_use_fsib, my_use_hsib, my_filter, get_model().is_x_linked(), my_use_x_pair[MM], my_use_x_pair[MF], my_use_x_pair[FF]);

        if( !sc.valid_pair_count() ) continue;

        my_sib_clusters.push_back(sc);

        if( sc.valid_fsib_pair_count() )
        {
          my_fsib_pair_count += sc.valid_fsib_pair_count();

          my_x_pair_count[MM] += sc.valid_fsib_mm_pair_count();
          my_x_pair_count[MF] += sc.valid_fsib_mf_pair_count();
          my_x_pair_count[FF] += sc.valid_fsib_ff_pair_count();
        }

        if( sc.valid_hsib_pair_count() )
        {
          my_hsib_pair_count += sc.valid_hsib_pair_count();
        }
      }
    }
    else
    {
      sib_cluster sc(c_iter->fsib_pairs, c_iter->hsib_pairs, pairs, my_use_fsib, my_use_hsib, my_filter, get_model().is_x_linked(), my_use_x_pair[MM], my_use_x_pair[MF], my_use_x_pair[FF]);

      if( !sc.valid_pair_count() ) continue;

      my_sib_clusters.push_back(sc);

      if( sc.valid_fsib_pair_count() )
      {
        my_fsib_pair_count += sc.valid_fsib_pair_count();

        my_x_pair_count[MM] += sc.valid_fsib_mm_pair_count();
        my_x_pair_count[MF] += sc.valid_fsib_mf_pair_count();
        my_x_pair_count[FF] += sc.valid_fsib_ff_pair_count();
      }

      if( sc.valid_hsib_pair_count() )
      {
        my_hsib_pair_count += sc.valid_hsib_pair_count();
      }
    }
  }

#if 0
if( get_model().is_x_linked() )
{
  cout << "my_use_fsib = " << my_use_fsib << ", my_use_hsib = " << my_use_hsib << endl;
  cout << "my_use_mm = " << my_use_x_pair[MM]
       << ", my_use_mf = " << my_use_x_pair[MF]
       << ", my_use_ff = " << my_use_x_pair[FF] << endl;
  cout << "my_fsib_pair_count = " << my_fsib_pair_count << endl;
  cout << "my_hsib_pair_count = " << my_hsib_pair_count << endl;

  cout << "fsib_mm_pair_count = " << my_x_pair_count[MM] << endl;
  cout << "fsib_mf_pair_count = " << my_x_pair_count[MF] << endl;
  cout << "fsib_ff_pair_count = " << my_x_pair_count[FF] << endl;
}
#endif
#if 0
  for( size_t i = 0; i < get_sib_clusters().size(); ++i )
  {
    const sib_cluster& sc = get_sib_clusters()[i];

    cout << endl << "sib cluster " << i << " ";
    sc.dump();

    cout << " All valid pairs:" << endl;
    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];
      cout << "  " << j << " " << this_pair.pair_number()
           << "(" << this_pair.rels().pair.first->name()
           << "," << this_pair.rels().pair.second->name()
           << ")";
    }
    cout << endl;
  }
  cout << endl;
#endif

  return;
}

void
TraitRegression::build()
{
  if( get_model().valid() && regression_method() && built() )
    return;

  // Clears counts
  invalidate();

  if( !regression_method() )
    set_regression_method( getRegressionVariant(*this) );

  if( !regression_method() )
  {
    errors << priority(error) << "Invalid or unsupported regression type: "
           << get_model().get_regression_method_name() << endl;
    return;
  }

  // Clear the simulation vector
  my_simulation_map.resize(0);
  my_group_permutation_vector.resize(0);

  dependent_variable& current_trait = get_model().get_trait();

  current_trait.info.clear();
  current_trait.info.sample_adjustment(1);

  current_trait.mm_pair_info.clear();
  current_trait.mm_pair_info.sample_adjustment(1);

  current_trait.mf_pair_info.clear();
  current_trait.mf_pair_info.sample_adjustment(1);

  current_trait.ff_pair_info.clear();
  current_trait.ff_pair_info.sample_adjustment(1);

  current_trait.trait_used_sibs_info.clear();
  current_trait.trait_used_sibs_info.sample_adjustment(1);

  current_trait.p_correlation[1]           = current_trait.p_correlation[0]             = QNAN;
  current_trait.p_empirical_correlation[1] = current_trait.p_empirical_correlation[0]   = QNAN;
  current_trait.p_sum_diff_correlation[2] = QNAN;
  current_trait.p_sum_diff_correlation[1] = QNAN;
  current_trait.p_sum_diff_correlation[0] = QNAN;
  current_trait.valid = false;

  regression_model::parameter_iterator pi;
  for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
  {
    pi->valid = false;

    pi->info.clear();
    pi->info.sample_adjustment(1);

    pi->mm_pair_info.clear();
    pi->mm_pair_info.sample_adjustment(1);

    pi->mf_pair_info.clear();
    pi->mf_pair_info.sample_adjustment(1);

    pi->ff_pair_info.clear();
    pi->ff_pair_info.sample_adjustment(1);

    for( size_t i = 0; i < pi->covariates.size(); ++i )
      pi->covariates[i].info.clear();
  }

  // Normalize marker effect
  normalize_effects();

  invalidate();

  SimpleCorrelationInfo  trait_fisher_correlations;

  // Iterate over all sib_cluster
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];

      if( !( this_pair.is_fsib_pair() || this_pair.is_hsib_pair() ))
      {
        errors << priority(error)
               << "Unexpected sib pair type in trait regression build" << endl;
      }

      pair<double,double> tvalues = pair_traits(this_pair, current_trait);

      const double& t1 = tvalues.first;
      const double& t2 = tvalues.second;

      trait_fisher_correlations.add(t1, t2);
      trait_fisher_correlations.add(t2, t1);

      for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
        for( size_t i = 0; i < pi->covariates.size(); ++i )
        {
          size_t cc = pi->covariates[i].covariate_index;
          const double cv1 = member_trait( this_pair.rels().pair.first,  cc);
          const double cv2 = member_trait( this_pair.rels().pair.second, cc);
          pi->covariates[i].info += cv1;
          pi->covariates[i].info += cv2;
        }
    }

    sib_set::const_iterator si = sc.get_valid_sibs().begin();
    for( ; si != sc.get_valid_sibs().end(); ++si )
      current_trait.trait_used_sibs_info += member_trait( *si, current_trait.trait_index );
  }
  
  estimate_independent_variable_parameters();

  current_trait.fisher_correlation = trait_fisher_correlations.correlation();

  if( my_use_fsib && !my_use_hsib )
    current_trait.fisher_fsib_correlation = trait_fisher_correlations.correlation();
  else if( !my_use_fsib && my_use_hsib )
    current_trait.fisher_hsib_correlation = trait_fisher_correlations.correlation();

  double trait_variance = current_trait.trait_used_sibs_info.variance();

  if( finite(trait_variance) && trait_variance > 10e-5 )
  {
    current_trait.valid = true;
    ++my_valid_trait_count;
  }
  else
  {
    errors << priority(error) << "No variation in trait '"
           << current_trait.name(pairs)
           << "'.  Trait will be skipped." << endl;

    return;
  }

  for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
  {
    double variance = pi->info.variance();

    if( finite(variance) && variance > 10e-5 )
    {
      pi->valid = true;
      ++my_valid_parameter_count;
    }
    else
      errors << priority(error) << "No variation in parameter '"
             << pi->name(pairs) << " (" << pi->effect_name() << ")"
             << "'.  Parameter will be skipped." << endl;
  }
  
  if( !my_valid_parameter_count )
    return;

  if( my_valid_parameter_count != get_model().get_parameter_count() )
    get_model().clear_invalid_parameters();

  assert( my_valid_parameter_count == get_model().get_parameter_count() );

  // Compute the sib trait correlations
  estimate_trait_sum_diff_correlation();

  if( my_use_fsib && !my_use_hsib )
  {
    estimate_trait_correlation();

    if(    get_model().get_analysis_options().sibship_mean
        || get_model().get_analysis_options().blup_mean )
    {
      estimate_sibship_blup_mean();
    }
  }

  // Allow regression variant to setup its own state so that we can safely
  // access dependent variable information
  regression_method()->build();

  get_model().validate();

  my_built = true;
}

void
TraitRegression::copy_build(const TraitRegression* treg)
{
  my_sib_clusters = treg->my_sib_clusters;
  my_filter = treg->my_filter;

  set_use_pairs( treg->get_use_pairs() );
  set_pair_counts( treg->get_pair_counts() );
  set_use_x_pairs( treg->get_use_x_pairs() );
  set_x_pair_counts( treg->get_x_pair_counts() );
  set_x_types( treg->get_x_types() );

  my_fsib_intercept = treg->my_fsib_intercept;
  my_hsib_intercept = treg->my_hsib_intercept;
  my_x_intercept = treg->my_x_intercept;
  my_x_intercept_types = treg-> my_x_intercept_types;

  if( !regression_method() )
    set_regression_method( getRegressionVariant(*this) );

  if( !regression_method() )
  {
    errors << priority(error) << "Invalid or unsupported regression type: "
           << get_model().get_regression_method_name() << endl;
    return;
  }
#if 0
  cout << "copied from : " << treg->get_model().get_regression_method_name() << endl;
  cout << "copied to   : " << get_model().get_regression_method_name() << endl;
#endif
  get_model().set_trait( treg->get_model().get_trait() );
  get_model().get_trait().valid = true;
  ++my_valid_trait_count;

  for( size_t j = 0; j < treg->get_model().get_parameter_count(); ++j )
  {
    get_model().add_parameter( treg->get_model().get_parameter(j) );
    ++my_valid_parameter_count;
  }

  my_sib_mean_map = treg->my_sib_mean_map;
  my_male_sib_mean_map = treg->my_male_sib_mean_map;
  my_female_sib_mean_map = treg->my_female_sib_mean_map;

  regression_method()->build();

  get_model().validate();

  my_built = true;

  return;
}

void
TraitRegression::normalize_effects()
{
  regression_model::parameter_iterator pi;
  regression_model::parameter_iterator pi2;
  for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
  {
    // Step 1 to normalize effect types -- set all ADDITIVE effects to TOTAL
    for( size_t i = 0; i < pi->markers.size(); ++i )
      if( pi->markers[i].effect == marker_type::ADDITIVE )
        pi->markers[i].effect = marker_type::TOTAL;

    // Step 2 -- sort marker and covariate vectors
    std::sort(pi->covariates.begin(), pi->covariates.end());
    std::sort(pi->markers.begin(),    pi->markers.end());
  }

  for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
  {
    // Step 3 -- find marker to normalize
    bool dom_adjust = false;
    for( size_t i = 0; i < pi->markers.size(); ++i )
      if( pi->markers[i].effect == marker_type::DOMINANCE )
      {
        dom_adjust = true;
        break;
      }

    if( !dom_adjust )
      continue;

    // Step 4: Look for any other parameters that match everything but
    //         TOTAL/ADDITIVE flags
    for( pi2 = get_model().parameter_begin(); pi2 != get_model().parameter_end(); ++pi2 )
    {
      if( &*pi == &*pi2 || !params_dom_equal(*pi, *pi2) )
        continue;

      // Step 5: Flip all TOTAL markers to ADDITIVE if the reference
      //         parameter includes the corresponding DOMINANCE term.
      for( size_t i = 0; i < pi->markers.size(); ++i )
        if(  pi->markers[i].effect  == marker_type::DOMINANCE &&
             pi2->markers[i].effect == marker_type::TOTAL)
          pi2->markers[i].effect = marker_type::ADDITIVE;
    }
  }
}

void 
TraitRegression::estimate_independent_variable_parameters()
{
  regression_model::parameter_iterator pi;

  // Iterate over all sib_clusters
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];

      for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
      {
        double par_val = parameter_value(this_pair, *pi, false, false);
        pi->info += par_val;

        if( is_mm_pair(this_pair.rels().pair) )
          pi->mm_pair_info += par_val;
        else if( is_mf_pair(this_pair.rels().pair) )
          pi->mf_pair_info += par_val;
        else if( is_ff_pair(this_pair.rels().pair) )
          pi->ff_pair_info += par_val;
        else if( get_model().is_x_linked() )
          errors << priority(error)
                 << "Unexpected pair type found in estimate_independent_variable_parameters." << endl;

      }
    }
  }

#if 0
  size_t i = 0;
  for(pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi, ++i)
  {
    cout << "param " << i
         << " : mean = " << pi->info.mean() << ", var = " << pi->info.variance() << endl
         << "  mm_mean = " << pi->mm_pair_info.mean() << ", var = " << pi->mm_pair_info.variance() << endl
         << "  mf_mean = " << pi->mf_pair_info.mean() << ", var = " << pi->mf_pair_info.variance() << endl
         << "  ff_mean = " << pi->ff_pair_info.mean() << ", var = " << pi->ff_pair_info.variance() << endl;
  }
#endif

  return;
}
void
TraitRegression::estimate_trait_sum_diff_correlation()
{
  // Trait sums and differences
  KahanAdder<double>    ts;
  KahanAdder<double>    td;
  KahanAdder<double>   ts2;
  KahanAdder<double>   td2;
  KahanAdder<double>  tstd;

  KahanAdder<double> t1_info;
  KahanAdder<double> t2_info;

  size_t      total_n = 0;

  dependent_variable& current_trait = get_model().get_trait();

  // Iterate over all sib_cluster to compute pass 1.
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    total_n += sc.valid_pair_count();

    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];

      double mean = current_trait.get_mean();

      pair<double,double> tvalues = pair_traits(this_pair, current_trait);
      const double& t1 = tvalues.first  - mean;
      const double& t2 = tvalues.second - mean;

      double std  = -(t1 - t2)*(t1 - t2);
      double sts  =  (t1 + t2)*(t1 + t2);

      ts   += sts;
      td   += std;
      ts2  += sts*sts;
      td2  += std*std;
      tstd += sts*std;

      t1_info += tvalues.first;
      t2_info += tvalues.second;
    }
  }

#if 0
  KahanAdder<double> cp;

  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    for( size_t j = 0; j < sc.valid_pair_count(); ++j )
    {
      const sib_pair this_pair = sc[j];

      double t1_m = t1_info / (size_t)total_n;
      double t2_m = t2_info / (size_t)total_n;
      double mu   = (t1_m + t2_m) / 2.;

      pair<double,double> tvalues = pair_traits( this_pair, current_trait);
      const double& t1 = tvalues.first;
      const double& t2 = tvalues.second;

      cp += ((t1 - mu)*(t2 - mu));
    }
  }
#endif

  // Compute the correlation between the squared sum and difference
  td   /= pair_count();
  ts   /= pair_count();
  ts2  /= pair_count();
  td2  /= pair_count();
  tstd /= pair_count();

  double var_td   = td2  - td*td;
  double var_ts   = ts2  - ts*ts;
  double cov_tstd = tstd - ts*td;

  double corr = cov_tstd / sqrt( var_ts * var_td );

#if 0
  cout << "var_td   = " << var_td   << endl;
  cout << "var_ts   = " << var_ts   << endl;
  cout << "cov_tstd = " << cov_tstd << endl;
  cout << "info_corr = " << corr     << endl;

//  cout << "corr_new = " << cp / (total_n) << endl;
#endif

  current_trait.sum_diff_correlation = corr;

  return;
}

void
TraitRegression::estimate_trait_correlation()
{
  //
  //                 sum (v_i*(xb_i - xb_v)^2) - v_c/(N-k) *  sum(SS_i)
  // Let gamma_v =   ---------------------------------------------------
  //                 sum (v_i*(xb_i - xb_v)^2) + v_c(v_0-1)
  //                                             ---------- * sum(SS_i)
  //                                              N-K
  //
  //                 a     -         v_c/(N-K) * b
  //             =  -------------------------------
  //                 a     + v_c*(v_0-1)/(N-K) * b
  //
  //
  // Let       v = sum(v_i)
  //        xb_i = mean(trait for current sibship)  [computed in pass 1.1]
  //        SS_i = sum of mean corrected squares    [computed in pass 1]
  //        xb_v = 1/v * sum(xb_i * v_i)
  //         v_c = sum( v_i - v_i^2/v )/n_i
  //         v_0 = (v - sum(v_i^2/v)/v_c
  //
  //  where    k = # of sibships
  //         n_i = # of sibs in the i-th sibship
  //
  // case 1:  v1_i = 1                            [assumed]
  //  then    v1   = k                            [known from build]
  //         xb1_v = 1/k * sum( xb_i )            [computed in pass 1]
  //          v1_c = sum(1/n_i - 1/(n_i*k))       [computed in pass 1]
  //          v1_0 = (k - 1) / v_c                [known after pass 1]
  //            a1 = sum((xb_i - xb1_v)^2)        [computed in pass 2]
  //             b = sum(SS_i)                    [computed in pass 1]
  //
  //                   a1 -          v1_c/(N-k) * b
  //        gamma1 = --------------------------------
  //                   a1 + v1_c*(v1_0-1)/(N-k) * b
  //
  // case 2:  v2_i = s_i*(s_i-1)                       [computed in pass 1]
  //  then      v2 = sum( s_i*(s_i-1) )                [computed  ]
  //         xb2_v = 1/v * sum(xb_i * v_i)             [computed in pass 1]
  //          v2_c = (v_i - v_i^2/v)/n_i               [computed in pass 1]
  //          v2_0 = (v - 1/v * sum(v_i^2)) / v_c      [computed in pass 1]
  //            a2 = sum( v2_i*(xb_i - xb2_v)^2 )      [computed in pass 2]
  //             b = ...                               [as above]
  //
  //                   a2 -          v1_c/(N-k) * b
  //        gamma2 = --------------------------------
  //                   a2 + v2_c*(v2_0-1)/(N-k) * b
  //
  // Gamma = gamma2 / ( 1 + gamma2 - gamma1 )
  //

  KahanAdder<double>    a1;
  KahanAdder<double>    a2;
  KahanAdder<double>     b;
  KahanAdder<double>  v1_0;
  KahanAdder<double>  v2_0;
  KahanAdder<double>  v1_c;
  KahanAdder<double>  v2_c;
  KahanAdder<double>  xb_i;
  KahanAdder<double> xb1_v;
  KahanAdder<double> xb2_v;
  KahanAdder<double> v2_i2;   // Sum of v2_i^2

  // Set up trait and marker filter and create sibship iterators
  dependent_variable& current_trait = get_model().get_trait();

  sib_set::const_iterator si;
  double v1  = 0;
  double  k  = 0;
  double v2  = 0;
  double  N  = 0.0;

  // Init values for pass 1
  xb1_v = 0.0;
  xb2_v = 0.0;
  v1_c  = 0.0;
  v2_c  = 0.0;
  v2_i2 = 0.0;
  b     = 0.0;
  a1    = 0.0;
  a2    = 0.0;

  // Iterate over all sibships to compute pass 1.
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    const double s_i  = sc.get_valid_sibs().size();
    const double v1_i = 1;
    const double v2_i = s_i*(s_i-1.0);           // Case 2

    // Pass 1.0: Compute sib set

    // Accumulate other things
    N  += s_i;
    v1 += v1_i;
    v2 += v2_i;
    k  += 1;

    // Pass 1.1: Compute mean trait value for sibship xb_i
    xb_i = 0.0;

    for( si = sc.get_valid_sibs().begin(); si != sc.get_valid_sibs().end(); ++si )
      xb_i += member_trait( *si, current_trait.trait_index );

    xb_i /= s_i;

    // Compute things based on xb_i and v_i
    xb1_v += xb_i * v1_i;
    xb2_v += xb_i * v2_i;
    v2_i2 += v2_i*v2_i;
  }

  xb1_v /= v1;
  if(xb2_v != 0 && v2 != 0)
    xb2_v /= v2;

  // Iterate over all sibships to compute pass 2
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    const double s_i  = sc.get_valid_sibs().size();
    const double v1_i = 1;                       // Case 1
    const double v2_i = s_i*(s_i-1.0);           // Case 2

    // Pass 1.0: Compute sib set

    // Pass 2.1: Compute mean trait value for sibship xb_i
    xb_i = 0.0;
    for( si = sc.get_valid_sibs().begin(); si != sc.get_valid_sibs().end(); ++si )
      xb_i += member_trait( *si, current_trait.trait_index );

    xb_i /= s_i;

    // Compute things based on xb_i
    for( si = sc.get_valid_sibs().begin(); si != sc.get_valid_sibs().end(); ++si )
    {
      double tv = member_trait( *si, current_trait.trait_index );
      b += (tv - xb_i)*(tv - xb_i);
    }

    a1 += v1_i*(xb_i - xb1_v)*(xb_i - xb1_v);
    a2 += v2_i*(xb_i - xb2_v)*(xb_i - xb2_v);
    v1_c += v1_i*(1.0 - v1_i/v1)/s_i;
    if(v2 != 0)
      v2_c += v2_i*(1.0 - v2_i/v2)/s_i;
  }

  // Pass 2.3: Compute final elements
  if( v1_c != 0 )
    v1_0 = (v1 -      1.0) / v1_c;
  else
    v1_0 = 0;

  if(v2 != 0 && v2_c != 0)
    v2_0 = (v2 - v2_i2/v2) / v2_c;
  else
    v2_0 = 0;

  // Compute the sib intra-class correlation

  current_trait.correlation = QNAN;

  double gamma11 = a1 - v1_c/(N-k)*b;
  double gamma12 = a1 + v1_c*(v1_0 - 1.0)/(N-k)*b;

  double gamma1 = 0;
  if(gamma12 != 0)
    gamma1 = gamma11 / gamma12;

  double gamma21 = a2 - v2_c/(N-k)*b;
  double gamma22 = a2 + v2_c*(v2_0 - 1.0)/(N-k)*b;

  double gamma2 = 0;
  if(gamma22 != 0)
     gamma2 = gamma21/gamma22;

  if( fabs(1.0 + gamma2 - gamma1) > 10*EPS )
  {
    double r = gamma2 / (1.0 + gamma2 - gamma1);
    double V = current_trait.trait_all_sibs_info.variance();

    current_trait.correlation = r;

    double rr = max(0., r);

    current_trait.var_b = V * rr;
    current_trait.var_r = V * (1.0 - rr);      
  }

#if 0
  cout << "N = " << N << endl;
  cout << "k = " << k << endl;
  cout << "v1 = " << v1 << endl;
  cout << "v2 = " << v2 << endl;
  cout << "gamma11 = " << gamma11 << endl;
  cout << "gamma12 = " << gamma12 << endl;
  cout << "gamma21 = " << gamma21 << endl;
  cout << "gamma22 = " << gamma22 << endl;
  cout << "b = " << b << endl;
  cout << "r     = " << current_trait.correlation << endl;
  cout << "Var_b = " << current_trait.var_b << endl;
  cout << "Var_r = " << current_trait.var_r << endl;
#endif

  return;
}

void
TraitRegression::estimate_sibship_blup_mean()
{
  my_sib_mean_map.clear();
  my_male_sib_mean_map.clear();
  my_female_sib_mean_map.clear();

  size_t t = get_model().get_trait().trait_index;

  // Set up trait and marker filter and create sibship iterators
  sib_set::const_iterator si;

  // Iterate over all sibships to compute pass 1.
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    SampleInfo  sibship_mean_info;
    SampleInfo  sibship_male_mean_info;
    SampleInfo  sibship_female_mean_info;

    for( si = sc.get_valid_sibs().begin(); si != sc.get_valid_sibs().end(); ++si )
    {
      sibship_mean_info += member_trait(*si, t);

      if( (*si)->is_male() )
        sibship_male_mean_info += member_trait(*si, t);
      else
        sibship_female_mean_info += member_trait(*si, t);
    }

    double xp_mean = sibship_mean_info.mean();
    double x_mean  = get_model().get_trait().trait_all_sibs_info.mean();
    double var_b   = get_model().get_trait().var_b;
    double var_r   = get_model().get_trait().var_r;
    double np      = sc.get_valid_sibs().size();

    double w = var_b / (var_b + (var_r / np));
    double blup_mean = w * xp_mean + (1. - w) * x_mean;

    double male_xp_mean = sibship_male_mean_info.mean();
    double male_x_mean  = get_model().get_trait().trait_male_sibs_info.mean();
    double male_blup_mean = w * male_xp_mean + (1. - w) * male_x_mean;

    double female_xp_mean = sibship_female_mean_info.mean();
    double female_x_mean  = get_model().get_trait().trait_female_sibs_info.mean();
    double female_blup_mean = w * female_xp_mean + (1. - w) * female_x_mean;

#if 0
  cout << "xp_mean = " << xp_mean << endl;
  cout << "x_mean  = " << x_mean  << endl;
  cout << "var_b   = " << var_b   << endl;
  cout << "var_r   = " << var_r   << endl;
  cout << "np      = " << np      << endl;
  cout << "w       = " << w       << endl;
  cout << "blup_mean = " << blup_mean << endl;
#endif

#if 0
  cout << "w       = " << w       << endl;
  cout << "male_xp_mean   = " << male_xp_mean << endl;
  cout << "male_x_mean    = " << male_x_mean  << endl;
  cout << "male_blup_mean = " << male_blup_mean << endl;
  cout << "female_xp_mean   = " << female_xp_mean << endl;
  cout << "female_x_mean    = " << female_x_mean  << endl;
  cout << "female_blup_mean = " << female_blup_mean << endl;
#endif

    for( si = sc.get_valid_sibs().begin(); si != sc.get_valid_sibs().end(); ++si )
    {
      my_sib_mean_map[*si] = make_pair(xp_mean, blup_mean);

      if( (*si)->is_male() )
        my_male_sib_mean_map[*si] = make_pair(male_xp_mean, male_blup_mean);
      else
        my_female_sib_mean_map[*si] = make_pair(female_xp_mean, female_blup_mean);
    }
  }

#if 0
  map< const MPED::member_base*, pair< double, double > >::const_iterator mi;

  for( mi = my_sib_mean_map.begin(); mi!= my_sib_mean_map.end(); ++mi )
  {
    cout << mi->first->name() << " : " << mi->second.first << "," << mi->second.second << endl;
  }

  for( mi = my_male_sib_mean_map.begin(); mi!= my_male_sib_mean_map.end(); ++mi )
  {
    cout << mi->first->name() << " : " << mi->second.first << "," << mi->second.second << endl;
  }

  for( mi = my_female_sib_mean_map.begin(); mi!= my_female_sib_mean_map.end(); ++mi )
  {
    cout << mi->first->name() << " : " << mi->second.first << "," << mi->second.second << endl;
  }
#endif

  return;
}

void
TraitRegression::regress_markers(GLS3& gls, bool use_empirical_correlations,
                                 map<const sib_cluster*, matrix>& As,
                                 map<const sib_cluster*, matrix>& Ws,
                                 map<const sib_cluster*, matrix>& Ys,
                                 double fv, double hv)
{
  map<const sib_cluster*, matrix> temp_As;
  map<const sib_cluster*, matrix> temp_Ws;
  map<const sib_cluster*, matrix> temp_Ys;

  matrix W;                    // weight matrix
  matrix A;                    // design matrix
  matrix y;                    // trait vector

  const trait_parameter& t = get_model().get_trait();

  bool  standardize = get_model().get_analysis_options().standardize_parameter;

  gls.reset();

  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    // Compute basis matrices

    weight_status_type status = BESTW;

    weight_matrix(sc, W, t, use_empirical_correlations, status, fv, hv);

    if( get_model().is_x_linked() )
      design_matrix_x(sc, A, standardize);
    else
      design_matrix(sc, A, standardize);

    trait_vector(sc, y, t, use_empirical_correlations, true);

#if 0
  cout << "A =" << endl;
  print_matrix(A, cout);
  cout << "y =" << endl;
  print_matrix(y, cout);
  cout << "W =" << endl;
  print_matrix(W, cout);
#endif

    gls.add_block(y, W, A, status);

    temp_As[&sc] = A;
    temp_Ws[&sc] = W;
    temp_Ys[&sc] = y;
  }

  As = temp_As;
  Ws = temp_Ws;
  Ys = temp_Ys;

  if( gls.observation_count - (valid_intercept_count()+valid_parameter_count()) <= 0 )
  {
    errors << priority(error) << "Not enough valid pairs." << endl;
    return;
  }

  gls.compute();

#if 0
  cout << "my_gls.beta :" << endl;
  print_matrix(my_gls.beta, cout);

  cout << "my_gls.Variance :" << endl;
  print_matrix(my_gls.Variance, cout);

  cout << "my_gls.residual_variance :" << endl;
  print_matrix(my_gls.residual_variance, cout);

  cout << "my_gls.total_variance :" << endl;
  print_matrix(my_gls.total_variance, cout);
#endif

  return;
}

bool
TraitRegression::build_residuals(GLS3& gls, map<const sib_cluster*, matrix>& Rs,
                                 double* fv, double* hv)
{
  if( !gls.beta || !gls.Variance )
    return false;

  map<const sib_cluster*, matrix> temp_Rs;

  matrix A;                    // design matrix
  matrix y;                    // trait vector
  matrix residuals;            // y-Ab

  size_t vi = valid_intercept_count();
  size_t x_type_count = my_x_types.size();

  size_t m = x_type_count * valid_parameter_count();

  m += vi;

  size_t fsib_pair_count = 0;
  size_t hsib_pair_count = 0;

  const trait_parameter& t = get_model().get_trait();
  bool         standardize = get_model().get_analysis_options().standardize_parameter;

  SimpleCorrelationInfo          pool_correlation;
  SimpleCorrelationInfo          pool_fsib_correlation;
  SimpleCorrelationInfo          pool_hsib_correlation;
  SimpleCorrelationInfo          pool_fh_correlation;

  vector<SimpleCorrelationInfo>  p_correlation(2);
  vector<SimpleCorrelationInfo>  p_fsib_correlation(2);
  vector<SimpleCorrelationInfo>  p_hsib_correlation(2);
  vector<SimpleCorrelationInfo>  p_fh_correlation(2);

  double fsib_residual_sum = 0.0;
  double hsib_residual_sum = 0.0;

  SampleInfo residual_info;
  SampleInfo residual_square_info;

  residual_info.clear();
  residual_info.sample_adjustment(m);

  residual_square_info.clear();
  residual_square_info.sample_adjustment(m);

  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    fsib_pair_count += sc.valid_fsib_pair_count();
    hsib_pair_count += sc.valid_hsib_pair_count();

    if( get_model().is_x_linked() )
      design_matrix_x(sc, A, standardize);
    else
      design_matrix(sc, A, standardize);

    trait_vector(sc, y, t, true, true);

    gls.build_residuals(y, A, residuals);

#if 0
    if( my_use_fsib && my_use_hsib )
    {
      cout << "A =" << endl;
      print_matrix(A, cout);
      cout << "y =" << endl;
      print_matrix(y, cout);
      cout << "residuals =" << endl;
      print_matrix(residuals, cout);
    }
#endif

    for( size_t i = 0; i < sc.valid_pair_count(); ++i )
    {
      residual_info        += residuals(i,0);
      residual_square_info += (residuals(i,0) * residuals(i,0));

      if( sc[i].is_fsib_pair() )
        fsib_residual_sum += (residuals(i,0) * residuals(i,0));
      else if( sc[i].is_hsib_pair() )
        hsib_residual_sum += (residuals(i,0) * residuals(i,0));
      else
        errors << priority(error)
               << "Unexpected pair type found in residuals calculation." << endl;

      for( size_t j = i+1; j < sc.valid_pair_count(); ++j )
      {
        switch( sibs_shared( sc[i].rels().pair, sc[j].rels().pair ) )
        {
          case 2:  break;  // Correlation is 1

          case 1:  // One sib in common

            if( sc[i].is_fsib_pair() && sc[j].is_fsib_pair() )
            {
              p_fsib_correlation[1].add(residuals(i,0), residuals(j,0));
              p_fsib_correlation[1].add(residuals(j,0), residuals(i,0));

              pool_fsib_correlation.add(residuals(i,0), residuals(j,0));
              pool_fsib_correlation.add(residuals(j,0), residuals(i,0));
            }
            else if( sc[i].is_hsib_pair() && sc[j].is_hsib_pair() )
            {
              p_hsib_correlation[1].add(residuals(i,0), residuals(j,0));
              p_hsib_correlation[1].add(residuals(j,0), residuals(i,0));

              pool_hsib_correlation.add(residuals(i,0), residuals(j,0));
              pool_hsib_correlation.add(residuals(j,0), residuals(i,0));
            }
            else
            {
              p_fh_correlation[1].add(residuals(i,0), residuals(j,0));
              p_fh_correlation[1].add(residuals(j,0), residuals(i,0));

              pool_fh_correlation.add(residuals(i,0), residuals(j,0));
              pool_fh_correlation.add(residuals(j,0), residuals(i,0));
            }

            p_correlation[1].add(residuals(i,0), residuals(j,0));
            p_correlation[1].add(residuals(j,0), residuals(i,0));

            pool_correlation.add(residuals(i,0), residuals(j,0));
            pool_correlation.add(residuals(j,0), residuals(i,0));

            break;

          case 0:  // No sibs in common

            if( sc[i].is_fsib_pair() && sc[j].is_fsib_pair() )
            {
              p_fsib_correlation[0].add(residuals(i,0), residuals(j,0));
              p_fsib_correlation[0].add(residuals(j,0), residuals(i,0));

              pool_fsib_correlation.add(residuals(i,0), residuals(j,0));
              pool_fsib_correlation.add(residuals(j,0), residuals(i,0));
            }
            else if( sc[i].is_hsib_pair() && sc[j].is_hsib_pair() )
            {
              p_hsib_correlation[0].add(residuals(i,0), residuals(j,0));
              p_hsib_correlation[0].add(residuals(j,0), residuals(i,0));

              pool_hsib_correlation.add(residuals(i,0), residuals(j,0));
              pool_hsib_correlation.add(residuals(j,0), residuals(i,0));
            }
            else
            {
              p_fh_correlation[0].add(residuals(i,0), residuals(j,0));
              p_fh_correlation[0].add(residuals(j,0), residuals(i,0));

              pool_fh_correlation.add(residuals(i,0), residuals(j,0));
              pool_fh_correlation.add(residuals(j,0), residuals(i,0));
            }

            p_correlation[0].add(residuals(i,0), residuals(j,0));
            p_correlation[0].add(residuals(j,0), residuals(i,0));

            pool_correlation.add(residuals(i,0), residuals(j,0));
            pool_correlation.add(residuals(j,0), residuals(i,0));

            break;

          default: // Unrelated pairs
            break;
        }
      }
    }

    temp_Rs[&sc] = residuals;
  }

  Rs = temp_Rs;
  
  size_t fsib_df = fsib_pair_count - valid_parameter_count();
  size_t hsib_df = hsib_pair_count - valid_parameter_count();

  double fsib_residual_var = fsib_residual_sum / fsib_df;
  double hsib_residual_var = hsib_residual_sum / hsib_df;

  *fv = fsib_residual_var;
  *hv = hsib_residual_var;

#if 0
  cout << endl
       << "residual: mean = " << my_residual_info.mean()
       << ", variance = "     << my_residual_info.variance()
       << ", skewness = "     << my_residual_info.skewness()
       << ", kurtosis = "    << my_residual_info.kurtosis()-3.0
       << endl;

  if( my_use_fsib && my_use_hsib )
  {
    cout << endl
         << "fsib_residual_sum = " << fsib_residual_sum << endl
         << "hsib_residual_sum = " << hsib_residual_sum << endl
         << "fsib_df           = " << fsib_df << endl
         << "hsib_df           = " << hsib_df << endl
         << "fsib_residual_var = " << fsib_residual_var << endl
         << "hsib_residual_var = " << hsib_residual_var << endl
         << endl;

    cout << "fsib residual: p0 = " << p_fsib_correlation[0].correlation()
         <<              ", p1 = " << p_fsib_correlation[1].correlation()
         <<        "pooled_cor = " << pool_fsib_correlation.correlation() << endl;
    cout << "hsib residual: p0 = " << p_hsib_correlation[0].correlation()
         <<              ", p1 = " << p_hsib_correlation[1].correlation()
         <<        "pooled_cor = " << pool_hsib_correlation.correlation() << endl;
    cout << "fh   residual: p0 = " << p_fh_correlation[0].correlation()
         <<              ", p1 = " << p_fh_correlation[1].correlation()
         <<        "pooled_cor = " << pool_fh_correlation.correlation() << endl;
    cout << "tot  residual: p0 = " << p_correlation[0].correlation()
         <<              ", p1 = " << p_correlation[1].correlation()
         <<        "pooled_cor = " << pool_correlation.correlation() << endl;
  }
#endif

  double p0 = p_correlation[0].correlation();
  double p1 = p_correlation[1].correlation();

  double restricted_p0 = max(0.,p0);
  double restricted_p1 = max(0.,max(p0,p1));

  double pooled_cor = pool_correlation.correlation();

  if( p0 > p1 && get_model().get_analysis_options().pool_correlation )
    restricted_p0 = restricted_p1 = max(0.,pooled_cor);

  trait_parameter& trait = get_model().get_trait();

  double last_p0 = trait.p_empirical_correlation[0];
  double last_p1 = trait.p_empirical_correlation[1];

  trait.p_empirical_correlation[0] = restricted_p0;
  trait.p_empirical_correlation[1] = restricted_p1;

  double p0_ff = p_fsib_correlation[0].correlation();
  double p1_ff = p_fsib_correlation[1].correlation();

  double restricted_p0_ff = max(0.,p0_ff);
  double restricted_p1_ff = max(0.,max(p0_ff,p1_ff));

  double pooled_ff_cor = pool_fsib_correlation.correlation();

  if( p0_ff > p1_ff && get_model().get_analysis_options().pool_correlation )
  {
    restricted_p0_ff = restricted_p1_ff = max(0.,pooled_ff_cor);
  }

  double p0_hh = p_hsib_correlation[0].correlation();
  double p1_hh = p_hsib_correlation[1].correlation();

  double restricted_p0_hh = max(0.,p0_hh);
  double restricted_p1_hh = max(0.,max(p0_hh,p1_hh));

  double pooled_hh_cor = pool_fsib_correlation.correlation();

  if( p0_hh > p1_hh && get_model().get_analysis_options().pool_correlation )
  {
    restricted_p0_hh = restricted_p1_hh = max(0.,pooled_hh_cor);
  }

  double p0_fh = p_fh_correlation[0].correlation();
  double p1_fh = p_fh_correlation[1].correlation();

  double restricted_p0_fh = max(0.,p0_fh);
  double restricted_p1_fh = max(0.,max(p0_fh,p1_fh));

  double pooled_fh_cor = pool_fh_correlation.correlation();

  if( p0_fh > p1_fh && get_model().get_analysis_options().pool_correlation )
  {
    restricted_p0_fh = restricted_p1_fh = max(0.,pooled_fh_cor);
  }

#if 0
  if( my_use_fsib && my_use_hsib )
  {
    cout << "last residual: p0 = " << last_p0 << ", p1 = " << last_p1 << endl;
    cout << "p0    = " << restricted_p0    << ", p1    = " << restricted_p1
         << ", pooled cor = " << pooled_cor << endl;
    cout << "p0_f  = " << restricted_p0_ff << ", p1_f  = " << restricted_p1_ff
         << ", pooled cor = " << pooled_ff_cor << endl;
    cout << "p0_h  = " << restricted_p0_hh << ", p1_h  = " << restricted_p1_hh
         << ", pooled cor = " << pooled_hh_cor << endl;
    cout << "p0_fh = " << restricted_p0_fh << ", p1_fh = " << restricted_p1_fh
         << ", pooled cor = " << pooled_fh_cor << endl;
  }
#endif

  trait.p_fsib_empirical_correlation[0] = restricted_p0_ff;
  trait.p_fsib_empirical_correlation[1] = restricted_p1_ff;
  trait.p_hsib_empirical_correlation[0] = restricted_p0_hh;
  trait.p_hsib_empirical_correlation[1] = restricted_p1_hh;
  trait.p_fh_empirical_correlation[0]   = restricted_p0_fh;
  trait.p_fh_empirical_correlation[1]   = restricted_p1_fh;

  my_reg_results.set_residual_info(residual_info);

  if( my_use_fsib && !my_use_hsib )
    my_reg_results.set_full_residual_info(residual_info);
  else if( !my_use_fsib && my_use_hsib )
    my_reg_results.set_half_residual_info(residual_info);

  if(    get_model().is_x_linked()
      && ( my_fix_x_pair[MM] || my_fix_x_pair[MF] || my_fix_x_pair[FF] ) )
    my_reg_results.set_residual_square_info_reduced(residual_square_info);
  else
    my_reg_results.set_residual_square_info(residual_square_info);


  double delta = 0;
  double max_delta = 0.02;

  if( finite(p0) )
  {
    if( finite(last_p0) && fabs(last_p0) > max_delta/2 )
      delta += fabs( (last_p0 - p0) / last_p0 );
    else
      delta += fabs(p0);
  }

  if( finite(p1) )
  {
    if( finite(last_p1) && fabs(last_p1) > max_delta/2 )
      delta += fabs( (last_p1 - p1) / last_p1 );
    else
      delta += fabs(p1);
  }

#if 0
  cout << "delta = " << delta << ", max_delta = " << max_delta << endl;
#endif

  if( delta > max_delta )
    return true;

  return false;
}

void
TraitRegression::weight_matrix(const sib_cluster&     ship,
                               matrix&                W,
                               const trait_parameter& t,
                               bool                   use_empirical_correlations,
                               weight_status_type&    status,
                               double                 fsib_var,
                               double                 hsib_var) const
{
  int n = ship.valid_pair_count();
  W.clear();

  if( !n )
  {
    W.setstate(matrix::badbit);
    return;
  }

  if( regression_method() )
    regression_method()->weight_matrix(ship, W, t, use_empirical_correlations, status, fsib_var, hsib_var);
  else
    W.setstate(matrix::badbit);

  return;
}

void
TraitRegression::trait_vector(const sib_cluster&     ship,
                              matrix&                y,
                              const trait_parameter& t,
                              bool                   use_empirical_correlations,
                              bool                   center) const
{
  if( regression_method() )
  {
    regression_method()->trait_vector(ship, y, t, use_empirical_correlations, center);
  }
  else
    errors << priority(error) << "No trait vector available!" << endl;

  return;
}

void
TraitRegression::design_matrix(const sib_cluster& ship,
                               matrix&            A,
                               bool               standardize) const
{
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    A.setstate(matrix::failbit);
    return;
  }

  size_t p  = get_model().get_parameter_count();
  size_t vi = valid_intercept_count();

  A.resize_nofill(pair_count, vi+p);

  for( size_t i = 0; i < pair_count; ++i )
  {
    size_t k = 0;
    for( size_t j = 0; j < p; ++j )
    {
      const independent_variable& param = get_model().get_parameter(j);

      double s = parameter_value(ship[i], param, false, true, false);

      if(standardize)
      {
        // FIXME:  The mean does not include the centered covariates!
        double mean = param.info.mean();
        double var  = param.info.variance();

        if(!finite(var) || var <= 0 )
          var = 1.0;
        s = (s-mean)/sqrt(var);
      }

      A(i, k) = s;
      ++k;
    }

    for( size_t c = 0; c < vi; ++c, ++k )
    {
      if( ship[i].is_fsib_pair() )
        A(i, k) = my_fsib_intercept[c];
      else
        A(i, k) = my_hsib_intercept[c];
    }
  }

  return;
}

void
TraitRegression::design_matrix_x(const sib_cluster& ship,
                                 matrix&            A,
                                 bool               standardize) const
{
  size_t pair_count = ship.valid_fsib_pair_count();

  if( !pair_count )
  {
    A.setstate(matrix::failbit);
    return;
  }

  size_t x_type_count = my_x_types.size();

  size_t p  = get_model().get_parameter_count();
  size_t vi = valid_intercept_count();

  A.resize_nofill(pair_count, x_type_count*p + vi);

  matrix no_fix_A    = eye<double>(3);
  matrix fix_A       = eye<double>(3);

  if( my_fix_x_pair[MM] )
    fix_A(0,0) = 0.0;

  if( my_fix_x_pair[MF] )
    fix_A(1,1) = 0.0;

  if( my_fix_x_pair[FF] )
    fix_A(2,2) = 0.0;

  for( size_t i = 0; i < pair_count; ++i )
  {
    pair_type pt = MM;
    size_t k = 0;
    for( size_t j = 0; j < p; ++j )
    {
      const independent_variable& param = get_model().get_parameter(j);

      bool fix_effect = false;
      if( param.markers.size() && !param.covariates.size() )
        fix_effect = true;

      double s = parameter_value(ship[i], param, false, true, false);

      pt = MM;
      double mean = param.mm_pair_info.mean();
      double var  = param.mm_pair_info.variance();
      size_t c = 0;

      if( ship[i].is_fsib_pair() && ship[i].is_mf_pair() )
      {
        pt = MF;
        mean = param.mf_pair_info.mean();
        var  = param.mf_pair_info.variance();

        if( !my_use_x_pair[MM] )
          c = 1;
      }
      else if( ship[i].is_fsib_pair() && ship[i].is_ff_pair() )
      {
        pt = FF;
        mean = param.ff_pair_info.mean();
        var  = param.ff_pair_info.variance();

        if( x_type_count == 2 )
          c = 1;
        else if ( x_type_count == 1 )
          c = 2;
      }

      if( !finite(var) || var <= 0 )
        var = 1.0;

      if( standardize )
      {
        s = (s-mean) / sqrt(var);
      }

      for( size_t st = 0; st < x_type_count; ++st, ++c )
      {
        if( fix_effect )
          A(i, k) = fix_A(pt, c) * s;
        else
          A(i, k) = no_fix_A(pt, c) * s;
        ++k;
      }
    }

    for( size_t v = 0; v < vi; ++v )
    {
      A(i, k) = my_x_intercept[pt][v];
      ++k;
    }
  }

  return;
}

void
TraitRegression::compute_robust_variance(GLS3& gls, bool use_empirical_correlations)
{
  const trait_parameter& t = get_model().get_trait();

  size_t vi = valid_intercept_count();
  size_t x_type_count = my_x_types.size();

  size_t m = x_type_count * valid_parameter_count();

  m += vi;

  bool standardize = get_model().get_analysis_options().standardize_parameter;

  if( gls.observation_count - m <= 0 )
  {
    return;
  }

  matrix W;                    // weight matrix
  matrix A;                    // design matrix
  matrix y;                    // trait vector

  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    // Compute basis matrices

    weight_status_type status = BESTW;

    if( get_model().is_x_linked() )
      design_matrix_x(sc, A, standardize);
    else
      design_matrix(sc, A, standardize);

    trait_vector(sc, y, t, use_empirical_correlations, true);
    weight_matrix(sc, W, t, use_empirical_correlations, status, 1.0, 1.0);

    gls.update_robust_variance_block(y, W, A, status);
  }

  gls.compute_robust_variance();
}

void
TraitRegression::compute_parameter_result(size_t j, double beta, double ss, size_t df, pair_type pt)
{
#if 0
  cout << "TraitRegression::compute_parameter_result.. " << endl
       << "j = " << j
       << ", beta = " << beta
       << ", ss = " << ss
       << ", df = " << df
       << ", pt = " << pt << endl;
#endif

  const independent_variable &param = get_model().get_parameter(j);

  reg_result current_result;

  current_result.set_parameter(param, j);
  current_result.set_pair_type(pt);
  current_result.estimate.set_degrees_of_freedom(df);
  current_result.estimate.set_expected_value(0.0);

  if( param.covariates.size() )
  {
    current_result.estimate.set_test_direction(parameter_estimate::dTWOSIDED);
  }
  else
  {
    current_result.estimate.set_test_direction(parameter_estimate::dRIGHTSIDE);
  }

  bool standardize = get_model().get_analysis_options().standardize_parameter;

  if( standardize )
  {
    double std = param.info.variance();

    if( get_model().is_x_linked() )
    {
      if( pt == MM )
        std = param.mm_pair_info.variance();
      else if( pt == MF )
        std = param.mf_pair_info.variance();
      else
        std = param.ff_pair_info.variance();
    }

    if( finite(std) )
    {
      beta /= sqrt(std);
      ss   /= std;
    }
  }

  current_result.estimate.set_value(beta);
  current_result.estimate.set_variance(ss);

  my_reg_results.add_result(current_result);

  if( my_use_fsib && !my_use_hsib )
  {
    my_reg_results.add_full_result(current_result);
  }
  else if( !my_use_fsib && my_use_hsib )
  {
    my_reg_results.add_half_result(current_result);
  }

#if 0
  cout << "After adjust: beta = " << beta << ", ss = " << ss << endl;
  //my_reg_results.dump(cout);
#endif

  return;
}

void
TraitRegression::compute_intercept_result(size_t r, double alpha, double ss, pair_type pt)
{
#if 0
  cout << "compute_intercept_result()..." << endl
       << "r = " << r
       << ", alpha = " << alpha
       << ", ss = " << ss
       << ", pt = " << pt << endl;
#endif

  bool standardize = get_model().get_analysis_options().standardize_parameter;

  alpha += get_model().get_trait().info.mean();

  if( standardize )
  {
    size_t p  = valid_parameter_count();
    size_t vt = valid_intercept_count();

    if( get_model().is_x_linked() )
      vt = my_x_types.size();

    for( size_t j = 0; j < p; ++j )
    {
      const independent_variable& param = get_model().get_parameter(j);

      double p_mean = param.info.mean();

      if( get_model().is_x_linked() )
      {
        if( pt == MM )
          p_mean = param.mm_pair_info.mean();
        else if( pt == MF )
          p_mean = param.mf_pair_info.mean();
        else
          p_mean = param.ff_pair_info.mean();

        double param_est = my_reg_results.get_result((j*vt)+pt).estimate.value();
        //cout << "p_mean = " << p_mean << ", less = " << param_est << endl;

        alpha -= param_est * p_mean;
      }
      else
      {
        double param_est = my_reg_results.get_result(j).estimate.value();
        //cout << "p_mean = " << p_mean << ", less = " << param_est << endl;

        alpha -= param_est * p_mean;
      }
    }
  }

  reg_result current_result;

  current_result.set_intercept();
  current_result.set_pair_type(pt);
  current_result.estimate.set_value(alpha);
  current_result.estimate.set_variance(ss);

  my_reg_results.add_result(current_result);

  if( my_use_fsib && !my_use_hsib )
  {
    my_reg_results.add_full_result(current_result);
  }
  else if( !my_use_fsib && my_use_hsib )
  {
    my_reg_results.add_half_result(current_result);
  }

#if 0
  cout << "After adjust:  alpha = " << alpha << ", ss = " << ss << endl;
//  my_reg_results.dump(cout);
#endif

  return;
}

double
TraitRegression::marker_sharing(const sib_pair&    pr,
                                const marker_type& param,
                                bool               permute,
                                bool               normalize) const
{
  size_t m = param.marker_index;
  bool dom = (param.effect == marker_type::DOMINANCE);

  double t = std::numeric_limits<double>::quiet_NaN();

  sib_pair pair = pr;

  if(permute)
    pair = permutation_pair(pr);

  if(dom)
    t = pair.prob_share(m,2);
  else
  {
    if( get_model().get_data_options().use_pairs == FSIB )
    {
      double w1 = get_model().get_analysis_options().w1;
      t = pair.weighted_share(m, 0.0, w1, 1.0);
    }
    else
      t = pair.avg_share(m);
  }

  if( pairs.marker_genotype_model(m) == MLOCUS::X_LINKED )
  {
    if( pair.is_fsib_pair() && pair.is_ff_pair() )
      t = pair.prob_share(m,1);
    else
      t = pair.prob_share(m,0);
  }

  return t;
}

double
TraitRegression::marker_sharing(const sib_pair&             pair,
                                const independent_variable& param,
                                bool                        permute,
                                bool                        normalize) const
{
  if( !param.markers.size() )
    return std::numeric_limits<double>::quiet_NaN();

  double m = 1.0;

  for(size_t i = 0; i < param.markers.size(); ++i)
    m *= marker_sharing(pair, param.markers[i], permute, normalize);

  return m;
}

double
TraitRegression::covariate_value(const sib_pair&       pair,
                                 const covariate_type& param,
                                 bool                  center) const
{
  size_t cc = param.covariate_index;
  double cv = numeric_limits<double>::quiet_NaN();

  const relative_pairs& pairs = *pair.pair_data();

  if( param.operation == covariate_type::pair )
  {
    if( cc >= pairs.pair_covariate_count() )
      return cv;

    cv = pairs.get_pair_covariate(pair.pair_number(), cc);

    return cv;
  }

  if( cc >= pairs.trait_count() )
    return cv;

  double t1;
  double t2;

  double mean = param.fixed_mean;

  if( center && !finite(mean) )
    mean = param.info.mean();

  if(!finite(mean))
    mean = 0.0;

  t1 = member_trait( pair.rels().pair.first,  cc) - mean;
  t2 = member_trait( pair.rels().pair.second, cc) - mean;

  switch( param.operation )
  {
    case covariate_type::sum:    cv = t1 + t2;           break;
    case covariate_type::diff:   cv = fabs(t1 - t2);     break;
    case covariate_type::prod:   cv = t1 * t2;           break;
    case covariate_type::single: cv = t1;                break;
    case covariate_type::avg:    cv = (t1 + t2)/2.;      break;
    case covariate_type::none: throw std::exception();
    default: break;
  }

  if(param.power != 1.0)
    cv = pow(cv, param.power);

  return cv;
}

double
TraitRegression::covariate_value(const sib_pair&             pair,
                                 const independent_variable& param,
                                 bool                        center) const
{
  if( !param.covariates.size() )
    return numeric_limits<double>::quiet_NaN();

  double t = 1.0;
  for(size_t i = 0; i < param.covariates.size(); ++i)
    t *= covariate_value(pair, param.covariates[i], center);

  return t;
}

double
TraitRegression::parameter_value(const sib_pair&             pair,
                                 const independent_variable& param,
                                 bool                        center_covariates,
                                 bool                        permute_markers,
                                 bool                        normalize_markers) const
{
  double c =  covariate_value(pair, param, center_covariates);
  double m =   marker_sharing(pair, param, permute_markers, normalize_markers);

  if( !finite(c) && !finite(m) )
    return numeric_limits<double>::quiet_NaN();

  if(!finite(c))
    return m;

  if(!finite(m))
    return c;

  return c*m;
}

//------------------------------------------------------------
// simulation functions
//------------------------------------------------------------

void
TraitRegression::simulate_univariate()
{
  if( !my_valid_regression )
  {
    errors << priority(error)
           << "Simulation requested before regression.  Possible program error."
           << endl;
    return;
  }

  const double ept = get_model().get_pvalue_options().threshold;
  if( finite(ept) && (ept <= 0.0 || ept > 1.0) )
  {
    errors << priority(error)
           << "Empirical p-value threshold out of bounds (" << ept << ").  "
           << "Skipping simulation." << endl;
    return;
  }

  do_simulation_univariate();

  return;
}

void
TraitRegression::do_simulation_univariate()
{
  result_vector& observed_results  = my_reg_results.get_results();
  F_result_type& observed_F_result = my_reg_results.get_F_result();

  for( size_t i = 0; i < observed_results.size(); ++i )
    observed_results[i].estimate.clear_simulation();

  TraitRegression simulation(pairs, errors);

  simulation.invalidate();
  simulation.set_model(get_model());
  simulation.get_model().invalidate();
  simulation.copy_build(this);
  simulation.build_simulation();
  simulation.set_reg_results(get_reg_results());
  simulation.validate();

  if( !simulation.valid() )
    return;

#if DEBUG_SIMULATION
  for( size_t i = 0; i < observed_results.size(); ++i )
  {
    // We only care about empirical p-values for parameters involving
    // marker sharing IBD.  Everything else is a nuisance parameter.
    if(    observed_results[i].type() != result::param
       || !observed_results[i].reg_param().markers.size() )
      continue;

    const double p_hat1 = observed_results[i].estimate.pvalue();
    const double p_hat2 = simulation.get_reg_results().get_result(i).estimate.pvalue();
    cout << "SANITY CHECK: " << p_hat1 << " =? " << p_hat2 << endl;
  }

  const double F_p1 = observed_F_result.get_F_result().get_F_pvalue();
  const double F_p2 = simulation.get_reg_results().get_F_result().get_F_pvalue();
    cout << "SANITY CHECK for F: " <<  F_p1 << " =? " << F_p2 << endl;
#endif

  // Get simulation parameters
  size_t fixed_replicates = get_model().get_pvalue_options().replicates;
  double max_replicates   = get_model().get_pvalue_options().max_replicates;
  double min_replicates   = get_model().get_pvalue_options().min_replicates;

  double alpha            = get_model().get_pvalue_options().confidence;
  double width            = get_model().get_pvalue_options().width;
  double threshold        = get_model().get_pvalue_options().threshold;

#if DEBUG_SIMULATION
  cout << "fixed_replicates = " << fixed_replicates << endl;
  cout << "max_replicates   = " << max_replicates   << endl;
  cout << "min_replicates   = " << min_replicates   << endl;
  cout << "alpha            = " << alpha            << endl;
  cout << "width            = " << width            << endl;
  cout << "threshold        = " << threshold        << endl;
#endif

  if(   !finite(width) || width < 100*EPS 
     || !finite(alpha) || alpha < 100*EPS )
  {
    width = 0.2;
    alpha = 0.95;
  }

  if( !finite(threshold) )
    threshold = 0.05;

  double min_p_hat = 1.0;

  for( size_t i = 0; i < observed_results.size(); ++i )
  {
    // We only care about empirical p-values for parameters involving
    // marker sharing IBD.  Everything else is a nuisance parameter.
    if(    observed_results[i].type() != reg_result::param
       || !observed_results[i].reg_param().markers.size() )
      continue;

    const double p_hat = observed_results[i].estimate.pvalue();
    if( finite(p_hat) )
      min_p_hat = min(p_hat, min_p_hat);
  }

  if( my_reg_results.get_F_result().is_valid() )
  {
    const double p_hat = observed_F_result.get_F_pvalue();
    if( finite(p_hat) )
      min_p_hat = min(p_hat, min_p_hat);
  }
#if DEBUG_SIMULATION
  cout << "min_p_hat = " << min_p_hat << ", ept = " << ept << endl;
#endif

  // Do not simulate if the minimum asymptotic p-value is greater than the
  // threshold value
  if( min_p_hat > threshold )
    return;

  // See binomial trails theory for explaination
  double precision = pow(inv_normal_cdf((1.+alpha)/2.)/width, 2);

  if( fixed_replicates )
  {
    precision = QNAN;
    min_replicates = fixed_replicates;
    max_replicates = fixed_replicates;
  }

#if DEBUG_SIMULATION
  cout << "fixed_replicates = " << fixed_replicates << endl;
  cout << "max_replicates   = " << max_replicates   << endl;
  cout << "min_replicates   = " << min_replicates   << endl;
  cout << "alpha            = " << alpha            << endl;
  cout << "width            = " << width            << endl;
  cout << "threshold        = " << threshold        << endl;
  cout << "precision        = " << precision        << endl; 
#endif

  do_simulation_univariate_replicates(simulation, precision,
                                      (size_t)min_replicates, (size_t)max_replicates);

  return;
}

void
TraitRegression::build_simulation()
{
#if DEBUG_SIMULATION
  cout << "TraitRegression::build_simulation()..." << endl;
#endif

  assert( get_model().valid() && built() );

  if( !get_model().valid() || !built() )
    return;

  if( get_model().get_pvalue_options().seed )
    randomizer.mt.reseed(get_model().get_pvalue_options().seed);

  // Clear the simulation vectors
  my_simulation_map.resize(0);
  my_group_permutation_vector.resize(0);

  simulation_map_type            sim_map(pairs.pair_count());
  group_permutation_vector_type  group_vector;

  // Iterate over all sib_cluster
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    size_t n = sc.valid_pair_count();

    if( group_vector.size() <= n )
      group_vector.resize(n+1);

    permutation_vector_type permutation_vector(n);

    size_t sibship_number = group_vector[n].size();

    for( size_t j = 0; j < n; ++j )
    {
      const sib_pair this_pair = sc[j];

      size_t m = this_pair.pair_number();

      // We initialize the data structures necessary to permute the IBD
      // sharing.  The simulation vector is a map from total pairs->valid
      // pairs.  Since valid pairs is a dense ordering, we construct the
      // initial permutation vector as an identity function.  To randomize
      // the IBD sharing, all we must do is permute the permutation vector.

      sim_map[m] = pair_index(n, sibship_number, j);
      permutation_vector[j] = m;
    }

    assert(permutation_vector.size() == n);
    group_vector[n].push_back(permutation_vector);
  }

# if DEBUG_SIMULATION
  cout << "dump simulation_map:" << endl;
  for( size_t mi = 0; mi < sim_map.size(); ++mi )
  {
    cout << mi << " : pair_index = (" << sim_map[mi].sibship_size
                               << "," << sim_map[mi].sibship_number
                               << "," << sim_map[mi].sib_number << ")" << endl;
  }
  cout << endl;

  cout << "dump group_permutation_vector : size = " << group_vector.size() << endl;
  for( size_t gi = 0; gi < group_vector.size(); ++gi )
  {
    cout << gi << endl;

    sibship_permutation_vector_type svector = group_vector[gi];

    for( size_t si = 0; si < svector.size(); ++si )
    {
      cout << si << " :";

      permutation_vector_type pvector = svector[si];
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        cout << " " << pvector[pi];
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;
#endif

  // Move the simulation vector in to place
  my_simulation_map.swap(sim_map);
  my_group_permutation_vector.swap(group_vector);
}

void
TraitRegression::do_simulation_univariate_replicates(TraitRegression& simulation,
                                                     double precision,
                                                     size_t min_replicates,
                                                     size_t max_replicates)
{
#if DEBUG_SIMULATION
  cout << endl
       << "*** vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv ***"
       << endl;
#endif

  result_vector& observed_results  = my_reg_results.get_results();
  F_result_type& observed_F_result = my_reg_results.get_F_result();

  for( size_t i = 0; i < max_replicates; ++i )
  {
    if( simulation.get_model().is_x_linked() )
      simulation.randomize_marker_data_x();
    else
      simulation.randomize_marker_data();

    simulation.simulating_build();

    if( !simulation.valid() )
      continue;

    assert( simulation.my_group_permutation_vector.size() );

    simulation.get_reg_results().clear();
    simulation.do_regress_univariate();

    if( !simulation.valid() )
      continue;

    if( simulation.get_model().is_x_linked() && observed_F_result.is_valid() )
      simulation.do_regress_F_test();

    const result_vector& empirical_results = simulation.get_reg_results().get_results();
    const F_result_type& empirical_F_result = simulation.get_reg_results().get_F_result();

    double min_p = 1.0;
    for( size_t j = 0; j < empirical_results.size(); ++j )
    {
      if(    empirical_results[j].type() != reg_result::param
         || !empirical_results[j].reg_param().markers.size() )
        continue;

      double beta = empirical_results[j].estimate.value();
      double ss   = empirical_results[j].estimate.variance();

      observed_results[j].estimate.add_simulation_result(beta, ss);

      //if(    observed_results[j].type() != reg_result::param
      //   || !observed_results[j].reg_param().markers.size() )
      //  continue;

      const double p = observed_results[j].estimate.empirical_pvalue();

      if( finite(p) )
        min_p = min(p, min_p);

#if DEBUG_SIMULATION
      //cout << empirical_results[j].estimate.pvalue() << endl;
      cout << "  REP " << i << " : "
           <<   "o.t = " << doub2str( observed_results[j].estimate.tvalue(), 7,4)
           << ", o.p = " << doub2str( observed_results[j].estimate.pvalue(), 7,4)
           << ", s.t = " << doub2str(empirical_results[j].estimate.tvalue(), 7,4)
           << ", s.p = " << doub2str(empirical_results[j].estimate.pvalue(), 7,4)
           << ", e.p = " << doub2str( observed_results[j].estimate.empirical_pvalue(), 7,4)
           << ", e.s = " << observed_results[j].estimate.significant_replicates()
           << ", e.e = " << observed_results[j].estimate.same_replicates()
           << ", e.n = " << observed_results[j].estimate.total_replicates()
           << endl;
#endif
    }

    if( observed_F_result.is_valid() && empirical_F_result.is_valid() )
    {
      observed_F_result.add_simulation_result(empirical_F_result.get_F_pvalue());

      const double p = observed_F_result.get_F_empirical_pvalue();

      if( finite(p) )
        min_p = min(p, min_p);

#if DEBUG_SIMULATION
      cout << "  F :"
           << "  o.p = " << doub2str( observed_F_result.get_F_pvalue(), 7,4)
           << ", s.p = " << doub2str(empirical_F_result.get_F_pvalue(), 7,4)
           << ", e.p = " << doub2str( observed_F_result.get_F_empirical_pvalue(), 7,4)
           << ", e.s = " << observed_F_result.get_significant_replicates()
           << ", e.s = " << observed_F_result.get_same_replicates()
           << ", e.n = " << observed_F_result.get_total_replicates()
           << endl;
#endif
    }

    double m = i*min_p/(1-min_p);

#if DEBUG_SIMULATION
    cout << "n=" << i+1 << ", p=" << min_p << ", m=" << m << ", prec=" << precision << endl;
#endif

    // Break if our results are precise enough
    if( i>min_replicates && finite(m) && m > precision )
      break; 
  }

#if DEBUG_SIMULATION
  cout << "*** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ ***"
       << endl;
#endif
}

void
TraitRegression::randomize_marker_data()
{
  size_t n = my_group_permutation_vector.size();

  for(size_t i = 0; i < n; ++i)
  {
    size_t sibship_count = my_group_permutation_vector[i].size();
    random_shuffle(my_group_permutation_vector[i].begin(),
                   my_group_permutation_vector[i].end(),
                   randomizer);

    for(size_t j = 0; j < sibship_count; ++j)
      random_shuffle(my_group_permutation_vector[i][j].begin(),
                     my_group_permutation_vector[i][j].end(),
                     randomizer);
  }

#if DEBUG_SIMULATION
  cout << "dump group_permutation_vector : size = " << my_group_permutation_vector.size() << endl;
  for( size_t gi = 0; gi < my_group_permutation_vector.size(); ++gi )
  {
    cout << gi << endl;

    sibship_permutation_vector_type svector = my_group_permutation_vector[gi];

    for( size_t si = 0; si < svector.size(); ++si )
    {
      cout << si << " :";

      permutation_vector_type pvector = svector[si];
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        cout << " " << pvector[pi];
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;
#endif

  return;
}

#if 0
void
dump_pair_data(const vector<size_t>& pairs)
{
  for( size_t i = 0; i < pairs.size(); ++i )
    cout << " " << pairs[i];
  cout << endl;
}
#endif

void
TraitRegression::randomize_marker_data_x()
{
#if 0
  cout << "TraitRegression::randomize_marker_data_x()..." << endl;
#endif

  group_permutation_vector_type group_vector;

  for( size_t gi = 0; gi < my_group_permutation_vector.size(); ++gi )
  {
    LinearDatamap<string, size_t> x_pair_type_map;

    sibship_permutation_vector_type svector = my_group_permutation_vector[gi];

    for( size_t si = 0; si < svector.size(); ++si )
    {
      permutation_vector_type pvector = svector[si];

      size_t mm_pair_cnt = 0;
      size_t mf_pair_cnt = 0;
      size_t ff_pair_cnt = 0;
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        sib_pair sp = sib_pair(pvector[pi], &pairs);

        if( sp.is_fsib_pair() && sp.is_mm_pair() )
          ++mm_pair_cnt;
        else if( sp.is_fsib_pair() && sp.is_mf_pair() )
          ++mf_pair_cnt;
        else if( sp.is_fsib_pair() && sp.is_ff_pair() )
          ++ff_pair_cnt;
      }

      stringstream s;
      s << mm_pair_cnt;
      s << mf_pair_cnt;
      s << ff_pair_cnt;

      x_pair_type_map.insert(make_pair(s.str(), si));
    }

    vector<size_t> org_svector_index;

    LinearDatamap<string, size_t>::type_iterator t;
    LinearDatamap<string, size_t>::value_iterator v;

    for( t = x_pair_type_map.type_begin(); t != x_pair_type_map.type_end(); ++t )
    {
      for( v = t->second.begin(); v != t->second.end(); ++v )
      {
        org_svector_index.push_back(*v);
      }
    }

#if 0
    cout << "dump before shuffle : " << endl;

    for( size_t si = 0; si < org_svector_index.size(); ++si )
      cout << org_svector_index[si] << " ";
    cout << endl;
#endif

    for( t = x_pair_type_map.type_begin(); t != x_pair_type_map.type_end(); ++t )
    {
      random_shuffle(t->second.begin(), t->second.end(), randomizer);
    }

#if 0
    cout << "dump after shuffle : " << endl;
    for( t = x_pair_type_map.type_begin(); t != x_pair_type_map.type_end(); ++t )
    {
      cout << t->first << ": ";
      for( v = t->second.begin(); v != t->second.end(); ++v )
        cout << *v << " ";
      cout << endl;
    }
    cout << endl;
#endif

    vector<size_t> new_svector_index(org_svector_index.size());
    size_t new_si = 0;
    for( t = x_pair_type_map.type_begin(); t != x_pair_type_map.type_end(); ++t )
    {
      for( v = t->second.begin(); v != t->second.end(); ++v )
        new_svector_index[org_svector_index[new_si++]] = *v;
    }

#if 0
    for( size_t si = 0; si < new_svector_index.size(); ++si )
      cout << new_svector_index[si] << " ";
    cout << endl;
#endif

    sibship_permutation_vector_type new_svector;

    for( size_t si = 0; si < new_svector_index.size(); ++si )
    {
      new_svector.push_back(svector[new_svector_index[si]]);
    }

    for( size_t si = 0; si < new_svector.size(); ++si )
    {
      permutation_vector_type pvector = new_svector[si];

      vector<size_t> mm_pairs;
      vector<size_t> mf_pairs;
      vector<size_t> ff_pairs;
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        sib_pair sp = sib_pair(pvector[pi], &pairs);

        if( sp.is_fsib_pair() && sp.is_mm_pair() )
          mm_pairs.push_back(pvector[pi]);
        else if( sp.is_fsib_pair() && sp.is_mf_pair() )
          mf_pairs.push_back(pvector[pi]);
        else if( sp.is_fsib_pair() && sp.is_ff_pair() )
          ff_pairs.push_back(pvector[pi]);
      }

#if 0
      cout << "dump before shuffle : " << endl;
      dump_pair_data(mm_pairs);
      dump_pair_data(mf_pairs);
      dump_pair_data(ff_pairs);
#endif

      random_shuffle(mm_pairs.begin(), mm_pairs.end(), randomizer);
      random_shuffle(mf_pairs.begin(), mf_pairs.end(), randomizer);
      random_shuffle(ff_pairs.begin(), ff_pairs.end(), randomizer);

#if 0
      cout << "dump after shuffle : " << endl;
      dump_pair_data(mm_pairs);
      dump_pair_data(mf_pairs);
      dump_pair_data(ff_pairs);
#endif

      size_t mmi = 0;
      size_t mfi = 0;
      size_t ffi = 0;

      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        sib_pair sp = sib_pair(svector[si][pi], &pairs);

        if( sp.is_fsib_pair() && sp.is_mm_pair() )
          new_svector[si][pi] = mm_pairs[mmi++];
        else if( sp.is_fsib_pair() && sp.is_mf_pair() )
          new_svector[si][pi] = mf_pairs[mfi++];
        else if( sp.is_fsib_pair() && sp.is_ff_pair() )
          new_svector[si][pi] = ff_pairs[ffi++];
      }
    }

    group_vector.push_back(new_svector);
  }

  my_group_permutation_vector.swap(group_vector);

/*
  for( size_t gi = 0; gi < my_group_permutation_vector.size(); ++gi )
  {
    vector<size_t> mm_pairs;
    vector<size_t> mf_pairs;
    vector<size_t> ff_pairs;

    sibship_permutation_vector_type svector = my_group_permutation_vector[gi];

    for( size_t si = 0; si < svector.size(); ++si )
    {
      permutation_vector_type pvector = svector[si];
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        sib_pair sp = sib_pair(pvector[pi], &pairs);

        if( sp.is_fsib_pair() && sp.is_mm_pair() )
          mm_pairs.push_back(pvector[pi]);
        else if( sp.is_fsib_pair() && sp.is_mf_pair() )
          mf_pairs.push_back(pvector[pi]);
        else if( sp.is_fsib_pair() && sp.is_ff_pair() )
          ff_pairs.push_back(pvector[pi]);
      }
    }

    cout << "dump before shuffle : " << endl;
    dump_pair_data(mm_pairs);
    dump_pair_data(mf_pairs);
    dump_pair_data(ff_pairs);

    random_shuffle(mm_pairs.begin(), mm_pairs.end(), randomizer);
    random_shuffle(mf_pairs.begin(), mf_pairs.end(), randomizer);
    random_shuffle(ff_pairs.begin(), ff_pairs.end(), randomizer);

    cout << "dump after shuffle : " << endl;
    dump_pair_data(mm_pairs);
    dump_pair_data(mf_pairs);
    dump_pair_data(ff_pairs);

    size_t mmi = 0;
    size_t mfi = 0;
    size_t ffi = 0;

    for( size_t si = 0; si < my_group_permutation_vector[gi].size(); ++si )
    {
      for( size_t pi = 0; pi < my_group_permutation_vector[gi][si].size(); ++pi )
      {
        sib_pair sp = sib_pair(my_group_permutation_vector[gi][si][pi], &pairs);

        if( sp.is_fsib_pair() && sp.is_mm_pair() )
          my_group_permutation_vector[gi][si][pi] = mm_pairs[mmi++];
        else if( sp.is_fsib_pair() && sp.is_mf_pair() )
          my_group_permutation_vector[gi][si][pi] = mf_pairs[mfi++];
        else if( sp.is_fsib_pair() && sp.is_ff_pair() )
          my_group_permutation_vector[gi][si][pi] = ff_pairs[ffi++];
      }
    }
  }
*/

#if 0
  cout << "dump group_permutation_vector : size = " << my_group_permutation_vector.size() << endl;
  for( size_t gi = 0; gi < my_group_permutation_vector.size(); ++gi )
  {
    cout << gi << endl;

    sibship_permutation_vector_type svector = my_group_permutation_vector[gi];

    for( size_t si = 0; si < svector.size(); ++si )
    {
      cout << si << " :";

      permutation_vector_type pvector = svector[si];
      
      for( size_t pi = 0; pi < pvector.size(); ++pi )
      {
        cout << " " << pvector[pi];
      }
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;
#endif

  return;
}

void
TraitRegression::simulating_build()
{
  if( !get_model().valid() || !regression_method() || !built() )
    return;

  dependent_variable& current_trait = get_model().get_trait();

  current_trait.p_empirical_correlation[1] = current_trait.p_empirical_correlation[0] = QNAN;

  // Allow regression variant to setup its own state so that we can safely
  // access dependent variable information
  regression_method()->update();
}

void
TraitRegression::copy_simulation_vectors(const TraitRegression& reg)
{
  my_simulation_map           = reg.my_simulation_map;
  my_group_permutation_vector = reg.my_group_permutation_vector;

  return;
}


/////////////////////////////////////////////////////////////////////////////

void
TraitRegression::print_optional_output()
{
  if( get_model().get_output_options().dump_data )
  {
    string dump_file_name = get_model().get_output_options().dump_data_file + ".qls";

    ofstream dump_file;
    dump_file.open( dump_file_name.c_str() );

    dump_data(dump_file);
  }

  if( get_model().get_output_options().print_design_matrix )
  {
    string design_file_name = get_model().get_output_options().design_matrix_file + ".design";

    ofstream design_file;
    design_file.open( design_file_name.c_str() );

    print_design_matrix_header(design_file);

    map<const sib_cluster*, matrix>::const_iterator mi = my_final_As.begin();
    for( size_t r = 0; mi != my_final_As.end() && r < get_model().get_output_options().design_matrix_rows; ++mi )
    {
      const sib_cluster& sc = *(mi->first);
      const matrix&       A = mi->second;

      print_sibship_data_matrix(sc, A, design_file, r, true);
      r += sc.valid_pair_count();
    }
  }

  if( get_model().get_output_options().print_correl_matrix )
  {
    string correl_file_name = get_model().get_output_options().correl_matrix_file + ".cor";

    ofstream correl_file;
    correl_file.open( correl_file_name.c_str() );

    map<const sib_cluster*, matrix>::const_iterator mi = my_final_Ws.begin();
    for( ; mi != my_final_Ws.end(); ++mi )
    {
      const sib_cluster& sc = *(mi->first);
      const matrix&       W = mi->second;

      my_correlation_matrix_map[make_pair(sc.valid_fsib_pair_count(), sc.valid_hsib_pair_count())] = W;
    }

    print_correlation_matrix(correl_file);
  }
}

double
TraitRegression::compute_b(const matrix& A, const matrix& W, const matrix& y)
{
  matrix Wi, AWi, AWiA, AWiA_i, AWiy, b;

  SAGE::SVD svd;

  svd.inverse_of(W, Wi);

  XTZ(A, Wi, AWi);
  multiply(AWi, A, AWiA);
  svd.inverse_of(AWiA, AWiA_i);
  multiply(AWi, y, AWiy);
  multiply(AWiA_i, AWiy, b);

  return b(0,0);
}
 
void
TraitRegression::dump_data(ostream& out)
{
  if( !get_model().valid() || !built() )
    return;

  const dependent_variable& current_trait = get_model().get_trait();

  out << "Ped\tSibship\tSib1\tSib2";
  out << '\t' << "Trait(" << current_trait.name(pairs) << ")";
  out << '\t' << get_model().get_regression_method_name()
      << "(" << current_trait.name(pairs) << ")";

  //regression_model::parameter_const_iterator pi;
  //for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
  //  out << '\t' << pi->name(pairs) << "_" << pi->effect_name();
  out << "\tQLS";
  out << endl;

  matrix Y, A, W;
  
  // Iterate over all sib_clusters
  for( size_t c = 0; c < get_sib_clusters().size(); ++c )
  {
    const sib_cluster& sc = get_sib_clusters()[c];

    //if( get_model().is_x_linked() )
    //  design_matrix_x(sc, A, true);
    //else
    //  design_matrix(sc, A, true);

    //trait_vector(sc, Y, current_trait, true, true);

    Y = my_final_Ys[&sc];
    A = my_final_As[&sc];
    W = my_final_Ws[&sc];

    double b = compute_b(A, W, Y);

    if( b != 0.0 )
      b /= A.rows();

    if( fabs(b) < 1.0e-11 )
      b = 0.0;

    for( size_t i = 0; i < A.rows(); ++i )
    {
      out << sc[i].rels().pair.first->pedigree()->name()  << '\t'
          << c+1                                     << '\t'
          << sc[i].rels().pair.first->name()              << '\t'
          << sc[i].rels().pair.second->name();

      std::pair<double,double> pt = pair_traits( sc[i], current_trait );

      const RPED::RefTraitInfo &info = pairs.fped_info().trait_info(current_trait.trait_index);

      if( info.type() == RPED::RefTraitInfo::continuous_trait )
        out << '\t' << fp(pt.first, 8, 5) << ',' << fp(pt.second, 8, 5);
      else if( info.type() == RPED::RefTraitInfo::binary_trait )
      {
        switch( int(pt.first + pt.second) )
        {
          case  0: out << "\tcon_unaff";  break;
          case  1: out << "\tdiscordant";               break;
          case  2: out << "\tcon_affect";    break;
          default: out << "\tUNKNOWN";                  break;
        }
      }
      else
        out << "\tUNKNOWN";

      out << '\t';
      if( fabs(Y(i,0)) < 10e-10 )
        out << fp(0.0, 8, 5);
      else
        out << fp(Y(i,0), 8, 5);

      out << '\t' << b;

/* Original dump

      if( get_model().is_x_linked() && get_x_types().size() > 1 )
      {

        for( size_t j = 0; j < A.cols()-valid_intercept_count(); )
        {
          out << '\t';
 
          double A_val = A(i,j);

          if( get_x_types().size() == 3 )
          {
            if( sc[i].is_mf_pair() )
              A_val = A(i,j+1);
            else if( sc[i].is_ff_pair() )
              A_val = A(i,j+2);
          }
          else
          {
            if( sc[i].is_mf_pair() && get_use_x_pairs()[MM] )
              A_val = A(i,j+1);              
            else if( sc[i].is_ff_pair() )
              A_val = A(i,j+1);
          }

          if( fabs(A_val) < 10e-10 )
            out << fp(0.0, 8, 5);
          else
            out << fp(A_val, 8 ,5);

          j += get_x_types().size();
        }
      }
      else
      {
        for( size_t j = 0; j < A.cols()-valid_intercept_count(); ++j )
        {
          out << '\t';
          if( fabs(A(i,j)) < 10e-10)
            out << fp(0.0, 8, 5);
          else
            out << fp(A(i,j), 8 ,5);
        }
      }
*/
      out << endl;
    }
  }
}

void
TraitRegression::print_design_matrix_header(ostream& out) const
{
  if( !get_model().valid() || !built() )
    return;

  regression_model::parameter_const_iterator pi;

  out << "Ped\tSib1\tSib2";

  if( valid_intercept_count() <= 1 )
  {
    for( size_t i = 0; i < valid_intercept_count(); ++i )
      out << "\tIntercept";

    for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
    {
      out << '\t' << pi->name(pairs);

      if( pi->markers.size() == 1 && pi->covariates.size() == 0 )
        out << "_" << pi->short_effect_name();
    }
  }
  else
  {
    if( get_model().is_x_linked() )
    {
      for( size_t i = 0; i < get_x_types().size(); ++i )
      {
        out << "\tIntercept";
        if( get_x_types()[i] == MM )
          out << "_BB";
        else if( get_x_types()[i] == MF )
          out << "_BS";
        else if( get_x_types()[i] == FF )
          out << "_SS";
      }

      for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
      {
        for( size_t i = 0; i < get_x_types().size(); ++i )
        {
          out << '\t' << pi->name(pairs);

          if( pi->markers.size() == 1 && pi->covariates.size() == 0 )
            out << "_" << pi->short_effect_name();

          if( get_x_types()[i] == MM )
            out << "_BB";
          else if( get_x_types()[i] == MF )
            out << "_BS";
          else if( get_x_types()[i] == FF )
            out << "_SS";
        }
      }
    }
    else
    {
      out << "\tIntercept_FSIB\tIntercept_HSIB";

      for( pi = get_model().parameter_begin(); pi != get_model().parameter_end(); ++pi )
      {
        out << '\t' << pi->name(pairs);

        if( pi->markers.size() == 1 && pi->covariates.size() == 0 )
          out << "_" << pi->short_effect_name();
      }
    }
  }

  out << endl;
}

void
TraitRegression::print_correlation_matrix_header(ostream& out) const
{
  if( !get_model().valid() || !built() )
    return;

  out << "=================================================================================" << endl
      << "  Sibship Specific Correlation Matrices for Dependent Variable" << endl
      << endl
      << "    Note: The first 3 columns are the legend for the row." << endl
      << "          Column legends are the same as row legends in order within a sibship." << endl
      << "=================================================================================" << endl
      << endl;

  out << "Ped\tSib1\tSib2"
      << endl;
}

void
TraitRegression::print_sibship_data_matrix(const sib_cluster& ship,
                                           const matrix&      A,
                                           ostream&           out,
                                           size_t             r,
                                           bool               rows_restricted) const
{
  if( !get_model().valid() || !built() )
    return;

  size_t rows_to_print = get_model().get_output_options().design_matrix_rows;

  for( size_t i = 0; i < A.rows() && r < rows_to_print; ++i )
  {
    if( i )
      out << std::endl;

    out << ship[i].rels().pair.first->pedigree()->name()  << '\t'
        << ship[i].rels().pair.first->name()              << '\t'
        << ship[i].rels().pair.second->name();

    for( size_t j = A.cols()-valid_intercept_count(); j < A.cols(); ++j )
    {
      out << '\t';

      if( fabs(A(i,j)) < 10e-15)
        out << "0.000000";
      else
        out << fp(A(i,j),8,6);
    }
    for( size_t j = 0; j < A.cols()-valid_intercept_count(); ++j )
    {
      out << '\t';

      if( fabs(A(i,j)) < 10e-15)
        out << "0.000000";
      else
        out << fp(A(i,j),8,6);
    }  

    if( rows_restricted )
      ++r;
  }
  out << std::endl;
}

void
TraitRegression::print_correlation_matrix(ostream& out) const
{
  if( !get_model().valid() || !built() )
    return;

  correlation_matrix_map::const_iterator mi = my_correlation_matrix_map.begin();
  for( ; mi != my_correlation_matrix_map.end(); ++mi )
  {
    const matrix& w = mi->second;

    for( int i = 0; i < (int)w.rows(); ++i )
    {
      if( i )
        out << endl;

      for( int j = 0; j < (int)w.cols(); ++j )
      {
        if ( j )
          out << '\t';

        if( fabs(w(i,j)) < 10e-15)
          out << "0.000000";
        else
          out << fp(w(i,j),8,6);
      }  
    }
    out << endl << endl;
  }
}

#ifdef BROKEN_STL_RANDOM_SHUFFLE
#define random_shuffle myrandom_shuffle

template<typename RandomAccessIter, typename RandomNumberGenerator>
inline void
random_shuffle(RandomAccessIter first, RandomAccessIter last,
               RandomNumberGenerator& rand)
{
  if (first == last)
    return;
  for (RandomAccessIter i = first + 1; i != last; ++i)
    iter_swap(i, first + rand((i - first) + 1));
}
#endif

} //end of namespace SIBPAL
} //end of namespace SAGE
