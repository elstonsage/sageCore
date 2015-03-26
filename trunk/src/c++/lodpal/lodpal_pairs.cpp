//****************************************************************************
//* File:      lodpal_pairs.cpp                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   1. Initial implementation                         yjs         *
//*            2. bad_sib_pair storage & func added.             yjs Mar. 01 *
//*            3. rel_pair_map storage & func added.             yjs Mar. 01 *
//*            4. seperated from ARPTest.                        yjs Jul. 01 *
//*            5. weight added.                                  yjs Jan. 02 *
//*            6. X-linkage added.                               yjs May. 02 *
//*            7. new cov option "avg" added.                    yjs may. 02 *
//*                                                                          *
//* Notes:     This file implements lodpal_pairs class.                      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "numerics/matfunctions.h"
#include "lodpal/lodpal_pairs.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of lodpal_pairs (non-Inline)                   //
////////////////////////////////////////////////////////////////////////////

lodpal_pairs::lodpal_pairs(RelativePairs& p, lodpal_parameters& par, cerrorstream& err)
             : pairs(p), params(par), errors(err)
{
  invalidate_build_pairs_info();
}

lodpal_pairs::lodpal_pairs(lodpal_pairs& lp)
             : pairs(lp.relative_pairs()), params(lp.parameters()), errors(lp.errs())
{
  my_filter              = lp.filter();
  my_pair_count          = lp.pair_count();
  my_fsib_pair_count     = lp.fsib_pair_count();
  my_hsib_pair_count     = lp.hsib_pair_count();
  my_pairs_info          = lp.pairs_info();
  my_pairs_map           = lp.pairs_map();
  my_prior_x_ibd_vector  = lp.prior_x_ibds();
  my_built_pairs_info    = lp.built_pairs_info();
  my_built_pairs_info_x  = lp.built_pairs_info_x();
  my_re_built_pairs_info = lp.re_built_pairs_info();
  my_re_built_covariate  = lp.re_built_covariate();
  my_max_ped_name        = lp.max_ped_name_size();
  my_max_ind_name        = lp.max_ind_name_size();
  my_removed_fsib_pairs  = lp.removed_fsib_pairs();
  my_removed_hsib_pairs  = lp.removed_hsib_pairs();
  my_removed_other_pairs = lp.removed_other_pairs();

  my_mm_pair_count          = lp.mm_pair_count();
  my_mm_fsib_pair_count     = lp.mm_fsib_pair_count();
  my_mfm_hsib_pair_count    = lp.mfm_hsib_pair_count();
  my_mm_invalid_pair_count  = lp.mm_invalid_pair_count();

  my_mf_pair_count          = lp.mf_pair_count();
  my_mf_fsib_pair_count     = lp.mf_fsib_pair_count();
  my_mff_hsib_pair_count    = lp.mff_hsib_pair_count();
  my_mf_invalid_pair_count  = lp.mf_invalid_pair_count();

  my_ff_pair_count          = lp.ff_pair_count();
  my_ff_fsib_pair_count     = lp.ff_fsib_pair_count();
  my_fff_hsib_pair_count    = lp.fff_hsib_pair_count();
  my_ff_invalid_pair_count  = lp.ff_invalid_pair_count();

  my_parent_of_origin_allowed = lp.parent_of_origin_allowed();

  my_drp_ibd_info = lp.my_drp_ibd_info;
}

lodpal_pairs::~lodpal_pairs()
{}

void
lodpal_pairs::invalidate_build_pairs_info()
{
  my_pair_count      = 0;
  my_fsib_pair_count = 0;
  my_hsib_pair_count = 0;

  my_mm_pair_count         = 0;
  my_mm_fsib_pair_count    = 0;
  my_mfm_hsib_pair_count   = 0;
  my_mm_invalid_pair_count = 0;

  my_mf_pair_count         = 0;
  my_mf_fsib_pair_count    = 0;
  my_mff_hsib_pair_count   = 0;
  my_mf_invalid_pair_count = 0;

  my_ff_pair_count         = 0;
  my_ff_fsib_pair_count    = 0;
  my_fff_hsib_pair_count   = 0;
  my_ff_invalid_pair_count = 0;

  my_built_pairs_info    = false;
  my_built_pairs_info_x  = false;
  my_re_built_pairs_info = false;

  my_parent_of_origin_allowed = false;

  my_removed_fsib_pairs.resize(0);
  my_removed_hsib_pairs.resize(0);
  my_removed_other_pairs.resize(0);

  my_max_ped_name = 8;
  my_max_ind_name = 5;

  my_drp_ibd_info.resize(0);
  my_drp_x_ibd_info.resize(0);
}

void
lodpal_pairs::build_pairs_info()
{
  if( built_pairs_info() )
    return;

  // Clears counts
  invalidate_build_pairs_info();

  // Covariates init
  lodpal_parameters::covariate_iterator ci;
  for( ci = params.covariate_begin(); ci != params.covariate_end(); ++ci )
  {
    ci->valid = false;
    ci->variance = 0.;
    ci->variance_y = 0.;
    ci->info.clear();
    ci->info.sample_adjustment(1);
    ci->info_y.clear();
    ci->info_y.sample_adjustment(1);
    ci->info_unweighted.clear();
    ci->info_unweighted.sample_adjustment(1);
    ci->info_y_unweighted.clear();
    ci->info_y_unweighted.sample_adjustment(1);
  }

  // Weight init
  params.weight_parameter().info.clear();
  params.weight_parameter().info.sample_adjustment(1);

  // Set up trait, marker and covariate filter and create pair iterators
  make_filter();

  my_pairs_info.resize(0);

  vector<SampleInfo> weighted_minimum_covariates_info;
  weighted_minimum_covariates_info.resize(params.covariate_count());

  size_t parent_of_origin_pair = 0;

  // sib:0 hsib:1 granp:2 avunc:3 cousin:4
  vector< pair<SampleInfo, SampleInfo> > drp_ibds;
  vector< pair<SampleInfo, SampleInfo> > drp_x_ibds;
  drp_ibds.resize(5);
  drp_x_ibds.resize(22);

  // Check the number of drp pairs first
  //
  if( params.trait_parameters(0).pair_select == trait_parameter::contrast )
  {
    for( size_t pair = 0 ; pair < pairs.pair_count(); ++pair )
    {
      rel_pair p(pair, &pairs);

      assert(p.rels().pair.first->pedigree() == p.rels().pair.second->pedigree() );

      double cutpoint = params.trait_parameters(0).cutpoint;

      // Filter all invalid pairs.
      //

      // 0. Filter out unsexed pairs if sex-dependent analysis requested.
      //
      if(   (    pairs.is_x_linked(params.marker_parameters(0).marker)
              || params.autosomal_model().parent_of_origin )
          && !filter().valid_sex_info(p) )
        continue;

      // 1. Filter out invalid marker, covariate, subset, trait.
      //
      if(    !filter().valid_marker(p)
          || !filter().valid_subset(p)
          || !filter().valid_covariate(p) ) // (p, cutpoint) )
        continue;

      // 2. If not discordant pair, continue
      //
      if( !filter().is_discordant_pair(p, cutpoint) )
        continue;

      // 3. Pair type for X-linked marker
      //
      if( pairs.is_x_linked(params.marker_parameters(0).marker) )
      {
        if( p.is_mm_pair() && !params.use_mm_pair() )
          continue;

        if( p.is_ff_pair() && !params.use_ff_pair() )
          continue;

        if( p.is_mf_pair() && !params.use_mf_pair() )
          continue;
      }

      // 4. Decide the pair type & get drp info.
      //
      size_t pair_type   = get_pair_type(p);
      size_t x_pair_type = get_x_pair_type(p);    
#if 0
  cout << " pass drp pre-process ";

  cout << "** pair made to here : ";
  cout << pair << " " << pair_type << "-" << x_pair_type << " ";
  cout << p.rels().pair.first->pedigree()->name() << " ("
       << p.rels().pair.first->name()             << ","
       << p.rels().pair.second->name()            << ") = ("
       << member_trait(p.rels().pair.first, params.trait_parameters(0).trait) << ","
       << member_trait(p.rels().pair.second, params.trait_parameters(0).trait) << ")" << endl;
#endif
      if( pair_type == (size_t)-1 )
        continue;
      else if( pairs.is_x_linked(params.marker_parameters(0).marker) && x_pair_type == (size_t)-1 )
        continue;

      if( pairs.is_x_linked(params.marker_parameters(0).marker) )
      {
        drp_x_ibds[x_pair_type].first  += p.prob_share(params.marker_parameters(0).marker, 0);
        drp_x_ibds[x_pair_type].second += p.prob_share(params.marker_parameters(0).marker, 2);
      }
      else
      {
        drp_ibds[pair_type].first  += p.prob_share(params.marker_parameters(0).marker, 0);
        drp_ibds[pair_type].second += p.prob_share(params.marker_parameters(0).marker, 2);
      }
    }

    my_drp_ibd_info.resize(0);
    my_drp_x_ibd_info.resize(0);

    if( pairs.is_x_linked(params.marker_parameters(0).marker) )
    {
      for( size_t i = 0; i < drp_x_ibds.size(); ++i )
      {
        my_drp_x_ibd_info.push_back(make_pair(drp_x_ibds[i].first.mean(), drp_x_ibds[i].second.mean()));
      }
    }
    else
    {
#if 0
  cout << "drp_ibds:" << endl;
#endif
      for( size_t i = 0; i < drp_ibds.size(); ++i )
      {
        my_drp_ibd_info.push_back(make_pair(drp_ibds[i].first.mean(), drp_ibds[i].second.mean()));
#if 0
  cout << i << " " << drp_ibds[i].first.count() << ":" << drp_ibds[i].first.mean()
       << "  -  "  << drp_ibds[i].second.count() << ":" << drp_ibds[i].second.mean()
       << endl;
#endif
      }
    }
  }

  // Now, iterate over all pairs & build valid pair list
  //
  for( size_t pair = 0 ; pair < pairs.pair_count(); ++pair )
  {
    bool   arp = true;

    size_t pair_type = (size_t)-1;
    size_t x_pair_type = (size_t)-1;

    rel_pair p(pair, &pairs);

    assert(p.rels().pair.first->pedigree() == p.rels().pair.second->pedigree() );

#if 0
    cout << endl << "* pair entered : ";
    if( p.is_fsib_pair() )
      cout << "full ";
    else if( p.is_hsib_pair() )
      cout << "half ";
    else
      cout << "other ";
    cout << p.rels().pair.first->pedigree()->name() << " ("
         << p.rels().pair.first->name()             << ","
         << p.rels().pair.second->name()            << ") = ("
         << member_trait(p.rels().pair.first, params.trait_parameters(0).trait) << ","
         << member_trait(p.rels().pair.second, params.trait_parameters(0).trait) << ")" << endl;
#endif

    double cutpoint = params.trait_parameters(0).cutpoint;

    // Filter all invalid pairs.
    //

    // 0. Filter out unsexed pairs if sex-dependent analysis requested.
    //
    if(   (    pairs.is_x_linked(params.marker_parameters(0).marker)
            || params.autosomal_model().parent_of_origin )
        && !filter().valid_sex_info(p) )
      continue;

    // 1. Filter out invalid marker, covariate, subset, trait.
    //
    if(    !filter().valid_marker(p)
        || !filter().valid_subset(p)
        || !filter().valid_covariate(p) ) // (p, cutpoint) )
      continue;

#if 0
  cout << " pass 1 ";
#endif

    // 2. If concordantly affected pairs only, then filter out all invalid pairs.
    //
    bool valid_t = true;

    if( params.trait_parameters(0).pair_select == trait_parameter::conaff )
    {
#if 0
  cout << " conaff ";
#endif
      if( !filter().is_concordant_aff_pair(p, cutpoint) ) 
        valid_t = false;
    }

    // 3. If include concordant & discordant sibs(all sibs)
    //
    else if( params.trait_parameters(0).pair_select == trait_parameter::condisc )
    {
#if 0
  cout << " condisc ";
#endif
      if( !is_sib(p.rels().pair) && !filter().is_concordant_aff_pair(p, cutpoint) )
        valid_t = false;
    }

    // 5. If no concordantly unaffected sibs, then filter out concordantly unaffected
    //    sib pairs.
    //
    else if( params.trait_parameters(0).pair_select == trait_parameter::noconunaff )
    {
#if 0
  cout << " noconunaff ";
#endif
      if(    !(is_sib(p.rels().pair) || filter().is_concordant_aff_pair(p, cutpoint))
          ||  (is_sib(p.rels().pair) && filter().is_concordant_unaff_pair(p, cutpoint)) )
        valid_t = false;
    }
    else
    {
#if 0
  cout << " contrast ";
#endif
      if( !filter().is_concordant_aff_pair(p, cutpoint) )
        valid_t = false;
    }

    if( !valid_t )
      continue;

#if 0
  cout << " pass 2-3-4 ";
#endif

    // 5. Get pair weight.
    //
    double pair_weight = weight_value(p, params.weight_parameter());

    if( SAGE::isnan(pair_weight) )
      continue;

#if 0
  cout << " pass 5 ";
#endif

    // 6. Pair type for X-linked marker
    //
    if( pairs.is_x_linked(params.marker_parameters(0).marker) )
    {
      if( p.is_mm_pair() && !params.use_mm_pair() )
        continue;

      if( p.is_ff_pair() && !params.use_ff_pair() )
        continue;

      if( p.is_mf_pair() && !params.use_mf_pair() )
        continue;
    }

    pair_type   = get_pair_type(p);
    x_pair_type = get_x_pair_type(p);    

#if 0
  cout << " pass 6 ";

  cout << "** pair made to here : ";
  cout << pair_type << "-" << x_pair_type << " ";
  cout << p.rels().pair.first->pedigree()->name() << " ("
       << p.rels().pair.first->name()             << ","
       << p.rels().pair.second->name()            << ") = ("
       << member_trait(p.rels().pair.first, params.trait_parameters(0).trait) << ","
       << member_trait(p.rels().pair.second, params.trait_parameters(0).trait) << ")" << endl;
#endif

    // 7. Decide the pair type & skip if not enough drp for contrast
    //
    if( params.trait_parameters(0).pair_select == trait_parameter::contrast )
    {
      if( pairs.is_x_linked(params.marker_parameters(0).marker) )
      {
        if( x_pair_type != (size_t)-1 && !drp_x_ibds[x_pair_type].first.count() )
          continue;
      }
      else
      {
        if( pair_type != (size_t)-1 && !drp_ibds[pair_type].first.count() )
          continue;
      }
    }

    // 8. If sib pairs only, then filter out all other pairs.
    //
    if( !p.is_fsib_pair() && params.sib_pairs_only() ) 
      continue;

#if 0
  cout << " pass 8 ";
#endif

    // 9. Decide the aff status.
    if( filter().is_discordant_pair(p, cutpoint) )
    {
      arp = false;

      if( params.trait_parameters(0).pair_select == trait_parameter::contrast )
        continue;
    }

#if 0
  cout << " pass 9 ";
#endif

    params.weight_parameter().info += pair_weight;

    ++my_pair_count;

    if( p.is_fsib_pair() )
    {
      ++my_fsib_pair_count;

      double  f1   = p.prob_share(params.marker_parameters(0).marker, 1);
      double  f1mp = p.prob_share(params.marker_parameters(0).marker, 3);
      double  f1m  = (f1 + f1mp) / 2.0;
      double  f1p  = (f1 - f1mp) / 2.0;

      if( !SAGE::isnan(f1mp) && f1m != f1p )
        ++parent_of_origin_pair;
    }
    else if( p.is_hsib_pair() )
      ++my_hsib_pair_count;

#if 0
    cout << "*** good pair!!! : ";
    if( p.is_fsib_pair() )
      cout << "full No. " << my_fsib_pair_count << " ";
    else if( p.is_hsib_pair() )
      cout << "half No. " << my_hsib_pair_count << " ";
    else
      cout << "other No. " << my_pair_count << " ";
//    cout << endl;
    if( arp )  cout << "ARP ";
    else       cout << "DRP ";
    cout << pair_type << "-" << x_pair_type << " ";
    cout << p.rels().pair.first->pedigree()->name() << " ("
         << p.rels().pair.first->name()             << ","  
         << p.rels().pair.second->name()            << ") = ("
         << member_trait(p.rels().pair.first, params.trait_parameters(0).trait) << ","
         << member_trait(p.rels().pair.second, params.trait_parameters(0).trait) << ")" << endl;
#endif

    // Get pair-specific covariate value & add to SampleInfo.
    // Build my_pairs_info vector & add pair-specific covariate values to sinfo for future use.
    //
    vector<covariates_info> pair_covariates;

    for( size_t i = 0; i < params.covariate_count(); ++i )
    {
      double cv = covariate_value(p, params.covariate_parameters(i));

      params.covariate_parameters(i).info_unweighted += cv;
      params.covariate_parameters(i).info            += (cv * pair_weight);

      covariates_info ac;
      ac.pair_value       = cv;
      ac.ad_pair_value    = std::numeric_limits<double>::quiet_NaN();
      ac.re_ad_pair_value = std::numeric_limits<double>::quiet_NaN();

      pair_covariates.push_back(ac);

      if( pair_weight > 0. )
        weighted_minimum_covariates_info[i] += cv;
    }

    my_pairs_info.push_back( lodpal_pair_info(p, pair_covariates, pair_weight, arp, pair_type) );

    my_max_ped_name = std::max(p.rels().pair.first->pedigree()->name().size(), my_max_ped_name);

    my_max_ind_name = std::max(p.rels().pair.first ->name().size(), my_max_ind_name);
    my_max_ind_name = std::max(p.rels().pair.second->name().size(), my_max_ind_name);

    if( pairs.is_x_linked(params.marker_parameters(0).marker) )
    {
      lodpal_pair_info& current_pair = my_pairs_info[my_pairs_info.size()-1];

      current_pair.prior_x_ibd_index = x_pair_type;

      if( x_pair_type == (size_t)-1 )
      {
        current_pair.removed_x = true;

        if( current_pair.lodpal_pair.is_mm_pair() )
          my_mm_invalid_pair_count += 1;
        else if( current_pair.lodpal_pair.is_mf_pair() )
          my_mf_invalid_pair_count += 1;
        else if( current_pair.lodpal_pair.is_ff_pair() )
          my_ff_invalid_pair_count += 1;
      }
      else if( current_pair.lodpal_pair.is_mm_pair() )
      {
        if( params.use_mm_pair() )
        {
          my_mm_pair_count += 1;
          current_pair.removed_x = false;

          if( current_pair.lodpal_pair.is_fsib_pair() )
            my_mm_fsib_pair_count += 1;

          else if( current_pair.lodpal_pair.is_hsib_pair() )
            my_mfm_hsib_pair_count += 1;
        }
        else
          current_pair.removed_x = true;
      }
      else if( current_pair.lodpal_pair.is_mf_pair() )
      {
        if( params.use_mf_pair() )
        {
          my_mf_pair_count += 1;
          current_pair.removed_x = false;

          if( current_pair.lodpal_pair.is_fsib_pair() )
            my_mf_fsib_pair_count += 1;

          else if( current_pair.lodpal_pair.is_hsib_pair() )
            my_mff_hsib_pair_count += 1;
        }
        else
          current_pair.removed_x = true;
      }
      else if( current_pair.lodpal_pair.is_ff_pair() )
      {
        if( params.use_ff_pair() )
        {
          my_ff_pair_count += 1;
          current_pair.removed_x = false;

          if( current_pair.lodpal_pair.is_fsib_pair() )
            my_ff_fsib_pair_count += 1;

          else if( current_pair.lodpal_pair.is_hsib_pair() )
            my_fff_hsib_pair_count += 1;
        }
        else
          current_pair.removed_x = true;
      }
    }
  }

  assert( pair_count() == my_pairs_info.size() );

  if( parent_of_origin_pair >= 10 )
    my_parent_of_origin_allowed = true;

#if 0
  cout << "pairs.pair_count() = " << pairs.pair_count() << endl;
  cout << "after build, valid pair count()  = " << pair_count() << endl;
  cout << "             valid my pair count = " << my_pair_count << endl;
  cout << "             valid my fsib_pair count = " << my_fsib_pair_count << endl;
  cout << "             valid my hsib_pair count = " << my_hsib_pair_count << endl;
  cout << "             my_pairs_info.size()       = " << my_pairs_info.size() << endl;

  for(size_t i = 0; i < params.covariate_count(); ++i)
  {
    cout << params.covariate_parameters(i).name(pairs) << " unweighted - count : " 
         << params.covariate_parameters(i).info_unweighted.count() << ", mean : "
         << params.covariate_parameters(i).info_unweighted.mean() << ", sum : "
         << params.covariate_parameters(i).info_unweighted.sum() << ", variance : "
         << params.covariate_parameters(i).info_unweighted.variance() << endl;

    cout << params.covariate_parameters(i).name(pairs) << "            - count : " 
         << params.covariate_parameters(i).info.count() << ", mean : "
         << params.covariate_parameters(i).info.mean() << ", sum : "
         << params.covariate_parameters(i).info.sum() << ", variance : "
         << params.covariate_parameters(i).info.variance() << endl;
  }
  cout << "sum of weight = " << params.weight_parameter().info.sum() << endl;

  for( size_t i = 0; i < 5; ++i )
  {
    cout << "drp_ibds[" << i << "], size = " << drp_ibds[i].first.count()  << " mean = " << drp_ibds[i].first.mean() << endl;
    cout << "drp_ibds[" << i << "], size = " << drp_ibds[i].second.count() << " mean = " << drp_ibds[i].second.mean() << endl;
  }

  cout << "my_drp_ibd_info size = " << my_drp_ibd_info.size() << endl;
  for(size_t i = 0; i < my_drp_ibd_info.size(); ++i)
    cout << my_drp_ibd_info[i].first << ", " << my_drp_ibd_info[i].second << endl;

  for( size_t i = 0; i < 22; ++i )
  {
    cout << "drp_x_ibds[" << i << "], size = " << drp_x_ibds[i].first.count()  << " mean = " << drp_x_ibds[i].first.mean() << endl;
    cout << "drp_x_ibds[" << i << "], size = " << drp_x_ibds[i].second.count() << " mean = " << drp_x_ibds[i].second.mean() << endl;
  }

  cout << "my_drp_x_ibd_info size = " << my_drp_x_ibd_info.size() << endl;
  for(size_t i = 0; i < my_drp_x_ibd_info.size(); ++i)
    cout << my_drp_x_ibd_info[i].first << ", " << my_drp_x_ibd_info[i].second << endl;
#endif

  if( pair_count() == 0 )
  {
    invalidate_build_pairs_info();
    errors << priority(error) <<  "No valid pairs found!" << endl;
    return;         // Error no pairs
  }

  // Get covariate adjust method, then add ajusted covariate to SampleInfo & my_pairs_info.
  //
  for( size_t i = 0; i < params.covariate_count(); ++i )
  {
    double mu_w = params.covariate_parameters(i).info.sum() / params.weight_parameter().info.sum();
    double ad_val = mu_w;
    SampleInfo variance_x;
    SampleInfo variance_x_one_pass;

    for( size_t p = 0 ; p < my_pairs_info.size(); ++p )
    {
      double w = my_pairs_info[p].lodpal_weight;
      double x = my_pairs_info[p].lodpal_cov[i].pair_value;
      double y = 0.;

      if( params.covariate_parameters(i).adjust == covariate_type::minimum )
      {
        ad_val = weighted_minimum_covariates_info[i].min();
        y      = x - ad_val;  // minimum adjusted
      }
      else if(params.covariate_parameters(i).adjust == covariate_type::mean )
      {
        if( SAGE::isnan(params.covariate_parameters(i).adjust_value) )
        {
          // old way 1: ad_val = params.covariate_parameters(i).info.mean();
          // old way 2
          //if( w )
          //  y = w * (x - mu_w); // weighted_mean centered
          y = x - mu_w; // weighted_mean centered
        }
        else
        {
          ad_val = params.covariate_parameters(i).adjust_value;
          y      = x - ad_val; // mean centered
        }
      }
      else if(params.covariate_parameters(i).adjust == covariate_type::prop )
      {
        if( SAGE::isnan(params.covariate_parameters(i).adjust_value) )
          ad_val = 0.8;        // default proportion.
        else
          ad_val = params.covariate_parameters(i).adjust_value;

        // adjusted pair-value is still mean-centered.
        // old way : mean = params.covariate_parameters(i).info.mean()
        // old way 2
        //if( w )
        //  y = w * (x - mu_w); // weighted_mean centered
        y = x - mu_w; // weighted_mean centered
      }
      else if(params.covariate_parameters(i).adjust == covariate_type::dsp )
        y        = x;
      else
        errors << priority(warning) <<  "Unknown covariate type!!" << endl;

      my_pairs_info[p].lodpal_cov[i].ad_pair_value = y;

      params.covariate_parameters(i).info_y_unweighted += y;
      params.covariate_parameters(i).info_y            += y * w;

      variance_x += (w*((x-mu_w)*(x-mu_w)));
      variance_x_one_pass += (w*(x*x));

#if 0
  cout << p << " : x = " << setw(8) << x << ", w = " << setw(8) 
       << w << ", mu_w = " << setw(8) << mu_w << ", ad_val = "
       << setw(8) << ad_val << ", y = " << setw(8) << y << endl;
#endif

    }
    params.covariate_parameters(i).adjust_value = ad_val;
    params.covariate_parameters(i).variance = variance_x.sum()/(params.weight_parameter().info.sum()-1.0);

#if 0
  cout << endl
       << i << " : variance_x = " << variance_x.sum()/(params.weight_parameter().info.sum()-1.0) << endl;
  cout << "one-pass jane = "
       << ( variance_x_one_pass.sum()
            - ((params.covariate_parameters(i).info.sum()*params.covariate_parameters(i).info.sum())
                /params.weight_parameter().info.sum()) )
        / ( params.weight_parameter().info.sum() - 1.0 ) << endl;
  cout << "one-pass mine = "
       << ( params.covariate_parameters(i).info.sum2()
            -((params.covariate_parameters(i).info.sum()*params.covariate_parameters(i).info.sum())
               /params.weight_parameter().info.sum()) )
        / ( params.weight_parameter().info.sum() - 1.0 ) << endl;
#endif
  }

  // New addition, check design_matrix
  //
  Matrix2D<double> d_matrix(pair_count(), params.covariate_count() + 1, 1.0);

  for( size_t p = 0 ; p < my_pairs_info.size(); ++p )
    for( size_t i = 0; i < params.covariate_count(); ++i )
      d_matrix(p, i+1) = my_pairs_info[p].lodpal_cov[i].ad_pair_value;
  
  Matrix2D<double> xtx;
  xtx = XTX(d_matrix, xtx);

  Matrix2D<double> xtx_I;
  Inverse(xtx, xtx_I, 10e-10);

  double det = Det(xtx);

#if 0
  cout << endl << "d_matrix :" << endl;
  print_matrix(d_matrix, cout);

  cout << endl << "XTX :" << endl;
  print_matrix(xtx, cout);

  cout << "det = " << det << endl;

  cout << endl << "XTX_I :" << endl;
  print_matrix(xtx_I, cout);
#endif

  if( !finite(det) || !xtx_I )
  {
    errors << priority(error) << "Singular design matrix, not of full rank!!" << endl;
    return;
  }

  for(size_t i = 0; i < params.covariate_count(); ++i)
  {

#if 0
    cout << params.covariate_parameters(i).name(pairs) << " unweighted - count : " 
         << params.covariate_parameters(i).info_y_unweighted.count() << ", mean : "
         << params.covariate_parameters(i).info_y_unweighted.mean() << ", sum : "
         << params.covariate_parameters(i).info_y_unweighted.sum() << ", variance : "
         << params.covariate_parameters(i).info_y_unweighted.variance() << endl;

    cout << params.covariate_parameters(i).name(pairs) << "            - count : " 
         << params.covariate_parameters(i).info_y.count() << ", mean : "
         << params.covariate_parameters(i).info_y.mean() << ", sum : "
         << params.covariate_parameters(i).info_y.sum() << ", variance : "
         << params.covariate_parameters(i).info_y.variance() << endl;
#endif

    double mu_w = params.covariate_parameters(i).info_y.sum() / params.weight_parameter().info.sum();
    SampleInfo variance_y;
    SampleInfo variance_y_one_pass;

    for( size_t p = 0 ; p < my_pairs_info.size(); ++p )
    {
      double y = my_pairs_info[p].lodpal_cov[i].ad_pair_value;
      double w = my_pairs_info[p].lodpal_weight;

      variance_y += (w*((y-mu_w)*(y-mu_w)));
      variance_y_one_pass += (w*(y*y));
    }

    params.covariate_parameters(i).variance_y = variance_y.sum()/(params.weight_parameter().info.sum()-1.0);

#if 0
    cout << endl
         << i << " : variance_y = " << variance_y.sum()/(params.weight_parameter().info.sum()-1.0) << endl;
    cout << "one-pass jane = "
         << ( variance_y_one_pass.sum()
              - ((params.covariate_parameters(i).info_y.sum()*params.covariate_parameters(i).info_y.sum())
                  /params.weight_parameter().info.sum()) )
          / ( params.weight_parameter().info.sum() - 1.0 ) << endl;
    cout << "one-pass = "
         << ( params.covariate_parameters(i).info_y.sum2()
              -((params.covariate_parameters(i).info_y.sum()*params.covariate_parameters(i).info_y.sum())
                 /params.weight_parameter().info.sum()) )
          / ( params.weight_parameter().info.sum() - 1.0 ) << endl;
#endif
  }

  for( ci = params.covariate_begin(); ci != params.covariate_end(); ++ci )
    ci->valid = true;

  my_built_pairs_info = true;

  if( pairs.is_x_linked(params.marker_parameters(0).marker) )
  {
    build_pairs_info_x();

    if( !built_pairs_info_x() )
    {
      invalidate_build_pairs_info();
      return;
    }

    bool use_mm_pair = params.use_mm_pair();
    if( use_mm_pair && !mm_pair_count() )
      use_mm_pair = false;

    bool use_mf_pair = params.use_mf_pair();
    if( use_mf_pair && !mf_pair_count() )
      use_mf_pair = false;

    bool use_ff_pair = params.use_ff_pair();
    if( use_ff_pair && !ff_pair_count() )
      use_ff_pair = false;

    params.set_x_linkage_pair_type(use_mm_pair, use_mf_pair, use_ff_pair);
  }
}

void
lodpal_pairs::build_pairs_info_x()
{
  if( !built_pairs_info() )
    return;

  if( built_pairs_info_x() )
    return;

  my_prior_x_ibd_vector.resize(22);

  // sib
  my_prior_x_ibd_vector[0] = prior_x_ibd_type("M,M", 1./2., 1./2.,    0.);
  my_prior_x_ibd_vector[1] = prior_x_ibd_type("M,F", 1./2., 1./2.,    0.);
  my_prior_x_ibd_vector[2] = prior_x_ibd_type("F,F",    0., 1./2., 1./2.);

  // half sib
  my_prior_x_ibd_vector[3] = prior_x_ibd_type("M,F,M", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[4] = prior_x_ibd_type("M,F,F", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[5] = prior_x_ibd_type("F,F,F", 1./2., 1./2., 0.);

  // grandparental
  my_prior_x_ibd_vector[6] = prior_x_ibd_type("MFM", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[7] = prior_x_ibd_type("MFF", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[8] = prior_x_ibd_type("FFM", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[9] = prior_x_ibd_type("FFF", 1./2., 1./2., 0.);

  // avuncular
  my_prior_x_ibd_vector[10] = prior_x_ibd_type("M,FM", 3./4., 1./4., 0.);
  my_prior_x_ibd_vector[11] = prior_x_ibd_type("M,FF", 3./4., 1./4., 0.);
  my_prior_x_ibd_vector[12] = prior_x_ibd_type("F,FM", 1./4., 3./4., 0.);
  my_prior_x_ibd_vector[13] = prior_x_ibd_type("F,FF", 1./4., 3./4., 0.);
  my_prior_x_ibd_vector[14] = prior_x_ibd_type("M,MF", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[15] = prior_x_ibd_type("F,MF", 1./2., 1./2., 0.);

  // cousin
  my_prior_x_ibd_vector[16] = prior_x_ibd_type("FM,MF", 1./2., 1./2., 0.);
  my_prior_x_ibd_vector[17] = prior_x_ibd_type("MF,FM", 5./8., 3./8., 0.);
  my_prior_x_ibd_vector[18] = prior_x_ibd_type("MF,FF", 5./8., 3./8., 0.);
  my_prior_x_ibd_vector[19] = prior_x_ibd_type("FF,FF", 5./8., 3./8., 0.);
  my_prior_x_ibd_vector[20] = prior_x_ibd_type("FM,FF", 3./4., 1./4., 0.);
  my_prior_x_ibd_vector[21] = prior_x_ibd_type("FM,FM", 3./4., 1./4., 0.);

  size_t valid_pair_count = 0;
  if(    params.use_mm_pair() )
    valid_pair_count += mm_pair_count();
  if(    params.use_mf_pair() )
    valid_pair_count += mf_pair_count();
  if(    params.use_ff_pair() )
    valid_pair_count += ff_pair_count();

  if( valid_pair_count == 0 )
  {
    invalidate_build_pairs_info();
    errors << priority(error) <<  "No pairs found for X-linkage analysis!" << endl;
    return;         // Error no pairs
  }

  // New addition, check design_matrix
  //
  Matrix2D<double> d_matrix(valid_pair_count, params.covariate_count() + 1, 1.0);

  for( size_t p = 0, v = 0; p < my_pairs_info.size() && v < valid_pair_count; ++p )
  {
    if( my_pairs_info[p].removed_x )
      continue;

    for( size_t i = 0; i < params.covariate_count(); ++i )
      d_matrix(v, i+1) = my_pairs_info[p].lodpal_cov[i].ad_pair_value;
    ++v;
  }
 
  Matrix2D<double> xtx;
  xtx = XTX(d_matrix, xtx);

  Matrix2D<double> xtx_I;
  Inverse(xtx, xtx_I, 10e-10);

  double det = Det(xtx);

#if 0
  cout << endl << "d_matrix :" << endl;
  print_matrix(d_matrix, cout);

  cout << endl << "XTX :" << endl;
  print_matrix(xtx, cout);

  cout << "det = " << det << endl;

  cout << endl << "XTX_I :" << endl;
  print_matrix(xtx_I, cout);
#endif

  if( !finite(det) || !xtx_I )
  {
    errors << priority(error) << "Singular design matrix, not of full rank!!" << endl;
    return;
  }

  my_built_pairs_info_x = true;
}

void
lodpal_pairs::build_pairs_map()
{
  //assert( params.marker_count() == 1);
  if( !my_built_pairs_info )
    return;

  my_pairs_map.clear();

  // Iterate over all pairs & build valid pair list
  for( size_t pair = 0 ; pair < my_pairs_info.size(); ++pair )
  {
    rel_pair p = my_pairs_info[pair].lodpal_pair;

    assert(p.rels().pair.first->pedigree() == p.rels().pair.second->pedigree() );

//    my_rel_pair_map[p] = 0.0001;
    my_pairs_map[p] = 0.45;
  }

  assert( my_pairs_info.size() == my_pairs_map.size() );
}

void
lodpal_pairs::re_build_pairs_info()
{
  assert( !re_built_pairs_info() );

  bool   bi_cov_exist = false;
  bool   re_cov_exist = false;
  size_t bi_cov = 0;
  size_t re_cov = 0;

  for( size_t i = 0; i < params.covariate_count(); ++i )
  {
    const RPED::RefTraitInfo& info = pairs.fped_info().trait_info(params.covariate_parameters(i).covariate);
    if( info.type() == RPED::RefTraitInfo::binary_trait )
    {
      bi_cov_exist = true;
      bi_cov = i;
      break;
    }
  }

  if(    bi_cov_exist
      && (    params.covariate_parameters(bi_cov).adjust == covariate_type::minimum
           || params.covariate_parameters(bi_cov).adjust == covariate_type::dsp ) )
  {
    re_cov = bi_cov;
    re_cov_exist = true;
  }
  else
  {
    for( size_t i = 0; i < params.covariate_count(); ++i )
    {
      if(    params.covariate_parameters(i).adjust == covariate_type::minimum
          || params.covariate_parameters(i).adjust == covariate_type::dsp )
      {
        re_cov = i;
        re_cov_exist = true;
        break;
      }
    }
  }

  assert(re_cov_exist);

  double max_y = params.covariate_parameters(re_cov).info_y.max();
  for( size_t p = 0 ; p < my_pairs_info.size(); ++p )
  {
    double y     = my_pairs_info[p].lodpal_cov[re_cov].ad_pair_value;
    double y_    = max_y - y;
    my_pairs_info[p].lodpal_cov[re_cov].re_ad_pair_value = y_;

//cout << "y = " << y << ", y* = " << y_ << endl;
  }

  my_re_built_covariate  = re_cov;
  my_re_built_pairs_info = true;
}

size_t
lodpal_pairs::remove_biggest_pair(bool remove_x)
{
  // Remove a pair with biggest lod_score.
  double max_lod = -std::numeric_limits<double>::infinity();

  size_t removed_pos = pairs_info().size();

  for( size_t i = 0; i < pairs_info().size(); ++i )
  {
    if( pairs_info()[i].removed )
      continue;

    if( remove_x && pairs_info()[i].removed_x )
      continue;

    if( pairs_info()[i].likelihood > max_lod )
    {
      max_lod      = pairs_info()[i].likelihood;
      removed_pos  = i;
    }
  }

  if( finite(max_lod) )
    pairs_info()[removed_pos].removed = true;

//cout << "removed max_lod = " << max_lod << ", pos = " << removed_pos <<endl;
  rel_pair p    = pairs_info()[removed_pos].lodpal_pair;

  if( p.is_fsib_pair() )
    my_removed_fsib_pairs.push_back(removed_pos);
  else if( p.is_hsib_pair() )
    my_removed_hsib_pairs.push_back(removed_pos);
  else
    my_removed_other_pairs.push_back(removed_pos);

  assert( removed_pos < pairs_info().size() );

  return removed_pos;
}

void
lodpal_pairs::reset_removed_pairs()
{
  for( size_t i = 0; i < my_removed_fsib_pairs.size(); ++i )
    pairs_info()[my_removed_fsib_pairs[i]].removed = false;
  
  for( size_t i = 0; i < my_removed_hsib_pairs.size(); ++i )
    pairs_info()[my_removed_hsib_pairs[i]].removed = false;

  for( size_t i = 0; i < my_removed_other_pairs.size(); ++i )
    pairs_info()[my_removed_other_pairs[i]].removed = false;

  my_removed_fsib_pairs.resize(0);
  my_removed_hsib_pairs.resize(0);
  my_removed_other_pairs.resize(0);
}

void
lodpal_pairs::make_filter()
{
  my_filter = pair_filter();

  pair_filter& pf = my_filter;

  //size_t p = params.parameter_count();
  size_t m = params.marker_count();
  size_t c = params.covariate_count();
  size_t t = params.trait_count();
  size_t s = params.subset_count();

  double tolerance = std::numeric_limits<double>::infinity();

  if( params.skip_uninformative_pairs() )
    tolerance = 0.001;

  for( size_t i = 0; i < m; ++i )
  {
    size_t mm = params.marker_parameters(i).marker;
    if( mm < pairs.marker_count() )
      pf.add_marker(mm, tolerance);
  }

  for( size_t i = 0; i < c; ++i )
  {
    size_t cc = params.covariate_parameters(i).covariate;
#if 0
  cout << "ar filter, cc = " << cc << ", pairs.trait_count() = " << pairs.trait_count() << endl;
  for( size_t j = 0; j < pairs.trait_count(); ++j )
    cout << j << " : " << pairs.trait_name(j) << endl;
  cout << endl; 
#endif

    if(    params.covariate_parameters(i).operation != covariate_type::pair
        && cc < pairs.trait_count() )
    {
//      cout << "cov name = " << pairs.trait_name(cc) << endl;
      pf.add_covariate(cc, -std::numeric_limits<double>::infinity());
    }
    else if(    params.covariate_parameters(i).operation == covariate_type::pair
        && cc < pairs.pair_covariate_count() )
    {
//      cout << "pair cov name = " << pairs.get_pair_covariate_name(cc) << endl;
      pf.add_pair_covariate(cc, -std::numeric_limits<double>::infinity());
    }
  }

  for( size_t j = 0; j < t; ++j )
  {
    size_t tt = params.trait_parameters(j).trait;
    if( tt < pairs.trait_count() )
      pf.add_trait(tt, (size_t)2, -std::numeric_limits<double>::infinity());
  }

  for(size_t j = 0; j < s; ++j)
  {
    size_t tt = params.subset_parameters(j).trait;
    if(tt < pairs.trait_count() )
      pf.add_subset(tt, (size_t)2, -std::numeric_limits<double>::infinity());
  }
}

double
lodpal_pairs::covariate_value(const rel_pair& pair, 
                              const covariate_type& param) const
{
  cout << flush;

  size_t  cc = param.covariate;
  double  cv = numeric_limits<double>::quiet_NaN();

  const RelativePairs& pairs = *pair.pair_data();

  if( param.operation == covariate_type::pair )
  {
    if( cc >= pairs.pair_covariate_count() )
      return cv;

    cv = pairs.get_pair_covariate(pair.pair_number(), cc);

    return cv;
  }

  if( cc >= pairs.trait_count() )
    return cv;

  double cv1 = member_trait(pair.rels().pair.first,  cc);
  double cv2 = member_trait(pair.rels().pair.second, cc);

  const double cutpoint = params.trait_parameters(0).cutpoint;
  switch( param.operation )
  {
    case covariate_type::sum:    cv = cv1 + cv2;            break;
    case covariate_type::none:   cv1 = (cv1 > cutpoint) ? 1 : 0;
                                 cv2 = (cv2 > cutpoint) ? 1 : 0;
                                 cv = (cv1 == cv2) ? 0 : 1; break;
    case covariate_type::diff:   cv = fabs(cv1 - cv2);      break;
    case covariate_type::prod:   cv = cv1 * cv2;            break;
    case covariate_type::single: cv = cv1;                  break;
    case covariate_type::avg:    cv = (cv1 + cv2)/2.;       break;
    default: break;
  }

  if( param.power != 1.0 )
    cv = pow(cv, param.power);

#if 0
  cout << "c_name = " << param.name(pairs) << endl;
  cout << "cc  = " << cc  << endl;
  cout << "cv1 = " << cv1 << endl;
  cout << "cv2 = " << cv2 << endl;
  cout << "cv  = " << cv  << endl;
#endif

  return cv;
}

double
lodpal_pairs::weight_value(const rel_pair& pair, 
                           const lodpal_weight_type& param) const
{
  size_t  ww = param.weight;
  double  wv = 1.0;

  const RelativePairs& pairs = *pair.pair_data();

  if( param.operation == lodpal_weight_type::pair )
  {
    if( ww >= pairs.pair_weight_count() )
      return wv;

    wv = pairs.get_pair_weight(pair.pair_number(), ww);

    if( 0 > wv || 1 < wv || !finite(wv) )
      return numeric_limits<double>::quiet_NaN();

    return wv;
  }

  if( ww >= pairs.trait_count() )
    return wv;

  wv = member_trait(pair.rels().pair.first,  ww);

#if 0
  cout << "w_name = " << param.name(pairs) << endl;
  cout << "ww  = " << ww  << endl;
  cout << "wv  = " << wv  << endl;
#endif

  if( 0 > wv || 1 < wv || !finite(wv) )
    return numeric_limits<double>::quiet_NaN();

  return wv;
}

size_t
lodpal_pairs::get_pair_type(const rel_pair& p) const
{
  pair_type r_t = p.rels().type;

  if( r_t == pair_generator::SIBSIB )
    return 0;
  else if( r_t == pair_generator::HALFSIB )
    return 1;
  else if( r_t == pair_generator::GRANDP )
    return 2;
  else if( r_t == pair_generator::AVUNC )
    return 3;
  else if( r_t == pair_generator::COUSIN )
    return 4;

  return (size_t)-1;
}

size_t
lodpal_pairs::get_x_pair_type(const rel_pair& p) const
{
  pair_x_type r_t = p.rels().x_type;

  if( r_t == INVALID )
    return (size_t)-1;

  return (size_t)r_t;
}

} // end of namespace LODPAL
} // end of namespace SAGE
