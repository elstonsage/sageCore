//****************************************************************************
//* File:      lodpal_analysis.cpp                                           *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Initial implementation                                        *
//*                                                                          *
//* Notes:     This file implements lodpal_analysis class.                   *
//*                                                                          *
//* Copyright (c) 2006 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_analysis.h"

namespace SAGE   {
namespace LODPAL {

////////////////////////////////////////////////////////////////////////////
//          Implementation of lodpal_analysis(non-Inline)                 //
////////////////////////////////////////////////////////////////////////////

lodpal_analysis::lodpal_analysis(cerrorstream& err)
           : errors(err)
{
  invalidate();
}

void
lodpal_analysis::invalidate()
{
  my_built           = false;
  my_method_D_called = 0;
}

void
lodpal_analysis::do_lodpal_test(RelativePairs &pairs, lodpal_parser& parser,
                                ostream &out, ostream &dig, ostream &xln)
{
  lodpal_parameters params = parser.parameters();

  bool wide_output  = parser.wide_output();
  bool csv_output   = parser.csv_output();
  bool pval_sci     = parser.get_pval_scientific_notation();

  // Set the locations to analyse.
  //
  if( !params.marker_count() )
  {
    //cout << "No markers listed -- using all markers" << endl;
    for(size_t m=0; m < pairs.marker_count(); ++m)
    {
      marker_type::inheritance_type in = marker_type::autosomal;

      if( pairs.is_x_linked(m) )
      {
        in = marker_type::x_linked;
        params.set_x_linked_marker_exist(true);
      }
      else
        params.set_autosomal_marker_exist(true);

      params.add_marker(m, in);
    }
  }

  // Get the maximum marker name length.
  //
  size_t max_marker = 10;

  for(size_t i = 0; i < params.marker_count(); ++i)
  {
    string marker_name = params.marker_parameters(i).name(pairs);

    max_marker = std::max(marker_name.size(), max_marker);
  }

  params.set_max_marker_name_size(max_marker);

  // Set the trait to analyse.
  //
  if( !params.trait_count() )
  {
    //out << "Using all traits since none were specified." << endl;

    for(size_t t=0; t < pairs.trait_count(); ++t)
    {
      const RPED::RefTraitInfo &info = pairs.fped_info().trait_info(t);

      if( info.usage() != RPED::RefTraitInfo::trait_variate)
        continue;

      else if(info.type() == RPED::RefTraitInfo::binary_trait)
        params.add_trait(t, trait_parameter::conaff);
    }
  }

  // Set the weight if any exist.
  //
  if(    params.weight_parameter().weight == (size_t)-1 
      && pairs.pair_weight_count() == 1 )
    params.add_weight(0, lodpal_weight_type::pair);

  // Initialize the output result file(s).
  //
  typedef boost::shared_ptr<lodpal_test_viewer> ro_ptr;
  ro_ptr test_output;
  ro_ptr diag_output;
  ro_ptr xlin_output;

  if( csv_output )
    test_output = ro_ptr( new lodpal_test_csvfile(out, wide_output, pval_sci)  );
  else
  {
    if( params.autosomal_marker_exist() )
      test_output = ro_ptr( new lodpal_test_textfile(out, wide_output, pval_sci) );

    if( params.x_linked_marker_exist() )
      xlin_output = ro_ptr( new lodpal_test_xfile(xln, wide_output, pval_sci) );
  }

  if( parser.diagnostic() )
    diag_output = ro_ptr( new lodpal_test_diagfile(dig) );

  // Trait by trait
  //
  for(size_t t = 0; t < params.trait_count(); ++t)
  {
    lodpal_parameters local_params = params;

    local_params.clear_parameters();
    local_params.set_trait( params.trait_parameters(t) );

    // Add covariates.
    //
    if(    (   params.trait_parameters(t).pair_select == trait_parameter::conaff
            || params.trait_parameters(t).pair_select == trait_parameter::contrast )
        && params.autosomal_model().model == autosomal_model_type::one_parameter
        && params.x_linkage_model().lambda2_fixed )
    {
      for( size_t i = 0; i < params.covariate_count(); ++i )
        local_params.add_parameter( params.covariate_parameters(i) );

      if( pairs.pair_covariate_count() )
        for( size_t i = 0; i < pairs.pair_covariate_count(); ++i )
          if( pairs.pair_covariate_info(i).get_usage() == PALBASE::pair_pheno_info::minimum )
            local_params.add_covariate(i, covariate_type::minimum, covariate_type::pair,
                                       pairs.pair_covariate_info(i).get_value());
          else
            local_params.add_covariate(i, covariate_type::mean, covariate_type::pair,
                                       pairs.pair_covariate_info(i).get_value());
    }

    // If include discordant sib pairs, user-defined covariates ignored.
    // A new covariate-concordance status created.
    //
    else if(    params.autosomal_model().model == autosomal_model_type::one_parameter
             && params.x_linkage_model().lambda2_fixed )
     local_params.add_covariate(params.trait_parameters(t).trait,
                                covariate_type::dsp, covariate_type::none);

    // If covariate exists, lambda can't be printed.
    //
    if( local_params.covariate_count() && local_params.print_lambda() )
    {
      errors << priority(warning)
             << "The parameter estimates can't be print as lambda values when covariate(s) are used. "
             << "'lambda' option will be ignored." << endl;

      local_params.set_print_lambda(false);
    }

    // Create the copies of param for later use.
    //
    lodpal_parameters previous_local_params = local_params;
    previous_local_params.clear_markers();

    lodpal_parameters original_local_params = previous_local_params;

    // Start the analysis, one marker by one marker
    //
    for(size_t i = 0; i < params.marker_count(); ++i)
    {
      bool valid_data = true;

      // Update the marker parameter
      //
      local_params.clear_markers();
      local_params.add_parameter(params.marker_parameters(i));

      // Create pairs_info.
      //
      lodpal_pairs  l_pairs(pairs, local_params, errors);

      l_pairs.build_pairs_info();
      if( !l_pairs.built_pairs_info() )
      {
        errors << priority(error) 
               << "Test not successful at '"
               << pairs.marker_name(local_params.marker_parameters(0).marker)
               << "'!" << endl;

        valid_data = false;
      }

      // Build relative pair map to store each pair likelihood.
      //
      l_pairs.build_pairs_map();

      // Create analysis object.
      //

      if( pairs.is_x_linked(local_params.marker_parameters(0).marker) )
      {
        if( local_params.x_linkage_model().lambda2_fixed )
          arp = analysis_ptr( new ARP_x_one_analysis(l_pairs, errors) );
        else
        {
          if(    (    local_params.x_linkage_model().lambda1_equal
                   && local_params.use_mm_pair() && params.use_mf_pair()
                   && l_pairs.mm_fsib_pair_count() + l_pairs.mf_fsib_pair_count() >= 15 )
              || (   !local_params.x_linkage_model().lambda1_equal
                   && local_params.use_ff_pair()
                   && (    l_pairs.ff_fsib_pair_count() >= 15
                        && l_pairs.ff_pair_count() - l_pairs.ff_fsib_pair_count() >= 15 ) ) )
            arp = analysis_ptr( new ARP_x_two_analysis(l_pairs, errors)  );
          else
          {
            errors << priority(warning) 
                   << "The parameter lambda2 can not be estimated with this data set at the location '"
                   << pairs.marker_name(local_params.marker_parameters(0).marker)
                   << "'.  The lambda2 is fixed." << endl;
            local_params.x_linkage_model().lambda2_fixed = true;
            arp = analysis_ptr( new ARP_x_one_analysis(l_pairs, errors)  );
          }
        }
      }
      else
      {
        if( local_params.autosomal_model().parent_of_origin )
        {
          if( l_pairs.parent_of_origin_allowed() )
          {
            if( local_params.autosomal_model().model == autosomal_model_type::one_parameter )
              arp = analysis_ptr( new ARP_po_one_analysis(l_pairs, errors) );
            else
              arp = analysis_ptr( new ARP_po_two_analysis(l_pairs, errors)  );
          }
          else
          {
            errors << priority(warning) 
                   << "The data set has less than 10 ASPs in which maternal and paternal "
                   << "ibd sharing are not equal at the location '"
                   << pairs.marker_name(local_params.marker_parameters(0).marker)
                   << "'.  The parent-of-origin test is not allowed.  Skipping..." << endl;

            valid_data = false;

            if( local_params.autosomal_model().model == autosomal_model_type::one_parameter )
              arp = analysis_ptr( new ARP_po_one_analysis(l_pairs, errors) );
            else
              arp = analysis_ptr( new ARP_po_two_analysis(l_pairs, errors)  );
          }
        }
        else
        {
          if( local_params.autosomal_model().model == autosomal_model_type::one_parameter )
            if(    params.trait_parameters(t).pair_select == trait_parameter::condisc
                || params.trait_parameters(t).pair_select == trait_parameter::noconunaff )
              arp = analysis_ptr( new DSP_one_analysis(l_pairs, errors) );
            else   
              arp = analysis_ptr( new ARP_one_analysis(l_pairs, errors) );
          else
            arp = analysis_ptr( new ARP_two_analysis(l_pairs, errors) );
        }
      }

      arp->set_parameters(local_params);

      if( valid_data )
      {
        // If multipoint, use parameter estimates from the previous point
        //
        if( arp->parameters().multipoint() )
          arp->set_previous_parameters(previous_local_params);

        bool good_bound = try_run();

        if( !good_bound )
        {
          arp->pairs_info().re_build_pairs_info();
          run();
        }

        if( arp->pairs_info().re_built_pairs_info() )
          arp->re_estimate_parameters();

        if( !arp->parameters().valid() )
        {
          errors << priority(error) 
                 << "Test not successful due to invalid parameters at '"
                 << pairs.marker_name(local_params.marker_parameters(0).marker)
                 << "'!" << endl;
        }

        // Set the previous parameter estimates.
        //
        //if( arp->parameters().multipoint() && arp->built() )
        if( arp->parameters().multipoint() && my_built )
          previous_local_params = arp->parameters();
      }

      // Print out the header.
      //
      if( i == 0 )
      {
        if( params.autosomal_marker_exist() )
          test_output->print_header(*arp);

        if( params.x_linked_marker_exist() )
          xlin_output->print_header(*arp);
      }

      // Print out the results.
      //
      if( valid_data )
      {
        if(    params.x_linked_marker_exist()
            && arp->relative_pairs().is_x_linked(arp->parameters().marker_parameters(0).marker) )
          xlin_output->print_results(*arp);
        else
          test_output->print_results(*arp);
      }
      else
      {
        if(    params.x_linked_marker_exist()
            && arp->relative_pairs().is_x_linked(arp->parameters().marker_parameters(0).marker) )
          xlin_output->print_skip_line(*arp);
        else
          test_output->print_skip_line(*arp);
      }

      if(    parser.diagnostic()
          && arp->parameters().marker_parameters(0).marker == arp->parameters().diagnostic_marker() )
      {
        diag_output->print_header(*arp);
        diag_output->print_results(*arp);
        diag_output->print_footer(*arp);
      }

      // Print out the footer
      //
      if( i == params.marker_count() - 1 )
      {
        if( params.autosomal_marker_exist() )
          test_output->print_footer(*arp);

        if( params.x_linked_marker_exist() )
          xlin_output->print_footer(*arp);
      }
    }
  }
}

void
lodpal_analysis::run()
{
  try_run();
}

bool
lodpal_analysis::try_run()
{
  if( !arp->built() )
    arp->build();

  if( !arp->built() )
    return true;

  if( !arp->valid_parameter_count() )
  {
    cout << "No valid parameters?!" << endl;
    return true;
  }

  assert( arp->parameters().marker_count() == 1 );

  arp->set_good_bound(true);

  run_model();

  return  arp->is_good_bound();
}

void
lodpal_analysis::run_model()
{
  invalidate();

  // Temporary place to store the results.            
  vector<inter_results>                temp_results;
  vector<lodpal_pairs::pair_info_type> temp_pair_info;
  vector<lodpal_pairs::pair_map_type>  temp_pair_map;
  vector<vector<size_t> >              temp_removed_fsib_pairs;
  vector<vector<size_t> >              temp_removed_hsib_pairs;
  vector<vector<size_t> >              temp_removed_other_pairs;
  vector<std::pair<bool, size_t> >     temp_built;

  temp_results.resize(0);
  temp_pair_info.resize(0);
  temp_pair_map.resize(0);
  temp_removed_fsib_pairs.resize(0);
  temp_removed_hsib_pairs.resize(0);
  temp_removed_other_pairs.resize(0);
  temp_built.resize(0);

  if( arp->previous_parameters().marker_count() )
  {
    temp_results.resize(4);
    temp_pair_info.resize(4);
    temp_pair_map.resize(4);
    temp_removed_fsib_pairs.resize(4);
    temp_removed_hsib_pairs.resize(4);
    temp_removed_other_pairs.resize(4);
    temp_built.resize(4);
  }
  else
  {
    temp_results.resize(3);
    temp_pair_info.resize(3);
    temp_pair_map.resize(3);
    temp_removed_fsib_pairs.resize(3);
    temp_removed_hsib_pairs.resize(3);
    temp_removed_other_pairs.resize(3);
    temp_built.resize(3);
  }

  // Reserve pairs_info.
  lodpal_pairs reserved_pairs = arp->pairs_info();

  vector< maxfun_param_mgr >   init_theta;
  encode_maxfun_params(init_theta);

  inter_results result;

  for( size_t ti = 0; ti < init_theta.size(); ++ti )
  {
    my_built = true;
    my_method_D_called = 0;

    inter_results resultA = run_method_B(init_theta[ti]);
    //inter_results resultA = run_method_A(init_theta[ti]);
    lodpal_pairs reserved_pairsA = arp->pairs_info();

    //inter_results resultB = run_method_B(init_theta[ti]);
    //lodpal_pairs reserved_pairsB = arp->pairs_info();

    double diffA = fabs(resultA.lod_score_with_cap - resultA.lod_score_without_cap);
    //double diffB = fabs(resultB.lod_score_with_cap - resultB.lod_score_without_cap);

#if 0
    cout << "diffA = " << resultA.lod_score_with_cap << "-" << resultA.lod_score_without_cap
         << " = " << fabs(resultA.lod_score_with_cap - resultA.lod_score_without_cap) << endl;
    //cout << "diffB = " << resultB.lod_score_with_cap << "-" << resultB.lod_score_without_cap
    //     << " = " << fabs(resultB.lod_score_with_cap - resultB.lod_score_without_cap) << endl;
#endif

    if( diffA > 0.75 )
    {
      if( !arp->parameters().turn_off_default() )
      {
        result = run_method_D(init_theta[ti]);

        // Find the best converged run to report instead of -----
        //
        if( !my_built && isnan(result.lod_score_without_cap) )
        {
          if( resultA.maxfun_info.getExitFlag() < 4 )
          {
            arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
            arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

            result = resultA;
            my_built = true;
          }
        }
        //
      }
      else
      {
        arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
        arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

        result = resultA;
      }      
    }
    else
    {
      arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
      arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

      result = resultA;
    }
/*
    if( diffA > 0.75 && diffB > 0.75 )
    {
      if( !arp->parameters().turn_off_default() )
      {
        result = run_method_D(init_theta[ti]);

        // Find the best converged run to report instead of -----
        //
        if( !my_built && isnan(result.lod_score_without_cap) )
        {
          if( resultA.maxfun_info.getExitFlag() < 4 && resultB.maxfun_info.getExitFlag() < 4 )
          {
            if( resultA.lod_score_without_cap >= resultB.lod_score_without_cap )
            {
              arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
              arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

              result = resultA;
              my_built = true;
            }
            else
            {
              arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
              arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();

              result = resultB;
              my_built = true;
            }
          } 
          else if( resultA.maxfun_info.getExitFlag() < 4 )
          {
            arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
            arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

            result = resultA;
            my_built = true;
          }
          else if( resultB.maxfun_info.getExitFlag() < 4 )
          {
            arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
            arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();

            result = resultB;
            my_built = true;
          }
        }
        //
      }
      else
      {
        if( resultA.lod_score_without_cap >= resultB.lod_score_without_cap )
        {
          arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
          arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

          result = resultA;
        }
        else
        {
          arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
          arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();

          result = resultB;
        }
      }      
    }
    else if( diffA <= 0.75 && diffB <= 0.75 )
    {
      if( fabs(resultA.lod_score_without_cap - resultB.lod_score_without_cap) <= 0.05 )
      {
        arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
        arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

        result = resultA;
      }
      else
      {
        if( resultA.lod_score_without_cap > resultB.lod_score_without_cap )
        {
          arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
          arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

          result = resultA;
        }
        else
        {
          arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
          arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();

          result = resultB;
        }
      }
    }
    else if( diffA <= 0.75 && diffB > 0.75 )
    {
      arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
      arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();

      result = resultA;
    }
    else
    {
      arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
      arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();

      result = resultB;
    }
*/
#if 0
    cout << "--result for init_theta[" << ti << "]" << endl;
    cout << "  lod with cap    = " << result.lod_score_with_cap << endl
         << "  lod without cap = " << result.lod_score_without_cap << endl;
    cout << "  param count     = " << result.maxfun_info.getParameterMgr().getParamCount() << endl;
    for( int p = 0; p < result.maxfun_info.getParameterMgr().getParamCount(); ++p )
    {
      cout << "  " << p 
           << ": " << result.maxfun_info.getParameterMgr().getParameter(p).getFinalEstimate();
      if( result.maxfun_info.getParameterMgr().getParameter(p).isDerivAvailable() )
        cout << " (" << result.maxfun_info.getParameterMgr().getParameter(p).getDeriv() << ") ";
      cout << endl;
    }
    if( my_built )
      cout << "  built true, pair-count = " << arp->pairs_info().pair_count() << endl;
    else
      cout << "  built false, pair-count = " << arp->pairs_info().pair_count() << endl;
#endif
 
    temp_results[ti]   = result;
    temp_pair_info[ti] = arp->pairs_info().pairs_info();
    temp_pair_map[ti]  = arp->pairs_info().pairs_map();
    temp_built[ti]     = make_pair(my_built, arp->pairs_info().pair_count());
    temp_removed_fsib_pairs[ti]  = arp->pairs_info().removed_fsib_pairs();
    temp_removed_hsib_pairs[ti]  = arp->pairs_info().removed_hsib_pairs();
    temp_removed_other_pairs[ti] = arp->pairs_info().removed_other_pairs();

    arp->pairs_info().pairs_info() = reserved_pairs.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairs.pairs_map();
    arp->pairs_info().reset_removed_pairs();
  }

  double max_lod = -std::numeric_limits<double>::infinity();
  size_t max_pos = temp_results.size();
  size_t max_pairs = reserved_pairs.pair_count();

  // Find the max_lod_score with least pairs removed.
  for( size_t ti = 0; ti < temp_results.size(); ++ti )
  {
    if( !temp_built[ti].first )
      continue;

    double temp_lod = temp_results[ti].lod_score_without_cap;

    if(    temp_lod > max_lod
        && temp_built[ti].second == max_pairs )
    {
      max_lod = temp_lod;
      max_pos = ti;
    }
  }

#if 0
  cout << "max_pos1 = " << max_pos << endl;
#endif

  if( max_pos == temp_results.size() )
  {
    for( size_t ti = 0; ti < temp_results.size(); ++ti )
    {
      if( !temp_built[ti].first )
        continue;

      double temp_lod = temp_results[ti].lod_score_without_cap;

      if( temp_lod > max_lod )
      {
        max_lod = temp_lod;
        max_pos = ti;
      }
    }
  }

#if 0
  cout << "max_pos2 = " << max_pos << endl;
#endif

  if(    max_pos < temp_results.size()
      && temp_results[max_pos].lod_score_without_cap >= 0. )
  {
    my_built = true;

    decode_final_result(temp_results[max_pos]);
    arp->pairs_info().pairs_info() = temp_pair_info[max_pos];
    arp->pairs_info().pairs_map()  = temp_pair_map[max_pos];

    if( temp_built[max_pos].second == max_pairs )
      arp->pairs_info().reset_removed_pairs();
    else
    {
      arp->pairs_info().removed_fsib_pairs() = temp_removed_fsib_pairs[max_pos];
      arp->pairs_info().removed_hsib_pairs() = temp_removed_hsib_pairs[max_pos];
      arp->pairs_info().removed_other_pairs() = temp_removed_other_pairs[max_pos];
    }
  }
  else
    my_built = false;

  return;
}

lodpal_analysis::inter_results
lodpal_analysis::run_maxfun(maxfun_seq_cfg& sc, maxfun_param_mgr& init_theta)
{
  maxfun_debug   db;

  inter_results result;

  result.maxfun_info = maxfun_maximizer::Maximize(sc, init_theta, *arp, db);

  double max_lod = result.maxfun_info.getFinalFunctionValue();

#if 0
  cout << "max_lod with cap = " << max_lod << endl;
  cout << "        lfl      = " << result.maxfun_info.getExitFlag() << endl
       << "        nfe      = " << arp->nfe << endl
       << "        ivfl     = " << result.maxfun_info.getCovMatrixStatus() << endl;
#endif

  if( fabs(max_lod) < 1.e-10 )
  {
    max_lod = 0.0;
    result.maxfun_info.setFinalFunctionValue(max_lod);
  }

  result.lod_score_with_cap = max_lod;

  decode_theta(result.maxfun_info, result.param_estimates);

  max_lod = arp->re_evaluate(result.param_estimates);

  if( fabs(max_lod) < 1.e-10 )
    max_lod = 0.0;

  result.lod_score_without_cap = max_lod;

#if 0
  cout << "max_lod without cap = " << max_lod << endl;
#endif

  if( arp->pairs_info().relative_pairs().is_x_linked(arp->current_marker().marker) )
  {
    if( arp->get_lod_mm() < 0. && arp->get_lod_mm() > -1e-10 )
      result.lod_score_mm = 0.0;
    else
      result.lod_score_mm = arp->get_lod_mm();

    if( arp->get_lod_mf() < 0. && arp->get_lod_mf() > -1e-10 )
      result.lod_score_mf = 0.0;
    else
      result.lod_score_mf = arp->get_lod_mf();

    if( arp->get_lod_ff() < 0. && arp->get_lod_ff() > -1e-10 )
      result.lod_score_ff = 0.0;
    else
      result.lod_score_ff = arp->get_lod_ff();
  }

  return result;
}

lodpal_analysis::inter_results
lodpal_analysis::run_method_A(maxfun_param_mgr& init_theta)
{
  //cout << "run_method_A..." << endl;

  // Geoff's standard optimum methods
  //int    max_mthds[] = { 2,     5,     2,     6 };
  //double max_epsc1[] = { 1e-3,  1e-3,  1e-4,  1e-12 };
  //double max_epsc2[] = { 1e-15, 1e-15, 1e-15, 1e-15 };
  //int    max_maxit[] = { 1,     20,    50,    20 };

  // Dr. Elston's new suggestion 09/20/2010
  //int    max_mthds[] = { 2,    5,    2,    6,    6 };
  //double max_epsc1[] = { 1e-2, 1e-3, 1e-4  1e-5, 1e-6 };
  //double max_epsc2[] = { 1e-8, 1e-8, 1e-8, 1e-9, 1e-9 };
  //int    max_maxit[] = { 1,    20,   50,   20,   20 };

  //maxfun_debug   db;
  maxfun_seq_cfg sc(maxfun_seq_cfg::USER_DEFINED, "A");

  sc.addRunCfg(maxfun_run_cfg::DIRECT_WITHOUT, 1);
  sc.getLatestRunCfg().epsilon1 = 1e-2;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_IDENTITY, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-3;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::DIRECT_WITHOUT, 50);
  sc.getLatestRunCfg().epsilon1 = 1e-4;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_ESTIMATE, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-5;
  sc.getLatestRunCfg().epsilon2 = 1e-9;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_ESTIMATE, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-6;
  sc.getLatestRunCfg().epsilon2 = 1e-9;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  inter_results resultA = run_maxfun(sc, init_theta);

  resultA.method = 'A';

  return resultA;
}

lodpal_analysis::inter_results
lodpal_analysis::run_method_B(maxfun_param_mgr& init_theta)
{
  //cout << "run_method_B..." << endl;

  //int    max_mthds[] = { 1 };
  //double max_epsc1[] = { 1e-3 };
  //double max_epsc2[] = { 1e-15 };
  //int    max_maxit[] = { 10 };

  maxfun_seq_cfg sc(maxfun_seq_cfg::USER_DEFINED, "B");

  sc.addRunCfg(maxfun_run_cfg::COMPLETE_DIRECT, 10);
  sc.getLatestRunCfg().epsilon1 = 1e-3;
  sc.getLatestRunCfg().epsilon2 = 1e-15;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;
/*
  sc.addRunCfg(maxfun_run_cfg::DIRECT_WITHOUT, 1);
  sc.getLatestRunCfg().epsilon1 = 1e-2;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_IDENTITY, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-3;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::DIRECT_WITHOUT, 50);
  sc.getLatestRunCfg().epsilon1 = 1e-4;
  sc.getLatestRunCfg().epsilon2 = 1e-8;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_ESTIMATE, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-5;
  sc.getLatestRunCfg().epsilon2 = 1e-9;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;

  sc.addRunCfg(maxfun_run_cfg::VAR_METRIC_ESTIMATE, 20);
  sc.getLatestRunCfg().epsilon1 = 1e-6;
  sc.getLatestRunCfg().epsilon2 = 1e-9;
  sc.getLatestRunCfg().var_cov = maxfun_run_cfg::FINAL;
*/
  inter_results resultB = run_maxfun(sc, init_theta);

  resultB.method = 'B';

  return resultB;
}

lodpal_analysis::inter_results
lodpal_analysis::run_method_D(maxfun_param_mgr& init_theta)
{
  my_method_D_called += 1;

  //cout << "run_method_D..." << endl;

  // Remove a pair with biggest lod_score.
  //cout << "removed 1 pair at method_D" << endl;
  size_t removed_pos = arp->pairs_info().remove_biggest_pair(true);

  lodpal_pairs reserved_pairs = arp->pairs_info();  // original my_pair - 1

  inter_results resultA = run_method_B(init_theta);
  //inter_results resultA = run_method_A(init_theta);
  lodpal_pairs reserved_pairsA = arp->pairs_info();

  //inter_results resultB = run_method_B(init_theta);
  //lodpal_pairs reserved_pairsB = arp->pairs_info();

  double diffA = fabs(resultA.lod_score_with_cap - resultA.lod_score_without_cap);
  //double diffB = fabs(resultB.lod_score_with_cap - resultB.lod_score_without_cap);

  if( diffA > 0.75 )
  {
    if( my_method_D_called >= 2 )
    {
      my_built = false;

      inter_results invalid_result;

      invalid_result.lod_score_with_cap    = std::numeric_limits<double>::quiet_NaN();
      invalid_result.lod_score_without_cap = std::numeric_limits<double>::quiet_NaN();

      return invalid_result;
    }

    arp->pairs_info().pairs_info() = reserved_pairs.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairs.pairs_map();

    return run_method_D(init_theta);
  }
  else
  {
    arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();
    arp->update_result(resultA.param_estimates, removed_pos);

    return resultA;
  }
/*
  if( diffA > 0.75 && diffB > 0.75 )
  {
    if( my_method_D_called >= 2 )
    {
      my_built = false;

      inter_results invalid_result;

      invalid_result.lod_score_with_cap    = std::numeric_limits<double>::quiet_NaN();
      invalid_result.lod_score_without_cap = std::numeric_limits<double>::quiet_NaN();

      return invalid_result;
    }

    arp->pairs_info().pairs_info() = reserved_pairs.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairs.pairs_map();

    return run_method_D(init_theta);
  }
  else if( diffA <= 0.75 && diffB <= 0.75 )
  {
    if( fabs(resultA.lod_score_without_cap - resultB.lod_score_without_cap) <= 0.05 )
    {
      arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
      arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();
      arp->update_result(resultA.param_estimates, removed_pos);

      return resultA;
    }
    else
    {
      if( resultA.lod_score_without_cap > resultB.lod_score_without_cap )
      {
        arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
        arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();
        arp->update_result(resultA.param_estimates, removed_pos);

        return resultA;
      }
      else
      {
        arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
        arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();
        arp->update_result(resultB.param_estimates, removed_pos);

        return resultB;
      }
    }
  }
  else if( diffA <= 0.75 && diffB > 0.75 )
  {
    arp->pairs_info().pairs_info() = reserved_pairsA.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairsA.pairs_map();
    arp->update_result(resultA.param_estimates, removed_pos);

    return resultA;
  }
  else
  {
    arp->pairs_info().pairs_info() = reserved_pairsB.pairs_info();
    arp->pairs_info().pairs_map()  = reserved_pairsB.pairs_map();
    arp->update_result(resultB.param_estimates, removed_pos);

    return resultB;
  }
*/
}

void
lodpal_analysis::encode_maxfun_params(vector<maxfun_param_mgr>& init_theta)
{
  assert( arp->parameters().marker_count() == 1 );

  if( arp->previous_parameters().marker_count() )
  {
    init_theta.resize(4);
    arp->encode_params(0, init_theta[0]);
    arp->encode_params(1, init_theta[1]);
    arp->encode_params(2, init_theta[2]);
    arp->encode_params(3, init_theta[3]);
  }
  else
  {
    init_theta.resize(3);
    arp->encode_params(0, init_theta[0]);
    arp->encode_params(1, init_theta[1]);
    arp->encode_params(2, init_theta[2]);
  }
  
  return;
}

void
lodpal_analysis::decode_theta(const maxfun_results& re, vector<double>& params_estimate)
{
  for( int i = 0; i < re.getParameterMgr().getParamCount(); ++i )
  {
    double par = re.getParameterMgr().getParameter(i).getFinalEstimate();

    if( fabs(par) < 1.0e-10 )
      par = 0.0;

    params_estimate.push_back(par);
  }
}

void
lodpal_analysis::decode_final_result(const inter_results& final_result)
{
  arp->get_lodpal_result().clear();

  vector<double> first_derivatives;

  for( int i = 0; i < final_result.maxfun_info.getParameterMgr().getParamCount(); ++i )
  {
    double dev = final_result.maxfun_info.getParameterMgr().getParameter(i).getDeriv();

    if( fabs(dev) < 1.0e-10 )
      dev = 0.0;

    first_derivatives.push_back(dev);
  }

  arp->decode_theta(final_result.param_estimates, first_derivatives);

  arp->get_lodpal_result().set_maxfun_result(final_result.maxfun_info);
  arp->get_lodpal_result().set_method(final_result.method);
  arp->get_lodpal_result().set_lod_score(make_pair(final_result.lod_score_with_cap,
                                                   final_result.lod_score_without_cap));

  if( arp->pairs_info().relative_pairs().is_x_linked(arp->current_marker().marker) )
  {
    arp->get_lodpal_result().set_lod_score_mm(final_result.lod_score_mm);
    arp->get_lodpal_result().set_lod_score_mf(final_result.lod_score_mf);
    arp->get_lodpal_result().set_lod_score_ff(final_result.lod_score_ff);
  }

#if 0
  cout << "my_lod_score = " << final_result.maxfun_info.getFinalFunctionValue() << endl
       << "my_lod_score_mm = " << final_result.lod_score_mm << endl
       << "my_lod_score_mf = " << final_result.lod_score_mf << endl
       << "my_lod_score_ff = " << final_result.lod_score_ff << endl;
#endif

  return;
}

} // end of namespace LODPAL
} // end of namespace SAGE
