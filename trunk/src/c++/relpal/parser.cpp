#include "relpal/parser.h"

namespace SAGE   {
namespace RELPAL {

relpal_parser::relpal_parser(const relative_pairs& relpairs, cerrorstream& err)
             : BasicParser(err), my_pairs(relpairs)
{
  clear();
}

void
relpal_parser::clear()
{
  my_traits.resize(0);
  my_ind_covariates.resize(0);
  my_ind_batch_covariates.resize(0);
  my_ped_null_covariates.resize(0);
  my_ped_test_covariates.resize(0);
  my_null_markers.resize(0);
  my_test_markers.resize(0);
  my_ind_interactions.resize(0);
  my_ped_null_interactions.resize(0);
  my_ped_test_interactions.resize(0);

  my_reg_type = STSM;

  my_ind_batch_test  = false;
}

void
relpal_parser::parse_symbols(const SymbolTable* syms)
{
  return;
}

void
relpal_parser::parse_parameter(const LSFBase* param)
{
  return;
}

void
relpal_parser::parse_test_parameter_section(const LSFBase* params)
{
  if( !params || !params->List() )
    return;

  // Step 1. Parse trait & model type first.
  //
  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !(*i) || !(*i)->name().size() )
      continue;

    parse_test_parameter(*i);
  }

  // At least 1 trait statement required.
  if( !my_traits.size() )
  {
    return;
  }

  if( my_traits.size() > 1 )
  {
    if( my_reg_type == STSM )
    {
      my_reg_type = MTSM;
    }
    else if( my_reg_type == STMM )
    {
      my_reg_type = MTMM;
    }
    else if( my_reg_type == STZM )
    {
      my_reg_type = MTZM;
    }
  }

  // Step 2. Parse effects parameters & others.
  //
  for( i = params->List()->begin(); i!=params->List()->end(); ++i )
  {
    if( !(*i) || !(*i)->name().size() )
      continue;

    string name = toUpper( (*i)->name() );

    if( name == "FIRST_LEVEL" || name == "FIRST" )
      parse_effect_parameters(*i, true);
    else if( name == "SECOND_LEVEL" || name == "SECOND" )
      parse_effect_parameters(*i, false);
    else if( name == "DATA_OPTIONS" || name == "DATA" )
      parse_data_options(*i);
    else if( name == "OUTPUT_OPTIONS" || name == "OUTPUT" )
      parse_output_options(*i);
    else if( name == "PVALUE_OPTIONS" || name == "PVALUE" )
      parse_pvalue_options(*i);
  }

  // Step 3. Check # of markers, covariates, interaction vs. model.
  //

  // More than one test interaction - ignore them.
  //
  if( my_ped_test_interactions.size() > 1 )
  {
    if( my_reg_type == STZM ||my_reg_type == MTZM )
    {
      errors << priority(warning)
             << "Only one covariate by covariate interaction term is allowed for zero_marker model.  "
             << "Please specify only one covariate by covariate interaction statement.  "
             << "Interaction terms are ignored..." << endl;
    }
    else if( my_reg_type == STSM || my_reg_type == MTSM )
    {
      errors << priority(warning)
             << "Only one covariate by marker interaction term is allowed for single_marker model.  "
             << "Please specify only one covariate by marker interaction statement.  "
             << "Interaction terms are ignored..." << endl;
    }
    else
    {
      errors << priority(warning)
             << "Only one covariate by marker interaction term or "
             << "marker by marker interaction term is allowed for multiple_marker model.  "
             << "Please specify only one test interaction statements.  "
             << "Interaction terms are ignored..." << endl;
    }

    my_ped_test_interactions.clear();
  }

  //=============================
  // 1. One test interaction
  //=============================
  if( my_ped_test_interactions.size() )
  {
    vector<interaction_type>  ped_test_interactions;

    //--------------------------------------------------
    // 1.1. zero_marker model - no interaction allowed
    //--------------------------------------------------
    //--------------------------------------------------
    // 1.2 single_marker model - only cm interaction
    //--------------------------------------------------
    if( my_reg_type == STSM || my_reg_type == MTSM )
    {
      // Case 1 - one cov only, the other from test_markers or ibd file
      if( my_ped_test_interactions[0].batch )
      {
        process_batch_MC_interaction(ped_test_interactions);
      }
      // Case 2 - one cov & one marker
      else
      {
        ped_test_interactions.push_back(my_ped_test_interactions[0]);
      }

      update_null_covariates();
    }
    //--------------------------------------------------
    // 1.3 multiple_marker model - cm or mm interaction
    //--------------------------------------------------
    else if( my_reg_type == STMM || my_reg_type == MTMM )
    {
      // Case 1 - one cov or marker only, the other from test_markers or ibd file
      if( my_ped_test_interactions[0].batch )
      {
        if( my_ped_test_interactions[0].covariates.size() )
        {
          process_batch_MC_interaction(ped_test_interactions);
        }
        else
        {
          process_batch_MM_interaction(ped_test_interactions);
        }
      }
      // Case 2 - one cov & one marker or one marker & one marker
      else
      {
        ped_test_interactions.push_back(my_ped_test_interactions[0]);
      }

      update_null_markers();
      update_null_covariates();
    }

    my_ped_test_interactions.clear();

    for( size_t i = 0; i < ped_test_interactions.size(); ++i )
      my_ped_test_interactions.push_back(ped_test_interactions[i]);
  }
  //=============================
  // 2. No test interaction
  //=============================
  else
  {
    //--------------------------------------------------------------
    // 2.1 zero_marker model - test covariate vs. no test covariate
    //--------------------------------------------------------------
    if( my_reg_type == STZM || my_reg_type == MTZM )
    {
      if( my_ped_test_covariates.size() > 1 )
      {
        errors << priority(warning)
               << "More than one test covariates are specified in the zero_marker model.  "
               << "It is required to specify only one test covariate statement "
               << "for the zero_marker model.  "
               << "Please specify a proper number of second level covariate statements.  "
               << "Will perform sigle_marker model..." << endl;

        switch_to_single_marker_model();

        add_test_markers();
      }
      else if( my_ped_test_covariates.size() == 0 )
      {
        if( my_ind_batch_test )
        {
          for( size_t t = 0; t < my_pairs.trait_count(); ++t )
          {
            if(   !is_in_first_level_model(t)
                && my_pairs.fped_info().trait_info(t).usage() == RPED::RefTraitInfo::trait_covariate
                && my_pairs.fped_info().trait_info(t).name() != "SEX_CODE"
                && my_pairs.fped_info().trait_info(t).name() != "FAMILIAL_INDICATOR"
                && my_pairs.fped_info().trait_info(t).name() != "FOUNDER_INDICATOR"
                && my_pairs.fped_info().trait_info(t).name() != "PEDIGREE_SIZE" )
            {
              covariate_type cov;

              cov.covariate_index = t;
              cov.adj_trait_index = (size_t)-1;
              cov.valid           = true;

              my_ind_batch_covariates.push_back(cov);
            }
          }
        }

        if( my_ind_covariates.size() + my_ind_batch_covariates.size() == 0 )
        {
          errors << priority(warning)
                 << "No test covariates nor first level covariates are specified "
                 << "in the zero_marker model.  "
                 << "Will perform sigle_marker model..." << endl;

          switch_to_single_marker_model();

          add_test_markers();
        }
      }
    }
    //--------------------------------------------------------------
    // 2.2 single_marker model - test marker vs. no test marker
    //--------------------------------------------------------------
    else if( my_reg_type == STSM || my_reg_type == MTSM )
    {
      if( !my_test_markers.size() )
      {
        errors << priority(information)
               << "No Markers are specified for the single_marker model.  "
               << "Using all markers from the ibd file..." << endl;

        add_test_markers();
      }
    }
    //--------------------------------------------------------------
    // 2.3 multiple_marker model - test marker vs. no test marker
    //--------------------------------------------------------------
    else
    {
      if( my_test_markers.size() + my_null_markers.size() == 0 )
      {
        errors << priority(warning)
               << "No Markers are specified for the multiple_marker model.  "
               << "It is required to specify "
               << "at least one marker for the multiple_marker model.  "
               << "Please specify a proper number of marker statements for the model.  "
               << "Will perform the sigle_marker model..." << endl;

        switch_to_single_marker_model();

        add_test_markers();
      }
      else if( my_test_markers.size() + my_null_markers.size() == 1 )
      {
        if( !my_test_markers.size() )
        {
          // Add markers from ibd file as test markers.
          for( size_t m = 0; m < my_pairs.marker_count(); ++m )
          {
            if( my_null_markers[0].marker_index == m )
              continue;

            marker_type mr(m, true);
            mr.valid = true;

            my_test_markers.push_back(mr);
          }      
        }
        else
        {
          // Add markers from ibd file as null markers.
          for( size_t m = 0; m < my_pairs.marker_count(); ++m )
          {
            if( my_test_markers[0].marker_index == m )
              continue;

            marker_type mr(m, false);
            mr.valid = true;

            my_null_markers.push_back(mr);
          }      
        }
      }
      else if( my_test_markers.size() > 1 )
      {
        errors << priority(warning)
               << "More than one test markers are specified in the multiple_marker model.  "
               << "It is required to specify only one test marker statement "
               << "for the multiple_marker model.  "
               << "Please specify a proper number of marker statements.  "
               << "Will perform sigle_marker model..." << endl;

        switch_to_single_marker_model();
      }
      else if( !my_test_markers.size() )
      {
        errors << priority(warning)
               << "No test markers are specified for the multiple_marker model.  "
               << "It is required to specify one marker statement as test marker "
               << "among multiple marker statements for the multiple_marker model.  "
               << "Please specify a test marker statement for the model.  "
               << "Will perform the sigle_marker model..." << endl;

        switch_to_single_marker_model();

        add_test_markers();
      }
    }
  }

  return;
}

void
relpal_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size() ) return;

  string name = toUpper( param->name() );

  if( name == "TRAIT" )
    parse_trait(param);
  else if( name == "MODEL" || name == "TYPE" )
    parse_model_type(param);

  return;
}

void
relpal_parser::parse_trait(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t     = my_pairs.trait_find(tname);

    if( t < my_pairs.trait_count() )
    {
      bool binary = false;

      if( my_pairs.fped_info().trait_info(t).type() == RPED::RefTraitInfo::continuous_trait )
        binary = false;
      else if( my_pairs.fped_info().trait_info(t).type() == RPED::RefTraitInfo::binary_trait)
        binary = true;

      dependent_variable dv(t, binary);
      my_traits.push_back(dv);
    }
    else
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait parameter is missing name.  Skipping..." << endl;

  return;
}

void
relpal_parser::parse_model_type(const LSFBase* param)
{
  my_reg_type = STSM;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    if( toUpper(a.String()) == "SINGLE_MARKER" )
      my_reg_type = STSM;
    else if( toUpper(a.String()) == "MULTIPLE_MARKER" )
      my_reg_type = STMM;
    else if( toUpper(a.String()) == "ZERO_MARKER" )
      my_reg_type = STZM;
    else
      errors << priority(error)
             << "Unknown value for parameter 'model'.  "
             << "Will perform sigle marker regression analysis..." << endl;
  }

  return;
}

bool
relpal_parser::parse_effect_parameters(const LSFBase* params, bool first)
{
  if( !params || !params->List() ) return false;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !(*i) || !(*i)->name().size() )
      continue;

    string name = toUpper( (*i)->name() );

    if( name == "MARKER" && !first && !(my_reg_type == STZM || my_reg_type == MTZM) )
      parse_marker(*i);
    else if( name == "COVARIATE" )
      parse_covariate(*i, first);
    else if( name == "INTERACTION" )
    {
      if( my_traits.size() > 1 && !first )
      {
        errors << priority(warning)
               << "Found an interaction in secone_level with a multivariate model.  "
               << "Interaction option in second_level is allowed only for an univariate analysis.  "
               << "Please specify a proper model and the number of trait to use this option.  "
               << "Skipping interaction block..." << endl;
      }
      else
        parse_interaction(*i, first);
    }
    else if( name == "BATCH" && first && (my_reg_type == STZM || my_reg_type == MTZM) )
      my_ind_batch_test = true;
    else if( (name == "TRANSFORM" || name == "TRANSFORM_RESIDUAL" ||
              name == "NORMALIZE" || name == "NORMALIZE_RESIDUAL") && first )
      parse_boolean(*i, my_analysis_opt.transform_residuals);
    else if( (name == "NAIVE_VAR" || name == "NAIVE_VARIANCE") && !first )
      parse_boolean(*i, my_analysis_opt.naive_variance);
    else if( (name == "ROBUST_VAR"   || name == "ROBUST_VARIANCE" ||
              name == "SANDWICH_VAR" || name == "SANDWICH_VARIANCE") && !first )
      parse_boolean(*i, my_analysis_opt.sandwich_variance);
    else if( (name == "ALTER_VAR"       || name == "ALTER_VARIANCE" ||
              name == "ALTERNATIVE_VAR" || name == "ALTERNATIVE_VARIANCE") && !first )
      parse_boolean(*i, my_analysis_opt.alternative_variance);
    else if( (name == "IBD_VAR" || name == "IBD_VARIANCE") && !first )
    {
      parse_boolean(*i, my_analysis_opt.IBD_variance);

      if( my_analysis_opt.IBD_variance )
      {
        AttrVal a = attr_value((*i), "file");

        if( !a.has_value() || !a.String().size() )
          a = attr_value((*i), "state_file");
        if( !a.has_value() || !a.String().size() )
          a = attr_value((*i), "ibd_state_file");
        if( a.has_value() && a.String().size() )
        {
          my_analysis_opt.state_file_name = a.String();
        }
      }
    }
  }

  return true;
}

void
relpal_parser::parse_marker(const LSFBase* param)
{
  marker_type mar;

  parse_marker_type(param, mar);

  if( mar.valid )
  {
    if( mar.test_variable )
      my_test_markers.push_back(mar);
    else
      my_null_markers.push_back(mar);
  }

  return;
}

void
relpal_parser::parse_covariate(const LSFBase* param, bool first)
{
  covariate_type cov;

  parse_covariate_type(param, cov);

  if( cov.valid )
  {
    if( first )
      my_ind_covariates.push_back(cov);
    else
    {
      cov.adj_trait_index = (size_t)-1;
      if( cov.test_variable )
        my_ped_test_covariates.push_back(cov);
      else
        my_ped_null_covariates.push_back(cov);
    }
  }

  return;
}

void
relpal_parser::parse_interaction(const LSFBase* params, bool first)
{
  if( !params->List() )
  {
    errors << priority(warning)
           << "Found an empty interaction block.  "
           << "Please check your choices.  "
           << "Skipping interaction block..." << endl;

    return;
  }

  interaction_type interaction;

  parse_interaction_type(params, interaction, first);

  if( interaction.valid )
  {
    if( first )
      my_ind_interactions.push_back(interaction);
    else
    {
      if(   !(my_reg_type == STZM || my_reg_type == MTZM)
          && interaction.type == interaction_type::CC )
        my_ped_null_interactions.push_back(interaction);
      else
        my_ped_test_interactions.push_back(interaction);
    }
  }

  return;
}

void
relpal_parser::parse_marker_type(const LSFBase* param, marker_type& mar)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string mname = strip_ws(a.String());
    size_t m     = my_pairs.marker_find(mname);

    if( m < my_pairs.marker_count() )
    {
      mar.marker_index = m;

      // If single_marker model - test marker by default.
      //
      if( has_attr(param, "test") || my_reg_type == STSM || my_reg_type == MTSM )
        mar.test_variable = true;

      mar.valid = true;
    }
    else
      errors << priority(error) << "Marker '" << mname << "' not found.  Skipping..." << endl;
  }
  else
    errors << priority(error)
           << "Marker parameter is missing name.  Skipping..." << endl;

  return;
}

void
relpal_parser::parse_covariate_type(const LSFBase* param, covariate_type& cov)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string cname = strip_ws(a.String());
    size_t c     = my_pairs.trait_find(cname);

    if( c < my_pairs.trait_count() )
    {
      size_t adjust_t = (size_t)-1;

      // A specific trait adjustment only allowed for first level cov in multivariate model.
      //
      if( my_traits.size() > 1 )
      {
        a = attr_value(param, "adj_trait");
        if( a.has_value() )
        {
          string tname = strip_ws(a.String());
          size_t t     = my_pairs.trait_find(tname);

          if( t < my_pairs.trait_count() )
            adjust_t = t;
          else
            errors << priority(error)
                   << "Trait '" << tname
                   << "' for covariate adjustment not found.  The 'adj_trait' is ignored..." << endl;
        }
      }

      cov.covariate_index = c;
      cov.adj_trait_index = adjust_t;

      // test covariate accepted only for zero_marker model.
      //
      if( has_attr(param, "test") && (my_reg_type == STZM || my_reg_type == MTZM) )
        cov.test_variable = true;

      cov.valid = true;
    }
    else
      errors << priority(error) << "Covariate '" << cname << "' not found.  Skipping..." << endl;
  }
  else
    errors << priority(error)
           << "Covariate parameter is missing name.  Skipping..." << endl;

  return;
}

void
relpal_parser::parse_interaction_type(const LSFBase* params, interaction_type& interaction, bool first)
{
  if( !params->List() )
    return;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i!=params->List()->end(); ++i )
  {
    string name = toUpper( (*i)->name() );
    if( name == "COVARIATE" )
    {
      covariate_type cov;

      parse_covariate_type(*i, cov);

      if( cov.valid )
        interaction.covariates.push_back(cov);
    }
    else if( name == "MARKER" && !first )
    {
      marker_type mar;

      parse_marker_type(*i, mar);

      if( mar.valid )
        interaction.markers.push_back(mar);
    }
  }

  // Do some basic checking of the interaction sub-sblock.

  // No valid interaction statement.
  //
  if( interaction.markers.size() + interaction.covariates.size() == 0 )
  {
    errors << priority(warning)
           << "Found no valid marker or covariate statements in the interaction block.  "
           << "Please check your choices.  "
           << "Skipping interaction block..." << endl;

    interaction.valid = false;

    return;
  }

  // Only two-way allowed.
  //
  if( interaction.markers.size() + interaction.covariates.size() > 2 )
  {
    errors << priority(warning)
           << "Found more than two marker or covariate statements in the interaction block.  "
           << "Only two-way interaction is allowed.  Please check your choices.  "
           << "Skipping interaction block..." << endl;

    interaction.valid = false;

    return;
  }

  // Only Cov * Cov interaction in first level.
  //
  if( first )
  {
    if( interaction.covariates.size() != 2 )
    {
      errors << priority(warning)
             << "Found invalid number of covariate statements in the first level "
             << "interaction block.  "
             << "Two covariate statements need to be specified.  "
             << "Skipping interaction block..." << endl;

      interaction.valid = false;
    }
    else if( interaction.covariates[0] == interaction.covariates[1] )
    {
      errors << priority(warning)
             << "Found the same covariate statements in the first level "
             << "interaction block.  "
             << "Two different statements need to be specified.  "
             << "Skipping interaction block..." << endl;

      interaction.valid = false;
    }
    else
    {
      interaction.type  = interaction_type::CC;
      interaction.valid = true;
    }

    return;
  }

  // Second_level interaction.
  //
  bool valid = true;

  if(    interaction.covariates.size() == 2
      && interaction.covariates[0] == interaction.covariates[1] )
  {
    errors << priority(warning)
           << "Found the same covariate statements in the second level "
           << "interaction block.  "
           << "Two different statements need to be specified.  "
           << "Skipping interaction block..." << endl;

    valid = false;
  }
  else if(    interaction.markers.size() == 2
           && interaction.markers[0] == interaction.markers[1] )
  {
    errors << priority(warning)
           << "Found the same marker statements in the second level "
           << "interaction block.  "
           << "Two different statements need to be specified.  "
           << "Skipping interaction block..." << endl;

    valid = false;
  }
  else if( my_reg_type == STZM || my_reg_type == MTZM )
  {
    errors << priority(warning)
           << "Second level interaction is not allowed in a zero_marker model.  "
           << "Please check your choices.  "
           << "Skipping interaction block..." << endl;

    valid = false;
  }
  else if( my_reg_type == STSM || my_reg_type == MTSM )
  {
    if( interaction.markers.size() == 2 )
    {
      errors << priority(warning)
             << "Found two marker statemants in the interaction block.  "
             << "Marker by marker interaction can be included only in a multiple_marker model.  "
             << "Please check your choices.  "
             << "Skipping interaction block..." << endl;

      valid = false;
    }
    else if( interaction.covariates.size() + interaction.markers.size() == 1 )
    {
      if( has_attr(params, "batch") )
      {
        interaction.batch = true;
      }
      else
      {
        errors << priority(warning)
               << "Found only one valid marker or covariate statements in the interaction block.  "
               << "Two statements need to be specified.  "
               << "Skipping interaction block..." << endl;

        valid = false;
      }
    }
  }
  else
  {
    if( interaction.markers.size() + interaction.covariates.size() == 1 )
    {
      if( has_attr(params, "batch") )
      {
        interaction.batch = true;
      }
      else
      {
        errors << priority(warning)
               << "Found only one valid marker or covariate statements in the interaction block.  "
               << "Two statements need to be specified.  "
               << "Skipping interaction block..." << endl;

        valid = false;
      }
    }
  }

  if( valid )
  {
    interaction.valid = true;

    if( interaction.covariates.size() == 2 )
      interaction.type = interaction_type::CC;
    else if( interaction.covariates.size() == 1 )
      interaction.type = interaction_type::MC;
    else
      interaction.type = interaction_type::MM;
  }
  else
    interaction.valid = false;

  return;
}

void
relpal_parser::parse_data_options(const LSFBase* params)
{
  if( !params || !params->List() ) return;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !(*i) || !(*i)->name().size() )
      continue;

    string name = toUpper( (*i)->name() );

    if( name == "SUBSET" )
      parse_subset(*i);
    else if( name == "USE_PAIRS" || name == "USE_PAIR" )
      parse_use_pairs(*i);
    else if( name == "USE_MEMBERS" || name == "USE_MEMBER" || name == "MEMBER" )
      parse_use_members(*i);
  }

  return;
}

void
relpal_parser::parse_subset(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t     s = my_pairs.trait_find(tname);

    if( s < my_pairs.trait_count() )
    {
      filtering_type ft(s, tname);
      my_data_opt.subsets.push_back(ft);
    }
    else
      errors << priority(error)
             << "Trait '" << tname << "' not found for subset.  Skipping..." << endl;
  }
  else
    errors << priority(error)
           << "Subset parameter is missing name.  Skipping..." << endl;
}

void
relpal_parser::parse_use_pairs(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());

    if( n == "FULL" || n == "FSIB" )
    {
      my_data_opt.use_pairs = FSIB;
    }
    else if( n == "HALF" || n == "HSIB" )
    {
      my_data_opt.use_pairs = HSIB;
    }
    else if( n == "BOTH" || n == "SIB" )
    {
      my_data_opt.use_pairs = SIB;
    }
    else
    {
      errors << priority(warning)
             << "Unknown value for parameter 'use_pairs'.  Using all pairs..." << endl;

      my_data_opt.use_pairs = ALL;
    }
  }
  else
  {
    errors << priority(warning)
           << "No value specified for parameter 'use_pairs'.  Using all pairs..." << endl;

    my_data_opt.use_pairs = ALL;
  }

  return;
}

void
relpal_parser::parse_use_members(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());

    if( n == "INFORMATIVE_LOCAL" || n == "INF_LOCAL" || n == "LOCAL" )
    {
      my_data_opt.use_members = INFORMATIVE_LOCAL;
    }
    else if( n == "INFORMATIVE_REGION" || n == "INF_REGION" || n == "REGION" )
    {
      my_data_opt.use_members = INFORMATIVE_REGION;;
    }
    else if( n == "ALL" || n == "ALL_MEMBER" || n == "EVERY" )
    {
      my_data_opt.use_members = EVERY;
    }
    else
    {
      errors << priority(warning)
             << "Unknown value for parameter 'use_members'.  Using every members..." << endl;

      my_data_opt.use_members = DEFAULT;
    }
  }
  else
  {
    errors << priority(warning)
           << "No value specified for parameter 'use_members'.  Using every members..." << endl;

    my_data_opt.use_members = DEFAULT;
  }

  return;
}

void
relpal_parser::parse_output_options(const LSFBase* params)
{
  if( !params || !params->List() ) return;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !(*i) || !(*i)->name().size() )
      continue;

    string name = toUpper( (*i)->name() );

    if( name == "DETAILED_OUT" || name == "DETAILED" )
      parse_boolean(*i, my_output_opt.detailed_out);
    else if( name == "CSV_OUT" || name == "EXPORT_OUTPUT" || name == "EXPORT_OUT" || name == "EXPORT" )
      parse_boolean(*i, my_output_opt.export_out);
    else if( name == "DATA" || name == "DATA_OUT" || name == "DUMP_DATA" )
      parse_boolean(*i, my_output_opt.data_out);
    else if( name == "DEBUG" || name == "DEBUG_OUT" )
      parse_boolean(*i, my_output_opt.debug_out);
    else if( name == "RESIDUAL" || name == "RESIDUAL_OUT" )
      parse_boolean(*i, my_output_opt.residual_out);
  }

  return;
}

void
relpal_parser::parse_pvalue_options(const LSFBase* params)
{
  if( !params || !params->List() ) return;

  AttrVal a=attr_value(params, 0);

  if( a.has_value() )
  {
    parse_boolean(params, my_pvalue_opt.valid);
  }
  else
    my_pvalue_opt.valid = true;


  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    const LSFBase* param = *i;

    if( !param || !param->name().size() )
      continue;

    string name = toUpper( param->name() );

    a = attr_value(param, 0);

    if(  a.has_value() )
    {
      if( name == "SEED" )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.seed = a.Int();
        else
          errors << priority(error)
                 << "Invalid seed value specified for parameter 'empirical_pvalues'" << endl;
      }
      else if( name == "REPLICATES" )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.replicates = (size_t)a.Real();
        else
          errors << priority(error)
                 << "Invalid number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }
      else if( name == "MIN_REPLICATES" )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.min_replicates = (size_t)a.Real();
        else
          errors << priority(error)
                 << "Invalid minimum number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }

      else if( name == "MAX_REPLICATES" )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.max_replicates = (size_t)a.Real();
        else
          errors << priority(error)
                 << "Invalid maximum number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }
      else if( name == "THRESHOLD" )
      {
        if( finite(a.Real()) && a.Real() <= 1.0 && a.Real() > 0.0 )
          my_pvalue_opt.threshold = a.Real();
        else
          errors << priority(error)
                 << "Invalid threshold value specified for parameter 'empirical_pvalues'" << endl;
      }
      else if( name == "WIDTH" )
      {
        if( finite(a.Real()) && a.Real() > 0.0 )
          my_pvalue_opt.width = a.Real();
        else
          errors << priority(error)
                 << "Invalid width value specified for parameter 'empirical_pvalues'" << endl;
      }
      else if( name == "CONFIDENCE" )
      {
        if( finite(a.Real()) && a.Real() < 1.0 && a.Real() > 0.0 )
          my_pvalue_opt.confidence = a.Real();
        else
          errors << priority(error)
                 << "Invalid confidence value specified for parameter 'empirical_pvalues'" << endl;
      }
    }
  }

  return;
}

void
relpal_parser::make_MM_interaction_type(interaction_type&  new_inter,
                                        const marker_type& mar1,
                                        const marker_type& mar2) const
{
  new_inter.markers.push_back(mar1);
  new_inter.markers.push_back(mar2);

  new_inter.type  = interaction_type::MM;
  new_inter.valid = true;

  return;
}

void
relpal_parser::make_MC_interaction_type(interaction_type&     new_inter,
                                        const marker_type&    mar,
                                        const covariate_type& cov) const
{
  new_inter.markers.push_back(mar);
  new_inter.covariates.push_back(cov);

  new_inter.type  = interaction_type::MC;
  new_inter.valid = true;

  return;
}

void
relpal_parser::process_batch_MM_interaction(vector<interaction_type>& ped_test_interactions) const
{
  if( my_test_markers.size() )
  {
    for( size_t m = 0; m < my_test_markers.size(); ++m )
    {
      if( my_test_markers[m] == my_ped_test_interactions[0].markers[0] )
        continue;

      interaction_type mm_interaction;

      make_MM_interaction_type(mm_interaction, my_test_markers[m],
                               my_ped_test_interactions[0].markers[0]);

      ped_test_interactions.push_back(mm_interaction);
    }
  }
  else
  {
    for( size_t m = 0; m < my_pairs.marker_count(); ++m )
    {
      marker_type mr(m, false);
      mr.valid = true;

      if( mr == my_ped_test_interactions[0].markers[0] )
        continue;

      interaction_type mm_interaction;

      make_MM_interaction_type(mm_interaction, mr,
                               my_ped_test_interactions[0].markers[0]);

      ped_test_interactions.push_back(mm_interaction);
    }
  }

  return;
}

void
relpal_parser::process_batch_MC_interaction(vector<interaction_type>& ped_test_interactions) const
{
  if( my_test_markers.size() )
  {
    for( size_t m = 0; m < my_test_markers.size(); ++m )
    {
      interaction_type mc_interaction;

      make_MC_interaction_type(mc_interaction, my_test_markers[m],
                               my_ped_test_interactions[0].covariates[0]);

      ped_test_interactions.push_back(mc_interaction);
    }
  }
  else
  {
    for( size_t m = 0; m < my_pairs.marker_count(); ++m )
    {
      marker_type mr(m, false);
      mr.valid = true;

      interaction_type mc_interaction;

      make_MC_interaction_type(mc_interaction, mr,
                               my_ped_test_interactions[0].covariates[0]);

      ped_test_interactions.push_back(mc_interaction);
    }
  }

  return;
}

void
relpal_parser::update_null_markers()
{
  list<marker_type>  ped_null_markers;

  for( size_t m = 0; m < my_null_markers.size(); ++m )
  {
    bool main_effect_marker = false;
    for( size_t mi = 0; mi < my_ped_test_interactions[0].markers.size(); ++mi )
      if( my_null_markers[m] == my_ped_test_interactions[0].markers[mi] )
      {
        main_effect_marker = true;
        break;
      }

    if( !main_effect_marker )
      ped_null_markers.push_back(my_null_markers[m]);
  }

  // Remove duplicated entiries.
  ped_null_markers.unique();

  my_null_markers.clear();

  list<marker_type>::const_iterator mi = ped_null_markers.begin();
  for( ; mi != ped_null_markers.end(); ++mi )
    my_null_markers.push_back(*mi);

  return;
}

void
relpal_parser::update_null_covariates()
{
  list<covariate_type>  ped_null_covariates;

  // Add main effect covariate from test_interaction.
  for( size_t c = 0; c < my_ped_test_interactions[0].covariates.size(); ++c )
  {
    ped_null_covariates.push_back(my_ped_test_interactions[0].covariates[0]);
  }

  // Add null covariates.
  for( size_t c = 0; c < my_ped_null_covariates.size(); ++c )
  {
    ped_null_covariates.push_back(my_ped_null_covariates[c]);
  }

  // Add main effect covariates from null_interaction.
  for( size_t c = 0; c < my_ped_null_interactions.size(); ++c )
  {
    ped_null_covariates.push_back(my_ped_null_interactions[c].covariates[0]);
    ped_null_covariates.push_back(my_ped_null_interactions[c].covariates[1]);
  }

  // Remove duplicated entiries.
  ped_null_covariates.unique();

  my_ped_null_covariates.clear();

  list<covariate_type>::const_iterator ci = ped_null_covariates.begin();
  for( ; ci != ped_null_covariates.end(); ++ci )
    my_ped_null_covariates.push_back(*ci);

  return;
}

void
relpal_parser::add_test_markers()
{
  for( size_t m = 0; m < my_pairs.marker_count(); ++m )
  {
    marker_type mr(m, true);
    mr.valid = true;

    my_test_markers.push_back(mr);
  }

  return;
}
void
relpal_parser::switch_to_single_marker_model()
{
  if( my_traits.size() > 1 )
    my_reg_type = STSM;
  else
    my_reg_type = MTSM;

  return;
}

bool
relpal_parser::is_in_first_level_model(size_t batch_cov) const
{
  bool is_in = false;

  for( size_t t = 0; t < my_traits.size(); ++t )
  {
    if( batch_cov == my_traits[t].trait_index )
    {
      is_in = true;
      break;
    }
  }

  if( is_in )
    return true;

  for( size_t c = 0; c < my_ind_covariates.size(); ++c )
  {
    if( batch_cov == my_ind_covariates[c].covariate_index )
    {
      is_in = true;
      break;
    }
  }

  return is_in;
}


void
relpal_parser::dump_parser(ostream &out) const
{
  out << endl
      << "=================================" << endl
      << "  Relpal Analysis Specification  " << endl
      << "=================================" << endl << endl;

  //my_pairs.dump_pairs(out);
  //out << endl;

  out << "model: " << get_regression_type_to_string(my_reg_type)
      << endl << endl;

  out << "trait:" << endl;
  for( size_t i = 0; i < my_traits.size(); ++i )
    out << "   trait " << i+1 << ": " << my_traits[i].name(my_pairs)
        //<< "(" << my_traits[i].trait_index << ")"
        << endl;

  out << endl;

  out << "first_level:" << endl;
  for( size_t i = 0; i < my_ind_covariates.size(); ++i )
  {
    out << "  covariate " << i+1 << ": "
        << my_ind_covariates[i].name(my_pairs)
        //<< "(" << my_ind_covariates[i].covariate_index << ")"
        << ", adj_t = ";
    if( my_ind_covariates[i].adj_trait_index != (size_t)-1 )
      out << my_ind_covariates[i].adj_trait_index;
    else
      out << "all";
    //if( my_ind_covariates[i].test_variable )
    //  out << ", as test covariate";
    out << endl;
  }

  for( size_t i = 0; i < my_ind_batch_covariates.size(); ++i )
  {
    out << "  batch covariate " << i+1 << ": "
        << my_ind_batch_covariates[i].name(my_pairs)
        //<< "(" << my_ind_batch_covariates[i].covariate_index << ")"
        << ", adj_t = ";
    if( my_ind_batch_covariates[i].adj_trait_index!= (size_t)-1 )
      out << my_ind_batch_covariates[i].adj_trait_index;
    else
      out << "all";
    //if( my_ind_batch_covariates[i].test_variable )
    //  out << ", as test covariate";
    out << endl;
  }

  for( size_t i = 0; i < my_ind_interactions.size(); ++i )
  {
    out << "  interaction " << i << ": "
        << my_ind_interactions[i].name(my_pairs);
    out << endl;
  }

  out << endl;

  if( my_analysis_opt.transform_residuals )
    out << "transform residuals from 1st level to 2nd level" << endl << endl;

  out << "second_level:" << endl;
  for( size_t i = 0; i < my_ped_null_covariates.size(); ++i )
  {
    out << "  null covariate " << i+1 << ": "
        << my_ped_null_covariates[i].name(my_pairs);
        //<< "(" << my_ped_null_covariates[i].covariate_index << ")"
        //<< ", adj_t = " << my_ped_null_covariates[i].adj_trait_index;
    out << endl;
  }

  for( size_t i = 0; i < my_ped_test_covariates.size(); ++i )
  {
    out << "  test covariate " << i+1 << ": "
        << my_ped_test_covariates[i].name(my_pairs);
        //<< "(" << my_ped_test_covariates[i].covariate_index << ")"
        //<< ", adj_t = " << my_ped_test_covariates[i].adj_trait_index;
    //if( my_ped_test_covariates[i].test_variable )
    //  out << ", as test covariate";
    out << endl;
  }

  for( size_t i = 0; i < my_null_markers.size(); ++i )
  {
    out << "  null marker " << i+1 << ": " << my_null_markers[i].name(my_pairs);
        //<< "(" << my_null_markers[i].marker_index << ")";
    //if( !my_null_markers[i].test_variable )
    //  out << ", as null marker";
    out << endl;
  }

  for( size_t i = 0; i < my_test_markers.size(); ++i )
  {
    out << "  test marker " << i+1 << ": " << my_test_markers[i].name(my_pairs);
        //<< "(" << my_test_markers[i].marker_index << ")";
    //if( my_test_markers[i].test_variable )
    //  out << ", as test marker";
    out << endl;
  }

  for( size_t i = 0; i < my_ped_null_interactions.size(); ++i )
  {
    out << "  null interaction " << i+1 << ": "
        << my_ped_null_interactions[i].name(my_pairs);
    out << endl;
  }

  for( size_t i = 0; i < my_ped_test_interactions.size(); ++i )
  {
    out << "  test interaction " << i+1 << ": "
        << my_ped_test_interactions[i].name(my_pairs);
    out << endl;
  }

  out << endl;

  out << "data options:" << endl;
  out << my_data_opt.dump("  ") << endl;

  out << "output options:" << endl;
  out << my_output_opt.dump("  ") << endl;

  out << "asymptotic pvalue options:" << endl;
  out << my_pvalue_opt.dump("  ") << endl;

  if( !do_first_level_test() )
  {
    out << "variance options:" << endl;
    if( my_analysis_opt.naive_variance )
      out << "  naive variance = true" << endl;
    else
      out << "  naive variance = false" << endl;
    if( my_analysis_opt.sandwich_variance )
      out << "  robust sandwich variance = true" << endl;
    else
      out << "  robust sandwich variance = false" << endl;
    if( my_analysis_opt.alternative_variance )
      out << "  alternative variance = true" << endl;
    else
      out << "  alternative variance = false" << endl;
    if( my_analysis_opt.IBD_variance )
      out << "  IBD variance conditional on trait = true" << endl;
    else
      out << "  IBD variance conditional on trait = false" << endl;
    out << endl;
  }

  out << endl;

  return;
}

} // end of namespace RELPAL
} // end of namespace SAGE
