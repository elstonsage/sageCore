#include "sibpal/parser.h"

namespace SAGE   {
namespace SIBPAL {

sibpal_parser::sibpal_parser(const relative_pairs& pairs, cerrorstream& err)
            : BasicParser(err), my_pairs(pairs)
{}

//
//------------------------------------------------------------------------
//

meantest_parser::meantest_parser(const relative_pairs& pairs, cerrorstream& err)
               : sibpal_parser(pairs, err)
{
  clear();
}

meantest_parser::meantest_parser(const relative_pairs& pairs, SymbolTable* syms, cerrorstream& err)
               : sibpal_parser(pairs, err)
{
  clear();
  parse_symbols(syms);
}

meantest_parser::meantest_parser(const relative_pairs& pairs, LSFBase* params, cerrorstream& err)
               : sibpal_parser(pairs, err)
{
  clear();
  parse_test_parameter_section(params);
}

meantest_parser::meantest_parser(const relative_pairs& pairs, SymbolTable* syms, const LSFBase* params, cerrorstream& err)
               : sibpal_parser(pairs, err)
{
  clear();
  parse_symbols(syms);
  parse_test_parameter_section(params);
}

void
meantest_parser::clear()
{
  sibpal_parser::clear();

  my_wide_output = false;
  my_csv_output = false;
}

void
meantest_parser::parse_symbols(const SymbolTable* syms)
{
  if(!syms)
    return;

  AttrVal v = syms->query("wide_out");
  if(v.has_value())
  {
    if( toUpper(v.String()) == "TRUE" )
      set_wide_output(true);
    else if( toUpper(v.String()) == "FALSE" )
      set_wide_output(false);
  }

  v = syms->query("csv_out");
  if(v.has_value())
  {
    if( toUpper(v.String()) == "TRUE" )
      set_csv_output(true);
    else if( toUpper(v.String()) == "FALSE" )
      set_csv_output(false);
  }
}

void
meantest_parser::parse_parameter(const LSFBase* param)
{
  if( !param || !param->name().size()) return;
  AttrVal a;

  string name = toUpper( param->name() );

  if     (name == "WIDE_OUT")
    parse_boolean(param, my_wide_output);
  else if(name == "CSV_OUT")
    parse_boolean(param, my_csv_output);
}

void
meantest_parser::parse_test_parameter_section(const LSFBase* params)
{
  if(!params)
    return;

  if(!params->List())
    return;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
    parse_test_parameter(*i);
}

void
meantest_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size()) return;
  AttrVal a;

  string name = toUpper( param->name() );

  if(name == "MARKER")
    parse_marker(param);
  else if(name == "SUBSET")
    parse_subset(param);
  else if(name == "TRAIT")
    parse_trait(param);
  else if(name == "WIDE_OUT")
    parse_boolean(param, my_wide_output);
  else if(name == "CSV_OUT" || name == "EXPORT_OUTPUT" || name == "EXPORT_OUT")
  {
    parse_boolean(param, my_csv_output);
    my_parameters.set_export_output(my_csv_output);
  }
  else if(name == "PVAL_SCIENTIFIC_NOTATION" || name == "PVALUES_SCIENTIFIC_NOTATION")
  {
    bool sci_notation = false;
    parse_boolean(param, sci_notation);
    my_parameters.set_pvalues_scientific_notation(sci_notation);
  }
  else if(name == "W" || name == "W1")
  {
    double w = 0.5;
    parse_real(param, w);

    if( w < 0.0 || w > 0.5 )
    {
      errors << priority(warning)
             << "Invalid value for 'w1'.  The value for 'w1' should be in the range of (0.0, 0.5).  "
             << "Will use the default value of 'w1'."
             << endl;

      w = 0.5;
    }

    my_parameters.set_w(w);
  }
}

void
meantest_parser::parse_marker(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string mname = strip_ws(a.String());
    size_t m = my_pairs.marker_find(mname);
    if(m < my_pairs.marker_count())
      my_parameters.add_marker(m);
    else
      errors << priority(error) << "Marker '" << mname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Marker parameter is missing name.  Skipping..." << endl;
}

void
meantest_parser::parse_trait(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t = my_pairs.trait_find(tname);
    if(t < my_pairs.trait_count())
    {
      if(my_parameters.trait_set())
        errors << priority(warning)
               << "Only one trait at a time may be used in the mean ibd test."
               << "  Trait '" << tname << "' ignored..." << endl;
      else
        my_parameters.set_trait(t);
    }
    else
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait parameter is missing name.  Skipping..." << endl;
}

void
meantest_parser::parse_subset(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);
  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t = my_pairs.trait_find(tname);
    if(t < my_pairs.trait_count())
    {
      size_t affecteds = (size_t)-1;

      a=attr_value(param,"affecteds");
      if( a.has_value() )
      {
        if(my_pairs.fped_info().trait_info(t).type() == RPED::RefTraitInfo::binary_trait)
        {
          double aff = a.Real();
          if( finite(aff) && int(aff) >= 0 && int(aff) <= 2 )
            affecteds = int(aff);
          else
            errors << priority(warning)
                   << "Invalid number of affecteds in subset selction.  "
                   << "The number must be 0, 1 or 2." << endl;
        }
        else
          errors << priority(warning)
                 << "Cannot subset by affecteds using the continuous trait '"
                 << tname << "'." << endl;
      }

      if (affecteds <= 2)
        my_parameters.add_subset(t, affecteds);
      else
        my_parameters.add_subset(t);
    }
    else
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait parameter is missing name.  Skipping..." << endl;
}

//
//------------------------------------------------------------------------
//

regression_parser::regression_parser(const relative_pairs& pairs, cerrorstream& err)
                 : sibpal_parser(pairs, err)
{
  clear();
}

regression_parser::regression_parser(const relative_pairs& pairs, SymbolTable* syms, cerrorstream& err)
                 : sibpal_parser(pairs, err)
{
  clear();
  parse_symbols(syms);
}

regression_parser::regression_parser(const relative_pairs& pairs, LSFBase* params, cerrorstream& err)
                 : sibpal_parser(pairs, err)
{
  clear();
  parse_test_parameter_section(params);
}


regression_parser::regression_parser(const relative_pairs& pairs, SymbolTable* syms,
                                     const LSFBase* params, cerrorstream& err)
                 : sibpal_parser(pairs, err)
{
  clear();
  parse_symbols(syms);
  parse_test_parameter_section(params);
}

void
regression_parser::clear()
{
  sibpal_parser::clear();

  my_traits.resize(0);
  my_subsets.resize(0);
  my_covariates.resize(0);
  my_markers.resize(0);
  my_interactions.resize(0);
  my_batch_interactions.resize(0);

  my_reg_type = SINGLE_MARKER;
  my_reg_method_name = "prod";
}

void
regression_parser::parse_symbols(const SymbolTable* syms)
{
  return;
}

void
regression_parser::parse_parameter(const LSFBase* param)
{
  return;
}

void
regression_parser::parse_test_parameter_section(const LSFBase* params)
{
  if( !params || !params->List() )
    return;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !(*i) || !(*i)->name().size())
      continue;

    if( toUpper((*i)->name()) == "PAIR_INFO_FILE" )
      my_pair_info_file.push_back(*i);
    else
      parse_test_parameter(*i);
  }

  return;
}

void
regression_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size()) return;
  AttrVal a;

  string name = toUpper( param->name() );

  if( name == "TRAIT" )
    parse_trait(param);
  else if( name == "SUBSET" )
    parse_subset(param);
  else if( name == "MARKER" )
    parse_marker(param);
  else if( name == "COVARIATE" )
    parse_covariate(param);
  else if( name == "INTERACTION" )
    parse_interaction(param);
  else if( name == "REGRESSION_METHOD" || name == "TRAIT_REGRESSION_METHOD" )
    parse_regression_method(param);
  else if( name == "IDENTITY_WEIGHTS" || name == "IDENTITY_WEIGHT" )
  {
    bool iw = false;
    parse_boolean(param, iw);
    my_analysis_opt.identity_weight = iw;
  }
  else if( name == "POOL_CORRELATIONS" || name == "POOL_CORRELATION" )
  {
    bool pc = false;
    parse_boolean(param, pc);
    my_analysis_opt.pool_correlation = pc;
  }
  else if( name == "ROBUST_VARIANCE" )
  {
    bool r = false;
    parse_boolean(param, r);
    my_analysis_opt.robust_variance = r;

    bool adjust = false;
    parse_boolean(param, string("leverage_adjustment"), adjust);
    parse_boolean(param, string("leverage"), adjust);

    my_analysis_opt.leverage_adjustment = adjust;
  }
  else if( name == "WIDE_OUT" )
    parse_boolean(param, my_output_opt.wide_out);
  else if( name == "CSV_OUT" || name == "EXPORT_OUTPUT" || name == "EXPORT_OUT" )
    parse_boolean(param, my_output_opt.export_out);
  else if( name == "COMPUTE_EMPIRICAL_PVALUES"   || name == "EMPIRICAL_PVALUES"
       || name == "COMPUTE_PERMUTATION_PVALUES" || name == "PERMUTATION_PVALUES" )
    parse_empirical_pvalues(param);
  else if( name == "SKIP_UNINFORMATIVE_PAIRS" )
  {
    bool b = true;
    parse_boolean(param, b);

    my_data_opt.skip_uninformative_pairs = b;
  }
  else if( name == "PVALUES_SCIENTIFIC_NOTATION" || name == "PVALUE_SCIENTIFIC_NOTATION" )
  {
    bool ze = false;
    parse_boolean(param, ze);

    my_output_opt.pvalues_scientific_notation = ze;
  }
  else if( name == "PRINT_DESIGN_MATRIX" || name == "PRINT_DESIGN" )
    parse_design_matrix(param);
  else if( name == "PRINT_CORRELATION_MATRIX" || name == "PRINT_CORRELATION" )
    parse_correlation_matrix(param);
  else if( name == "USE_PAIRS" )
    parse_use_pairs(param);
  else if( name == "W" || name == "W1" )
  {
    double w = 0.5;
    parse_real(param, w);

    if( w < 0.0 || w > 0.5 )
    {
      errors << priority(warning)
             << "Invalid value for 'w1'.  The value for 'w1' should be in the range of (0.0, 0.5).  "
             << "Will use the default value of 'w1'."
             << endl;

      w = 0.5;
    }

    my_analysis_opt.w1 = w;
  }

  //
  // options not in document
  //
  else if( name == "DUMP_REGRESSION_DATA" || name == "DUMP_DATA" || name == "DATA_DUMP" ||
           name == "PRINT_QLS" || name == "QLS" || name == "DUMP_QLS" )
    parse_dump_data(param);
  else if( name == "DUMP_DEBUG" || name == "DEBUG_DUMP" )
    parse_dump_data(param);
  else if( name == "X_LINKAGE" || name == "X-LINKAGE" || name == "X" )
    parse_x_linkage(param);

  return;
}

void
regression_parser::check_test_options()
{
  if( my_data_opt.use_pairs != FSIB && my_pvalue_opt.is_on )
  {
    errors << priority(warning)
           << "Empirical P-value option is not allowed when half sib is used.  Skipping..."
           << endl;

    my_pvalue_opt.is_on = false;
  }

  if( my_data_opt.use_pairs == HSIB )
  {
    my_analysis_opt.sibship_mean = false;
    my_analysis_opt.blup_mean = false;
  }

  // Check batch interaction validity.
  //
  if( my_batch_interactions.size() > 1 )
  {
    errors << priority(warning)
           << "Only one batch interaction is allowed.  Skipping..." << endl;

    my_batch_interactions.resize(0);
  }
  else if( my_batch_interactions.size() == 1 )
  {
    if( get_regression_type() == SINGLE_MARKER )
    {
      if( my_batch_interactions[0].markers.size() )
      {
        errors << priority(warning)
               << "Marker by marker interaction is allowed only in "
               << "multiple_marker regression.  Skipping..."
               << endl;

        my_batch_interactions.resize(0);
      }
    }
    else if( get_regression_type() == ZERO_MARKER )
    {
      errors << priority(warning)
             << "Batch interaction option is not allowed in zero_marker regression.  Skipping..."
             << endl;

      my_batch_interactions.resize(0);
    }
    else if( my_batch_interactions[0].markers.size() )
      my_batch_interactions[0].markers[0].effect = marker_type::TOTAL;
  }

  // Check # of traits, markers vs. model.
  //
  if( !get_traits().size() )
  {
    for( size_t t = 0; t < my_pairs.trait_count(); ++t )
    {
      const RPED::RefTraitInfo &info = my_pairs.fped_info().trait_info(t);

      if( info.usage() != RPED::RefTraitInfo::trait_variate )
        continue;

      trait_type tt(t, info.type() == RPED::RefTraitInfo::binary_trait);

      my_traits.push_back(tt);
    }
  }

  // Add pair_covariate if exist
  //
  if( my_pairs.pair_covariate_count() )
  {
    for( size_t i = 0; i < my_pairs.pair_covariate_count(); ++i )
    {
      covariate_type cov(i, covariate_type::pair);
      my_covariates.push_back(cov);
    }
  }

  if(    get_regression_type() != ZERO_MARKER
      && !get_markers().size() )
  {
    for( size_t m = 0; m < my_pairs.marker_count(); ++m )
    {
      marker_type mt(m, marker_type::TOTAL, my_pairs.marker_genotype_model(m) == MLOCUS::X_LINKED);

      my_markers.push_back(mt);
    }

    if( get_regression_type() == MULTIPLE_MARKER && my_markers.size() == 1 )
    {
      errors << priority(warning)
             << "Only one marker present in ibd file.  "
             << "Changed regression type to single_marker_regression."
             << endl;

      my_reg_type = SINGLE_MARKER;
    }
  }

  // Make sure all markers are x_linked or not, no mixed allowed.
  //
  if( get_regression_type() == MULTIPLE_MARKER )
  {
    vector<size_t> x_markers;
    x_markers.resize(0);
    for( size_t i = 0; i < my_markers.size(); ++i )
    {
      if( my_markers[i].x_linked )
        x_markers.push_back(i);
    }

    if( x_markers.size() && x_markers.size() != my_markers.size() )
    {
      errors << priority(error)
             << "Regression on autosomal marker(s) and X-linked marker(s) "
             << "mixed not allowed.  Only using x_linked markers..."
             << "'" << endl;

      vector<marker_type> new_markers;
      for( size_t i = 0; i < x_markers.size(); ++i )
        new_markers.push_back(my_markers[x_markers[i]]);

      my_markers = new_markers;
    }

    if( my_markers.size() > 1 )
    {
      errors << priority(warning)
             << "At least 10 pairs per parameter are suggested for multiple regression."
             << endl;
    }
    else if( my_markers.size() == 1 )
    {
      errors << priority(warning)
             << "Only one marker present in the analysis block.  "
             << "Changed regression type to single_marker_regression."
             << endl;

      my_reg_type = SINGLE_MARKER;
    }
  }

  return;
}

void
regression_parser::parse_trait(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t     = my_pairs.trait_find(tname);

    if( t < my_pairs.trait_count() )
    {
      bool binary = false;

      if( my_pairs.fped_info().trait_info(t).type() == RPED::RefTraitInfo::binary_trait )
        binary = true;

      trait_type tt(t, binary);

      double        mean = param->attrs()->FloatAttr("mean");
      double   male_mean = param->attrs()->FloatAttr("male_mean");
      double female_mean = param->attrs()->FloatAttr("female_mean");
      double correlation = param->attrs()->FloatAttr("correlation");

      if( finite(mean) )
      {
        tt.fixed_mean = mean;
      }
      else
      {
        AttrVal am = attr_value(param, "mean");
        if( am.has_value() )
        {
          string m_option = strip_ws(am.String());
          if( toUpper(m_option)  == "SIBSHIP" || toUpper(m_option)  == "SIBSHIP_MEAN" )
          {
           my_analysis_opt.sibship_mean = true;
          }
          else if( toUpper(m_option)  == "BLUP" || toUpper(m_option)  == "BLUP_MEAN" )
          {
           my_analysis_opt.blup_mean = true;
          }
          else if( toUpper(m_option)  == "SAMPLE" || toUpper(m_option)  == "SAMPLE_MEAN" )
          {
           my_analysis_opt.blup_mean = false;
          }
        }
      }

      if( finite(male_mean) )
        tt.fixed_male_mean = male_mean;

      if( finite(female_mean) )
        tt.fixed_female_mean = female_mean;

      if( finite(correlation) )
        tt.fixed_correlation = correlation;

      my_traits.push_back(tt);
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
regression_parser::parse_subset(const LSFBase* param)
{
  AttrVal a=attr_value(param,0);

  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t s     = my_pairs.trait_find(tname);

    if( s < my_pairs.trait_count() )
    {
      trait_type st(s, true);
      my_subsets.push_back(st);
    }
    else
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Trait subset parameter is missing name.  Skipping..." << endl;

  return;
}

void
regression_parser::parse_marker(const LSFBase* param)
{
  marker_type mar;

  parse_marker_type(param, mar);

  if( mar.marker_index != (size_t)-1 )
  {
    my_markers.push_back(mar);
  }

  return;
}

void
regression_parser::parse_covariate(const LSFBase* param)
{
  covariate_type cov;

  parse_covariate_type(param, cov);

  if( cov.covariate_index != (size_t)-1 )
  {
    my_covariates.push_back(cov);
  }

  return;
}

void
regression_parser::parse_interaction(const LSFBase* params)
{
  if(!params->List())
    return;

  independent_variable interaction;
  interaction.type = independent_variable::INTERACTION;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    string name = toUpper( (*i)->name() );

    if( name == "COVARIATE" )
    {
      covariate_type cov;

      parse_covariate_type(*i, cov);

      if( cov.covariate_index != (size_t)-1 )
        interaction.covariates.push_back(cov);
    }
    else if( name == "MARKER" )
    {
      marker_type mar;

      parse_marker_type(*i, mar, true);

      if( mar.marker_index != (size_t)-1 )
        interaction.markers.push_back(mar);
    }
  }

  if( interaction.markers.size() + interaction.covariates.size() > 1 )
    my_interactions.push_back(interaction);
  else if( interaction.markers.size() + interaction.covariates.size() == 1 )
    my_batch_interactions.push_back(interaction);

  return;
}

void
regression_parser::parse_marker_type(const LSFBase* param, marker_type& mar, bool interaction_marker)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string mname = strip_ws(a.String());
    size_t m     = my_pairs.marker_find(mname);

    if( m < my_pairs.marker_count() )
    {
      mar.marker_index = m;
      mar.x_linked = (my_pairs.marker_genotype_model(m) == MLOCUS::X_LINKED);

      if( has_attr(param,"dom") || has_attr(param,"dominance") )
      {
        mar.effect = marker_type::DOMINANCE;
      }
      else if( interaction_marker && (has_attr(param,"add") || has_attr(param,"additive")) )
      {
        mar.effect = marker_type::ADDITIVE;
      }
      else
      {
        mar.effect = marker_type::TOTAL;
      }
    }
    else
      errors << priority(error) << "Marker '" << mname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Marker parameter is missing name.  Skipping..." << endl;

  return;
}

void
regression_parser::parse_covariate_type(const LSFBase* param, covariate_type& cov)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string cname = strip_ws(a.String());
    size_t c     = my_pairs.trait_find(cname);

    if( c < my_pairs.trait_count() )
    {
      cov.covariate_index = c;

      double cpower = param->attrs()->FloatAttr("power");

      if( !finite(cpower) || cpower == 0.0 )
        cpower = 1.0;

      cov.power = cpower;
      cov.fixed_mean = param->attrs()->FloatAttr("mean");

      cov.operation = covariate_type::prod;

      if( param->attrs()->has_attr("sum") )
      {
        cov.operation = covariate_type::sum;
      }
      if( param->attrs()->has_attr("diff") )
      {
        cov.operation = covariate_type::diff;
      }
      if( param->attrs()->has_attr("avg") )
      {
        cov.operation = covariate_type::avg;
      }
      if( param->attrs()->has_attr("single") )
      {
        cov.operation = covariate_type::single;
      }
      if( param->attrs()->has_attr("all") )
      {
        cov.operation = covariate_type::both;
      }
    }
    else
      errors << priority(error) << "Covariate '" << cname << "' not found." << endl;
  }
  else
    errors << priority(error)
           << "Covariate parameter is missing name.  Skipping..." << endl;

  return;
}

void
regression_parser::parse_regression_method(const LSFBase* param)
{
  AttrVal a=attr_value(param, 0);

  if( a.has_value() )
  {
    if( a.String().size() )
      my_reg_method_name = a.String();
    else
      errors << priority(error) << "Unknown value for parameter 'regression_method'" << endl;
  }

  return;
}

void
regression_parser::parse_use_pairs(const LSFBase* param)
{
  AttrVal a=attr_value(param, 0);

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
    else if( n == "BOTH" || n == "ALL" )
    {
      my_data_opt.use_pairs = SIB;
    }
    else
      errors << priority(warning)
             << "Unknown value for parameter 'use_pairs'.  Using full sib-pairs..." << endl;
  }
  else
  {
    errors << priority(warning)
           << "No value specified for parameter 'use_pairs'.  Using full sib-pairs..." << endl;

    my_data_opt.use_pairs = SIB;
  }

  return;
}

void
regression_parser::parse_empirical_pvalues(const LSFBase* param)
{
  AttrVal a=attr_value(param, 0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());
    if( n == "TRUE" || n == "YES" )
      my_pvalue_opt.is_on = true;
    else if( n == "FALSE" || n == "NO" )
      my_pvalue_opt.is_on = false;
    else
      errors << priority(error) << "Unknown value for parameter 'compute_empirical_pvalues'" << endl;
  }

  if( my_pvalue_opt.is_on )
  {
    a = attr_value(param,"permutations");

    if( a.has_value() )
    {
      if( finite(a.Real()) )
        my_pvalue_opt.replicates = (size_t)a.Real();
      else
        errors << priority(error) << "Invalid number of permutations specified for parameter 'empirical_pvalues'" << endl;
    }
    else
    {
      a = attr_value(param,"replicates");
      if( a.has_value() )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.replicates = (size_t)a.Real();
        else
          errors << priority(error) << "Invalid number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }
    }

    a = attr_value(param,"min_permutations");
    if( a.has_value() )
    {
      if( finite(a.Real()) )
        my_pvalue_opt.min_replicates = (size_t)a.Real();
      else
        errors << priority(error) << "Invalid minimum number of permutations specified for parameter 'empirical_pvalues'" << endl;
    }
    else
    {
      a = attr_value(param,"min_replicates");
      if( a.has_value() )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.min_replicates = (size_t)a.Real() ;
        else
          errors << priority(error) << "Invalid minimum number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }
    }

    a = attr_value(param,"max_permutations");
    if( a.has_value() )
    {
      if( finite(a.Real()) )
        my_pvalue_opt.max_replicates = (size_t)a.Real();
      else
        errors << priority(error) << "Invalid maximum number of permutations specified for parameter 'empirical_pvalues'" << endl;
    }
    else
    {
      a = attr_value(param,"max_replicates");
      if( a.has_value() )
      {
        if( finite(a.Real()) )
          my_pvalue_opt.max_replicates = (size_t)a.Real();
        else
          errors << priority(error) << "Invalid maximum number of replicates specified for parameter 'empirical_pvalues'" << endl;
      }
    }

    a = attr_value(param,"threshold");
    if( a.has_value() )
    {
      if( finite(a.Real()) && a.Real() <= 1.0 && a.Real() > 0.0 )
        my_pvalue_opt.threshold = a.Real() ;
      else
        errors << priority(error) << "Invalid threshold value specified for parameter 'empirical_pvalues'" << endl;
    }

    a = attr_value(param,"seed");
    if( a.has_value() )
    {
      if(finite(a.Real()))
        my_pvalue_opt.seed = a.Int();
      else
        errors << priority(error) << "Invalid seed value specified for parameter 'empirical_pvalues'" << endl;
    }

    a = attr_value(param,"width");
    if( a.has_value() )
    {
      if( finite(a.Real()) && a.Real() > 0.0 )
        my_pvalue_opt.width = a.Real();
      else
        errors << priority(error) << "Invalid width value specified for parameter 'empirical_pvalues'" << endl;
    }

    a = attr_value(param,"confidence");
    if( a.has_value() )
    {
      if( finite(a.Real()) && a.Real() < 1.0 && a.Real() > 0.0 )
        my_pvalue_opt.confidence = a.Real();
      else
        errors << priority(error) << "Invalid confidence value specified for parameter 'empirical_pvalues'" << endl;
    }
  }

  return;
}

void
regression_parser::parse_dump_data(const LSFBase* param)
{
  if( my_reg_type == ZERO_MARKER )
  {
    my_output_opt.dump_data_file = "";
    my_output_opt.dump_data = false;
  }
  else if( my_reg_type == MULTIPLE_MARKER )
  {
    my_output_opt.dump_data_file = "traits";
    my_output_opt.dump_data = true;
  }
  else
  {
    AttrVal a = attr_value(param,0);

    if( a.has_value() )
    {
      string f = "";
      parse_string(param, f);

      size_t m = my_pairs.marker_find(f);
      if( m < my_pairs.marker_count() )
      {
        my_output_opt.dump_data_file = f;
        my_output_opt.dump_data = true ;
      }
      else
      {
        errors << priority(warning)
               << "Invalid value specified for parameter 'print_qls'.  Skipping..." << endl;

        my_output_opt.dump_data = false;
      }
    }
    else
    {
      errors << priority(warning)
             << "No value specified for parameter 'print_qls'.  Skipping..." << endl;

      my_output_opt.dump_data = false;
    }
  }

  return;
}

void
regression_parser::parse_design_matrix(const LSFBase* param)
{
  AttrVal a = attr_value(param,"rows");

  if( a.has_value() )
  {
    if( finite(a.Real()) )
      my_output_opt.design_matrix_rows = (size_t)a.Real();
    else if( toUpper(a.String()) == "ALL" )
      my_output_opt.design_matrix_rows = size_t(-1);
    else
    {
      errors << priority(warning)
             << "Invalid number of rows specified for parameter 'print_design_matrix'.  "
             << "Will print 10 rows by default." << endl;
    }
  }

  if( my_reg_type != SINGLE_MARKER )
  {
    my_output_opt.design_matrix_file = "traits" ;
    my_output_opt.print_design_matrix = true ;
  }
  else
  {
    a = attr_value(param,0);

    if( a.has_value() )
    {
      string f = "";
      parse_string(param, f);

      size_t m = my_pairs.marker_find(f);
      if( m < my_pairs.marker_count() )
      {
        my_output_opt.design_matrix_file = f;
        my_output_opt.print_design_matrix = true ;
      }
      else
      {
        errors << priority(warning)
               << "Invalid value specified for parameter 'print_design_matrix'.  "
               << "Will print for all locations." << endl;

        my_output_opt.design_matrix_file = "traits" ;
        my_output_opt.print_design_matrix = true ;
      }
    }
    else
    {
      errors << priority(warning)
             << "No value specified for parameter 'print_design_matrix'.  "
             << "Will print for all locations." << endl;

      my_output_opt.design_matrix_file = "traits" ;
      my_output_opt.print_design_matrix = true ;
    }
  }
}

void
regression_parser::parse_correlation_matrix(const LSFBase* param)
{
  if( my_reg_type != SINGLE_MARKER )
  {
    my_output_opt.correl_matrix_file = "traits" ;
    my_output_opt.print_correl_matrix = true ;
  }
  else
  {
    AttrVal a = attr_value(param,0);

    if( a.has_value() )
    {
      string f = "";
      parse_string(param, f);

      size_t m = my_pairs.marker_find(f);
      if( m < my_pairs.marker_count() )
      {
        my_output_opt.correl_matrix_file = f ;
        my_output_opt.print_correl_matrix = true ;
      }
      else
      {
        errors << priority(warning)
               << "Invalid value specified for parameter 'print_correlation_matrix'.  Skipping..." << endl;

        //my_output_opt.correl_matrix_file = "traits" ;
        my_output_opt.print_correl_matrix = false ;
      }
    }
    else
    {
      errors << priority(warning)
             << "No value specified for parameter 'print_correlation_matrix'.  Skipping..." << endl;

      //my_output_opt.correl_matrix_file = "traits" ;
      my_output_opt.print_correl_matrix = false ;
    }
  }
}

void
regression_parser::parse_x_linkage(const LSFBase* params)
{
  if(!params)
    return;
  
  if(!params->List())
    return;

  bool pair_type_chosen = false;
  bool mm_pair = false;
  bool mf_pair = false;
  bool ff_pair = false;

//  bool sex_specific_mean = false;
//  double male_mean = 0.0;
//  double female_mean = 0.0;

  LSFList::const_iterator i;
  for(i=params->List()->begin(); i!=params->List()->end(); ++i)
  {
    if( !(*i) || !(*i)->name().size())
      continue;

    const LSFBase* param = *i;

    if( toUpper((*i)->name()) == "PAIR_TYPE" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        pair_type_chosen = true;
          
        if( toUpper(a.String()) == "M_M" || toUpper(a.String()) == "M-M" || toUpper(a.String()) == "MM")
          mm_pair = true;
        else if( toUpper(a.String()) == "M_F" || toUpper(a.String()) == "M-F" || toUpper(a.String()) == "MF")
          mf_pair = true;
        else if( toUpper(a.String()) == "F_F" || toUpper(a.String()) == "F-F" || toUpper(a.String()) == "FF")
          ff_pair = true;
        else
        {
          mm_pair = true;
          mf_pair = true;
          ff_pair = true;
        }
      }
    }
/*
    else if(    toUpper((*i)->name()) == "TRAIT_MEAN_SEX_SPECIFIC"
             || toUpper((*i)->name()) == "MEAN_SEX_SPECIFIC"
             || toUpper((*i)->name()) == "SEX_SPECIFIC_MEAN" )
    {
      AttrVal a=attr_value(param,0);
      if( a.has_value() )
      {
        if( toUpper(a.String()) == "TRUE" || toUpper(a.String()) == "YES" )
        {
          sex_specific_mean = true;

          if( param->attrs()->has_attr("male") )
          {
            male_mean = param->attrs()->FloatAttr("male");
            if( !finite(male_mean) )
              male_mean = NaN;
          }

          if( param->attrs()->has_attr("female") )
          {
            female_mean = param->attrs()->FloatAttr("female");
            if( !finite(female_mean) )
              female_mean = NaN;
          }
        }
      }
    }
*/
  }

  if( !pair_type_chosen )
  {
    mm_pair = true;
    mf_pair = true;
    ff_pair = true;
  }

  my_data_opt.use_mm_pair = mm_pair;
  my_data_opt.use_mf_pair = mf_pair;
  my_data_opt.use_ff_pair = ff_pair;
}

void
regression_parser::dump_parser(ostream &out) const
{
  out << endl 
      << "===========================================" << endl
      << "  Trait Regression Analysis Specification  " << endl
      << "===========================================" << endl << endl;

  //my_pairs.dump_pairs(out);
  //out << endl;

  out << "regression type: " << get_regression_type_to_string(my_reg_type)
      << endl << endl;

  out << "regression method: " << get_regression_method_name()
      << endl << endl;

  out << "trait:" << endl;
  for( size_t i = 0; i < my_traits.size(); ++i )
    out << "   trait " << i+1 << ": " << my_traits[i].name(my_pairs)
        //<< "(" << my_traits[i].trait_index << ")"
        << endl;
  out << endl;

  if( my_subsets.size() )
  {
    for( size_t i = 0; i < my_subsets.size(); ++i )
      out << "  subset " << i << ": " << my_subsets[i].name(my_pairs) << endl;
    out << endl;
  }

  out << "independent params:" << endl;

  for( size_t i = 0; i < my_covariates.size(); ++i )
  {
    out << "  covariate " << i+1 << ": " << my_covariates[i].name(my_pairs)
        << "(" << my_covariates[i].short_effect_name() << ")";
    out << endl;
  }

  for( size_t i = 0; i < my_markers.size(); ++i )
  {
    out << "  marker " << i+1 << ": " << my_markers[i].name(my_pairs)
        << "(" << my_markers[i].short_effect_name() << ")";
    out << endl;
  }

  for( size_t i = 0; i < my_interactions.size(); ++i )
  {
    out << "  interaction " << i+1 << ": " << my_interactions[i].name(my_pairs)
        << "(" << my_interactions[i].short_effect_name() << ")";
    out << endl;
  }

  out << endl;

  out << "data options:" << endl;
  out << my_data_opt.dump("  ") << endl;

  out << "analysis options:" << endl;
  out << my_analysis_opt.dump("  ") << endl;

  out << "empirical pvalue options:" << endl;
  out << my_pvalue_opt.dump("  ") << endl;

  out << "output options:" << endl;
  out << my_output_opt.dump("  ") << endl;

  out << endl;

  return;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
