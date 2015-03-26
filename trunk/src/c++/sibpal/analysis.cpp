#include "sibpal/analysis.h"

namespace SAGE   {
namespace SIBPAL {

sibpal_analysis::sibpal_analysis(cerrorstream& err)
              : errors(err)
{ }

void
sibpal_analysis::run_mean_test(relative_pairs&    pairs,
                               meantest_parser&   mean_parser,
                               ostream&           out,
                               ostream&           csv)
{
  if(!out) return;

  meantest_parameters params = mean_parser.parameters();

  bool wide_output  = mean_parser.wide_output();
  bool csv_output   = mean_parser.csv_output();

#if 0
  bool default_rest = false;

  if( !params.marker_count()  && !params.trait_set() )
    default_test = true;
#endif

  if( !params.marker_count() )
  {
    for(size_t m=0; m < pairs.marker_count(); ++m)
      params.add_marker(m);
  }

  SibMeanTest mtest(pairs, errors);
  mtest.set_parameters(params);

  typedef boost::shared_ptr<meantest_viewer> ro_ptr;
  ro_ptr reg_output;
  ro_ptr expo_output;

  reg_output = ro_ptr( new meantest_textfile(out, pairs, wide_output, false) );

  if( csv_output )
    expo_output = ro_ptr( new meantest_csvfile(csv, pairs, wide_output, false)  );

  if( !mtest.parameters().trait_set() )
  {  
    if( pairs.fsib_pair_count() )
    {
      mtest.set_use_pairs(make_pair(true, false));
      mtest.regress();

      if(!mtest.parameters().valid())
      {
        errors << priority(error) << "Marker Regression: Regression not successfull." << endl;
        return;
      }

      reg_output->print_results( mtest.parameters() );

      if( csv_output )
        expo_output->print_results( mtest.parameters() );
    }

    if( pairs.hsib_pair_count() )
    {
      mtest.set_use_pairs(make_pair(false, true));
      mtest.regress();

      if(!mtest.parameters().valid())
      {
        errors << priority(error) << "Marker Regression: Regression not successfull." << endl;
        return;
      }

      reg_output->print_results( mtest.parameters() );

      if( csv_output )
        expo_output->print_results( mtest.parameters() );
    }
  }
  else
  {
    if( pairs.fsib_pair_count() )
    {
      mtest.set_use_pairs(make_pair(true, false));

      size_t t = mtest.parameters().trait().trait;

      mtest.parameters().set_trait(t);

      vector<meantest_parameters> aff_mreg_params;
      for( int a = 0; a <= 2; ++a )
      {
        mtest.parameters().set_affection_status(a);
        mtest.regress();
        aff_mreg_params.push_back( mtest.parameters() );
      }

      reg_output->print_results( aff_mreg_params, pairs.trait_name(t) );

      if( csv_output )
        expo_output->print_results( aff_mreg_params, pairs.trait_name(t) );
    }

    if( pairs.hsib_pair_count() )
    {
      mtest.set_use_pairs(make_pair(false, true));

      size_t t = mtest.parameters().trait().trait;

      mtest.parameters().set_trait(t);

      vector<meantest_parameters> aff_mreg_params;
      for( int a = 0; a <= 2; ++a )
      {
        mtest.parameters().set_affection_status(a);
        mtest.regress();
        aff_mreg_params.push_back( mtest.parameters() );
      }

      reg_output->print_results( aff_mreg_params, pairs.trait_name(t) );

      if( csv_output )
        expo_output->print_results( aff_mreg_params, pairs.trait_name(t) );
    }
  }
}

void
sibpal_analysis::run_regression_analysis(relative_pairs&    r_pairs,
                                         regression_parser& r_parser,
                                         ostream&           out,
                                         ostream&           det,
                                         ostream&           exp)
{
  if( !out ) return;

  my_pairs = &r_pairs;

  r_parser.check_test_options();

  build_models(r_parser);

#if 0
  cout << "number of models = " << my_analysis_models.size() << endl;
#endif

  if( !my_analysis_models.size() )
  {
    errors << priority(error)
           << "No valid regression analysis model exist!  "
           << "Skipping..."
           << endl;

    return;
  }

  // Output
  //
  typedef boost::shared_ptr<regression_outfile> to_ptr;

  to_ptr reg_output = to_ptr( new regression_outfile(out, det, exp, errors) );

  // Start regression
  //
  TraitRegression reg(r_pairs, errors);
  string previous_trait = "";
  size_t a_index = 0;

  for( size_t i = 0; i < my_analysis_models.size(); ++i, ++a_index )
  {
    my_current_model = my_analysis_models[i];

    //my_current_model.dump_model(r_pairs, cout);
    string current_trait = my_current_model.get_trait().name(r_pairs);

    reg.set_model(my_current_model);

    reg.regress();

    if( current_trait != previous_trait )
    {
      if( i )
        reg_output->print_footer();

      reg_output->print_header(reg);
      a_index = 0;
    }

    if( reg.get_svd_return_code() == 2 )
      reg_output->print_instability_note();

    reg_output->print_results(a_index+1);

    if( i == my_analysis_models.size()-1 )
      reg_output->print_footer();

    previous_trait = current_trait;
  }

  return;
}

//
//-------------------------------------------------------------------------
//

void
sibpal_analysis::build_models(const regression_parser& r_parser)
{
  my_analysis_models.resize(0);

  if( r_parser.get_regression_type() == ZERO_MARKER )
  {
    for( size_t t = 0; t < r_parser.get_traits().size(); ++t )
    {
      regression_model current_model;

      add_a_trait(current_model, r_parser.get_traits()[t]);

      add_covariates(current_model, r_parser);
      add_interactions(current_model, r_parser);
      add_other_options(current_model, r_parser);

      if( is_valid_model(current_model) )
      {
        current_model.validate();
        my_analysis_models.push_back(current_model);
      }
    }
  }
  else if( r_parser.get_regression_type() == MULTIPLE_MARKER )
  {
    for( size_t t = 0; t < r_parser.get_traits().size(); ++t )
    {
      if( r_parser.get_batch_interactions().size() )
      {
        for( size_t m = 0; m < r_parser.get_markers().size(); ++m )
        {
          regression_model current_model;

          add_a_trait(current_model, r_parser.get_traits()[t]);
          add_a_marker(current_model, r_parser.get_markers()[m]);

          add_covariates(current_model, r_parser);
          add_batch_interactions(current_model, r_parser, r_parser.get_markers()[m], false);
          add_other_options(current_model, r_parser);

          if( is_valid_model(current_model) )
          {
            current_model.validate();
            my_analysis_models.push_back(current_model);
          }
        }
      }
      else
      {
        regression_model current_model;

        add_a_trait(current_model, r_parser.get_traits()[t]);

        add_markers(current_model, r_parser);

        add_covariates(current_model, r_parser);
        add_interactions(current_model, r_parser);
        add_other_options(current_model, r_parser);

        if( is_valid_model(current_model) )
        {
          current_model.validate();
          my_analysis_models.push_back(current_model);
        }
      }
    }
  }
  else
  {
    for( size_t t = 0; t < r_parser.get_traits().size(); ++t )
    {
      for( size_t m = 0; m < r_parser.get_markers().size(); ++m )
      {
        regression_model current_model;

        add_a_trait(current_model, r_parser.get_traits()[t]);

        add_a_marker(current_model, r_parser.get_markers()[m]);

        add_covariates(current_model, r_parser);

        if( r_parser.get_batch_interactions().size() )
          add_batch_interactions(current_model, r_parser, r_parser.get_markers()[m], true);

        add_interactions(current_model, r_parser);
        add_other_options(current_model, r_parser);

        if( is_valid_model(current_model) )
        {
          current_model.validate();
          my_analysis_models.push_back(current_model);
        }
      }
    }
  }

  set_output_options(r_parser);

  return;
}


void
sibpal_analysis::add_markers(regression_model&        a_model,
                             const regression_parser& r_parser)
{
  for( size_t m = 0; m < r_parser.get_markers().size(); ++m )
    add_a_marker(a_model, r_parser.get_markers()[m]);

  return;
}

void
sibpal_analysis::add_covariates(regression_model&        a_model,
                                const regression_parser& r_parser)
{
  for( size_t c = 0; c < r_parser.get_covariates().size(); ++c )
    add_a_covariate(a_model, r_parser.get_covariates()[c]);

  return;
}

void
sibpal_analysis::add_interactions(regression_model&        a_model,
                                  const regression_parser& r_parser)
{
  for( size_t i = 0; i < r_parser.get_interactions().size(); ++i )
  {
    const independent_variable& inter = r_parser.get_interactions()[i];

    bool main_exist = check_main_effects(a_model, inter);

    if( main_exist )
      a_model.add_parameter(inter);
    else
    {
      errors << priority(warning)
             << "Interaction terms are allowed only when both corresponding main effects "
             << "are included in the model.  Skipping interaction term '"
             << inter.name(*my_pairs)
             << "(" << inter.effect_name()
             << ")'." << endl;
    }
  }

  return;
}

void
sibpal_analysis::add_batch_interactions(regression_model&        a_model,
                                        const regression_parser& r_parser,
                                        const marker_type&       mar,
                                        bool                     allow_dom)
{
  if( !allow_dom && mar.effect == marker_type::DOMINANCE )
    return;

  const independent_variable& b_inter = r_parser.get_batch_interactions()[0];

  independent_variable inter;

  inter.markers.push_back(mar);

  if( b_inter.markers.size() )
  {
    inter.markers.push_back(b_inter.markers[0]);

    if( !check_main_marker_effect(a_model, b_inter.markers[0]) )
    {
      add_a_marker(a_model, b_inter.markers[0]);
    }
  }
  else
  {
    inter.covariates.push_back(b_inter.covariates[0]);

    if( !check_main_covariate_effect(a_model, b_inter.covariates[0]) )
    {
      add_a_covariate(a_model, b_inter.covariates[0]);
    }
  }

  bool main_exist = check_main_effects(a_model, inter);

  if( main_exist )
  {
    a_model.add_parameter(inter);

    if( mar.effect == marker_type::DOMINANCE )
    {
      marker_type mar_add(mar.marker_index, marker_type::ADDITIVE, mar.x_linked);

      independent_variable inter_add;
      inter_add.markers.push_back(mar_add);
      inter_add.covariates.push_back(inter.covariates[0]);

      a_model.add_parameter(inter_add);
    }
  }
  else
  {
    errors << priority(warning)
           << "Interaction terms are allowed only when both corresponding main effects "
           << "are included in the model.  Skipping interaction term '"
           << inter.name(*my_pairs)
           << "(" << inter.effect_name()
           << ")'." << endl;
  }

  return;
}

void
sibpal_analysis::add_other_options(regression_model&        a_model,
                                   const regression_parser& r_parser)
{
  for( size_t s = 0; s < r_parser.get_subsets().size(); ++s )
    a_model.add_subset( r_parser.get_subsets()[s] );

  a_model.set_data_options(r_parser.get_data_options());
  a_model.set_analysis_options(r_parser.get_analysis_options());
  a_model.set_pvalue_options(r_parser.get_pvalue_options());

  a_model.set_regression_type(r_parser.get_regression_type());
  a_model.set_regression_method_name(r_parser.get_regression_method_name());

  return;
}

void
sibpal_analysis::add_a_trait(regression_model& a_model,
                             const trait_type& trt)
{
  dependent_variable dv;

  dv.trait_index       = trt.trait_index;
  dv.binary            = trt.binary;
  dv.fixed_mean        = trt.fixed_mean;
  dv.fixed_male_mean   = trt.fixed_male_mean;
  dv.fixed_female_mean = trt.fixed_female_mean;
  dv.fixed_correlation = trt.fixed_correlation;

  a_model.set_trait(dv);

  return;
}

void
sibpal_analysis::add_a_marker(regression_model&  a_model,
                              const marker_type& mar)
{
  if( mar.effect == marker_type::DOMINANCE )
  {
    independent_variable param_add;
    param_add.type = independent_variable::MARKER;
    
    marker_type mar_add(mar.marker_index, marker_type::ADDITIVE, mar.x_linked);

    param_add.markers.push_back(mar_add);

    a_model.add_parameter(param_add);
  }

  independent_variable param;

  param.type = independent_variable::MARKER;
  param.markers.push_back(mar);

  a_model.add_parameter(param);

  return;
}

void
sibpal_analysis::add_a_covariate(regression_model&     a_model,
                                 const covariate_type& cov)
{
  independent_variable param;

  param.type = independent_variable::COVARIATE;
  param.covariates.push_back(cov);

  a_model.add_parameter(param);

  return;
}

bool
sibpal_analysis::check_main_effects(const regression_model&      a_model,
                                    const independent_variable&  interaction) const
{
  bool main_exist = true;

  for( size_t m = 0; m < interaction.markers.size(); ++m )
  {
    bool m_main_exist = false;

    m_main_exist = check_main_marker_effect(a_model, interaction.markers[m]);

    if( !m_main_exist )
      main_exist = false;

  }

  if( !main_exist )
    return main_exist;

  for( size_t c = 0; c < interaction.covariates.size(); ++c )
  {
    bool c_main_exist = false;

    c_main_exist = check_main_covariate_effect(a_model, interaction.covariates[c]);

    if( !c_main_exist )
      main_exist = false;
  }

  return main_exist;
}

bool
sibpal_analysis::check_main_marker_effect(const regression_model& a_model,
                                          const marker_type&      mar) const
{
  string m_name     = mar.name(*my_pairs);
  string m_eff_name = mar.effect_name();

  regression_model::parameter_const_iterator pi = a_model.parameter_begin();
  for( ; pi != a_model.parameter_end(); ++pi )
  {
    if(    pi->type == independent_variable::MARKER
        && pi->markers[0] == mar )
      return true;
  }

  return false;
}

bool
sibpal_analysis::check_main_covariate_effect(const regression_model& a_model,
                                             const covariate_type&   cov) const
{
  string c_name     = cov.name(*my_pairs);
  string c_eff_name = cov.effect_name();

  regression_model::parameter_const_iterator pi = a_model.parameter_begin();
  for( ; pi != a_model.parameter_end(); ++pi )
  {
    if(    pi->type == independent_variable::COVARIATE
        && pi->covariates[0] == cov )
      return true;
  }

  return false;
}

void
sibpal_analysis::set_output_options(const regression_parser& r_parser)
{
  size_t param_name_max = 11;
  size_t param_effect_name_max = 11;

  for( size_t i = 0; i < my_analysis_models.size(); ++i )
  {
    regression_model& a_model = my_analysis_models[i];

    for( size_t i = 0; i < a_model.get_parameter_count(); ++i )
    {
      string param_name     = a_model.get_parameter(i).name(*my_pairs);
      string param_eff_name = a_model.get_parameter(i).effect_name();

      param_name_max        = std::max(param_name.size(), param_name_max);
      param_effect_name_max = std::max(param_eff_name.size(), param_effect_name_max);
    }
  }

  for( size_t i = 0; i < my_analysis_models.size(); ++i )
  {
    regression_model& a_model = my_analysis_models[i];

    output_options output_opt = r_parser.get_output_options();;
    output_opt.param_name_max_size = param_name_max;
    output_opt.param_effect_name_max_size = param_effect_name_max;

    if( output_opt.dump_data )
    {
      bool dump_data = false;
      for( size_t i = 0; i < a_model.get_parameter_count(); ++i )
      {
        string param_name = a_model.get_parameter(i).name(*my_pairs);

        if( param_name == output_opt.dump_data_file )
        {
          dump_data = true;
          break;
        }
      }

      if( output_opt.dump_data_file == "TRAITS" )
        dump_data = true;

      output_opt.dump_data = dump_data;
    }

    if( output_opt.print_design_matrix )
    {
      bool print_design = false;
      for( size_t i = 0; i < a_model.get_parameter_count(); ++i )
      {
        string param_name = a_model.get_parameter(i).name(*my_pairs);

        if( param_name == output_opt.design_matrix_file )
        {
          print_design = true;
          break;
        }
      }

      if( output_opt.design_matrix_file == "TRAITS" )
        print_design = true;

      output_opt.print_design_matrix = print_design;
    }

    if( output_opt.print_correl_matrix )
    {
      bool print_correl = false;
      for( size_t i = 0; i < a_model.get_parameter_count(); ++i )
      {
        string param_name = a_model.get_parameter(i).name(*my_pairs);

        if( param_name == output_opt.correl_matrix_file )
        {
          print_correl = true;
          break;
        }
      }

      if( output_opt.correl_matrix_file == "TRAITS" )
        print_correl = true;

      output_opt.print_correl_matrix = print_correl;
    }

    a_model.set_output_options(output_opt);
  }

  return;
}

bool
sibpal_analysis::is_valid_model(const regression_model& a_model) const
{
  size_t m_count = 0;
  size_t c_count = 0;
  for( size_t c = 0; c < a_model.get_parameter_count(); ++c )
  {
    if( a_model.get_parameter(c).type == independent_variable::MARKER )
      ++m_count;

    if( a_model.get_parameter(c).type == independent_variable::COVARIATE )
      ++c_count;
  }

  if( a_model.get_regression_type() == ZERO_MARKER && (!c_count || m_count) )
    return false;
  //if( a_model.get_regression_type() == SINGLE_MARKER )
  //  if( m_count != 1 &&  )
  //  return false;
  if( a_model.get_regression_type() == MULTIPLE_MARKER && m_count < 2 )
    return false;

  return true;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
