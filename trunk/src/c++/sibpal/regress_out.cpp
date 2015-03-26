#include "sibpal/regress_out.h"

namespace SAGE   {
namespace SIBPAL {

regression_outfile::regression_outfile(ostream& sum, ostream& det, ostream& exp, cerrorstream& err)
                  : my_sum_output(sum), my_det_output(det), my_exp_output(exp), errors(err)
{}

void
regression_outfile::print_header(const TraitRegression& reg)
{
  my_analysis = &reg;

  my_detailed_out = reg.get_model().get_output_options().detailed_out;
  my_export_out   = reg.get_model().get_output_options().export_out;

  my_param_name_max_size        = get_analysis()->get_model().get_output_options().param_name_max_size;
  my_param_effect_name_max_size = get_analysis()->get_model().get_output_options().param_effect_name_max_size;

  my_x_type_count = get_analysis()->get_x_types().size();

  print_common_header(sum_output_stream());

  if( my_detailed_out )
    print_common_header(det_output_stream(), true);

  if( my_export_out )
    print_exp_header();

  return;
}

void
regression_outfile::print_results(size_t test_number)
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  print_sum_results(test_number);

  if( my_detailed_out )
    print_det_results(test_number);

  if( my_export_out )
    print_exp_results(test_number);

  return;
}

void
regression_outfile::print_footer()
{
  print_double_line(sum_output_stream());

  sum_output_stream() << endl << endl;

  if( my_detailed_out )
  {
    print_double_line(det_output_stream(), true);

    det_output_stream() << endl << endl;
  }

  return;
}

void
regression_outfile::print_instability_note()
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  out << "     - Beware there MAY be numerical instability in the above result." << endl;

  return;
}

//--------------------------------------------------------------------
//
// Header
//
//--------------------------------------------------------------------

void
regression_outfile::print_common_header(ostream& out, bool det)
{
  if( !out )
    return;

  print_double_line(out, det);
  print_analysis_title(out);
  print_double_line(out, det);

  out << endl;

  print_trait_name(out);
  print_dependent_name(out);
  print_subset_info(out);  
  print_other_option_info(out);
  print_legend(out);

  if( !det )
    print_double_line(out, det);

  return;
}

void
regression_outfile::print_exp_header()
{
  ostream& out = exp_output_stream();

  if(!out)
    return;

  out << "Trait"          << "\t"
      << "TestNo"         << "\t"
      << "IndependentVar" << "\t";

  if( get_analysis()->get_pairs().valid_distance_exist() )
    out << "Distance" << "\t";

  if( my_x_type_count > 1 )
    out << "Type" << "\t";

  out << "Pairs"     << "\t"
      << "Estimate"  << "\t"
      << "Std.Error" << "\t";

  if( get_analysis()->get_model().get_output_options().wide_out )
    out << "T-value" << "\t";

  out << "Nom.P-value";

  if( get_analysis()->get_model().get_pvalue_options().is_on )
  {
    out << "\t" << "Emp.P-value"
        << "\t" << "Replicates";
  }

  out << "\t" << "RSS";

  out << endl;

  return;
}

//-----------------------------------
// Common
//-----------------------------------

void
regression_outfile::print_analysis_title(ostream& out)
{
  if( !out )
    return;

  string pair_type = "Full Sibs";
  if( get_analysis()->get_use_pairs().first && get_analysis()->get_use_pairs().second )
    pair_type = "Full and Half Sibs";
  else if( !get_analysis()->get_use_pairs().first )
    pair_type = "Half Sibs";

  string reg_type = get_regression_type_to_string(get_analysis()->get_model().get_regression_type());
  
  string line = "  Haseman-Elston Regression Analysis of " + pair_type;

  out << line << endl << "  - " + reg_type << endl;

  return;
}

void
regression_outfile::print_trait_name(ostream& out)
{
  if( !out )
    return;

  const relative_pairs& pairs = get_analysis()->get_pairs();

  size_t t = get_analysis()->get_model().get_trait().trait_index;

  const RPED::RefTraitInfo &info = pairs.fped_info().trait_info(t);

  string line = "";

  if( info.type() == RPED::RefTraitInfo::continuous_trait )
  {
    line = "    Continuous trait : " + pairs.trait_name(t);
  }
  else if( info.type() == RPED::RefTraitInfo::binary_trait )
  {
    line =   "    Binary trait : " + pairs.trait_name(t)
           + ", affected = '" + info.string_affected_code()
           + "' , unaffected = '" + info.string_unaffected_code()
           + "'";
  }

  out << line << endl << endl;

  return;
}

void
regression_outfile::print_dependent_name(ostream& out)
{
  if( !out )
    return;

  string line = "    Dependent variate : ";

  if( get_analysis()->regression_method() )
  {
    line += get_analysis()->regression_method()->description();
  }
  else
    line += "(unknown regression method)";

  out << line << endl << endl;

  return;
}


void
regression_outfile::print_subset_info(ostream& out)
{
  if( !out )
    return;

  if( get_analysis()->get_model().get_subset_count() )
  {
    const relative_pairs& pairs = get_analysis()->get_pairs();

    string line = "    Subset used : ";

    for( size_t i = 0; i < get_analysis()->get_model().get_subset_count(); ++i )
    {
      if( i )
        line += ", ";
      line += get_analysis()->get_model().get_subset(i).name(pairs);
    }
    out << line << endl << endl;
  }

  return;
}

void
regression_outfile::print_other_option_info(ostream& out)
{
  if( !out )
    return;

  out << "    w1 : " << get_analysis()->get_model().get_analysis_options().w1 << endl << endl;


  out << "    Other options used :" << endl;

  if( get_analysis()->get_model().get_analysis_options().identity_weight )
    out << "      Identity weights = yes"
        << endl;
  else
    out << "      Identity weights = no"
        << endl;

  if( get_analysis()->get_model().get_analysis_options().robust_variance )
  {
    out << "      Robust variance  = yes"
        << endl;
    if( get_analysis()->get_model().get_analysis_options().leverage_adjustment )
      out << "          Leverage adjustment       = yes" << endl;
  }
  else
    out << "      Robust variance  = no"
        << endl;

  if( get_analysis()->get_model().get_analysis_options().sibship_mean )
    out << "      Use sibship mean = yes"
        //<< " (threshold="
        //<< get_analysis()->get_model().sibship_mean_threshold()
        //<< " sibs)"
        << endl;
  else
    out << "      Use sibship mean = no"
        << endl;

  if( get_analysis()->get_model().get_analysis_options().blup_mean )
    out << "      Use BLUP mean    = yes"
        << endl;
  else
    out << "      Use BLUP mean    = no"
        << endl;

  if( get_analysis()->get_model().get_data_options().skip_uninformative_pairs )
    out << "      Uninformative pairs skipped = yes"
        << endl;

  out << endl;

  return;
}

void
regression_outfile::print_legend(ostream& out)
{
  if( !out )
    return;

  out << "    Legend :" << endl
      << "      Note - kurtosis = coefficient of kurtosis - 3" << endl
      << "      *   - significance  .05 level;" << endl
      << "      **  - significance  .01 level;" << endl
      << "      *** - significance .001 level;" << endl
      << "      #   - Negative intercept set to 0";

  if( get_analysis()->get_use_pairs().first && get_analysis()->get_use_pairs().second )
    out << " or a common value is estimated;" << endl;
  else
    out << ";" << endl;

  if( get_analysis()->get_model().get_pvalue_options().is_on )
  {
    out << "      !!  - number of replicates reached the boundary value of " << endl
        << "            "
        << get_analysis()->get_model().get_pvalue_options().max_replicates
        << " replicates;" << endl;
  }

  out << endl;

  return;
}

void
regression_outfile::print_double_line(ostream& out, bool det)
{
  if( !out )
    return;

  string line = "==============================================================================";

  if( det )
  {
    out << line << "=======" << endl;
    return;
  }

  size_t pre_w = 64;

  if( my_param_name_max_size > 13 )
  {
    for( size_t i = 13; i < my_param_name_max_size; ++i )
    {
      line += "=";
      ++pre_w;
    }
  }

  if( my_x_type_count > 1 )
  {
    line += "=====";
    pre_w += 5;
  }

  if( my_param_effect_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_effect_name_max_size; ++i )
    {
      line += "=";
      ++pre_w;
    }
  }

  if( get_analysis()->get_model().get_output_options().wide_out )
  {
    line += "==========";
    pre_w += 10;
  }

  my_pre_width = pre_w;

//  if( print_adjusted_pvalues() )
//    out << "==================";
  
  if( get_analysis()->get_model().get_pvalue_options().is_on )
    line += "========================";

  out << line << endl;

  return;
}

//--------------------------------------------------------------------
//
// Results
//
//--------------------------------------------------------------------

//-----------------------------------
// Sum file
//-----------------------------------

void
regression_outfile::print_sum_results(size_t test_number)
{
  if( test_number == 1 )
  {
    print_sum_table_heading();
    print_sum_single_line();
  }

  if( my_x_type_count > 1 && test_number != 1 )
    print_sum_single_line();

  for( size_t i = 0; i < get_analysis()->get_reg_results().get_result_count(); ++i )
  {
    switch( get_analysis()->get_reg_results().get_result(i).type() )
    {
      case reg_result::none:
      case reg_result::intercept:
        break;

      case reg_result::param:
        print_sum_result(test_number, i);
        break;
    }
  }

  if(    get_analysis()->get_model().is_x_linked()
      && get_analysis()->get_reg_results().get_F_result().is_valid() )
  {
    print_joint_F_test();
  }

  return;
}

void
regression_outfile::print_sum_table_heading()
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  //Heading column count = 6+18+6+12+12+10+14 = 78
  // test = 6
  // independent + distance = 12 + 6 = 18
  // + if x_linked, type = 5
  // pairs = 6
  // parameter = 12
  // estimate = 12
  // std error = 10
  // nominal p-value = 10 + 4
  // + if wide_out, t-value = 10
  // + if emp p-value, p-value = 10 + 4, replicate = 10
  //Then, add extra for longer independent name

  string line = "Test  Independent       ";

  if( my_param_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_name_max_size; ++i )
      line += " ";
  }

  if( my_x_type_count > 1 )
    line += "     ";

  line +=  "                                        ";

  if( my_param_effect_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_effect_name_max_size; ++i )
      line += " ";
  }

  if( get_analysis()->get_model().get_output_options().wide_out )
    line += "          ";

  line += "  Nominal     ";
  if( get_analysis()->get_model().get_pvalue_options().is_on )
    line += "  Emp.        Number of ";

  out << line << endl;

  line = "No.   variable   ";

  if( my_param_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_name_max_size; ++i )
      line += " ";
  }

  if( get_analysis()->get_pairs().valid_distance_exist() )
    line += "   cM  ";
  else
    line += "       ";

  if( my_x_type_count > 1 )
    line += "Type ";

  line += "Pairs ";

  if( my_param_effect_name_max_size > 11 )
    for( size_t i = 11; i < my_param_effect_name_max_size; ++i )
      line += " ";

  line += "  Parameter    Estimate Std Error ";

  if( get_analysis()->get_model().get_output_options().wide_out )
    line += "  T-value ";

  line += "  P-value     ";

//  if( print_adjusted_pvalues() )
//    out << "     Adj. P-value ";

  if( get_analysis()->get_model().get_pvalue_options().is_on )
    line += "  P-value     Replicates";

  out << line << endl;

  return;
}

void
regression_outfile::print_sum_single_line()
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  string line = "----- -----------";

  if( my_param_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_name_max_size; ++i )
      line += "-";
  }

  if( get_analysis()->get_pairs().valid_distance_exist() )
    line += " ----- ";
  else
    line += "------ ";

  if( my_x_type_count > 1 )
    line += "---- ";

  line += "----- -----------";

  if( my_param_effect_name_max_size > 11 )
  {
    for( size_t i = 11; i < my_param_effect_name_max_size; ++i )
      line += "-";
  }

  line += " ----------- --------- ";
  if( get_analysis()->get_model().get_output_options().wide_out )
    line += "--------- ";

  line += "---------     ";

//  if( print_adjusted_pvalues() )
//    out << "     -------------";
  
  if( get_analysis()->get_model().get_pvalue_options().is_on )
    line += "---------     ----------";

  out << line << endl;

  return;
}

void
regression_outfile::print_sum_result(size_t test_number, size_t i)
{
  ostream& out = sum_output_stream();
  if(!out)
    return;

  const relative_pairs& pairs = get_analysis()->get_pairs();

  reg_result result = get_analysis()->get_reg_results().get_result(i);

  const independent_variable& param = get_analysis()->get_model().get_parameter(result.index());

  if( !param.valid )
  {
    cerr << "Internal error: invalid parameter in marker results." << endl;
    return;
  }

  string name = param.name(pairs);
  string effect_name = param.effect_name();

  if( i == 0 )
    out << left << setw(5) << test_number << " ";
  else
    out << setw(6) << " ";

  out << left << setw(my_param_name_max_size) << name << " ";

  if(    pairs.valid_distance_exist() 
      && effect_name.substr(0, 3) != "Cov" )
    out << fp(pairs.marker_distance(pairs.marker_find(name)), 5, 1) << " ";
  else
    out << setw(6) << " ";
    

  if( my_x_type_count > 1 )
  {
    if( result.get_pair_type() == MM )
      out << " BB  ";
    else if( result.get_pair_type() == MF )
      out << " BS  ";
    else if( result.get_pair_type() == FF )
      out << " SS  ";
  }

  if( i < my_x_type_count )
  {
    if( result.get_pair_type() == MM )
      out << right << setw(5) << get_analysis()->mm_pair_count() << " ";
    else if( result.get_pair_type() == MF )
      out << right << setw(5) << get_analysis()->mf_pair_count() << " ";
    else if( result.get_pair_type() == FF )
      out << right << setw(5) << get_analysis()->ff_pair_count() << " ";
    else
      out << right << setw(5) << get_analysis()->pair_count() << " ";
  }
  else
    out << setw(6) << " ";

  out << left << setw(my_param_effect_name_max_size) << effect_name << " ";

  out << fp(result.estimate.value(), 10);

  size_t i_cnt     = get_analysis()->get_reg_results().get_intercept_count();
  size_t i_cnt_int = get_analysis()->get_reg_intercept_results().get_intercept_count();

  if( i_cnt != i_cnt_int )
    out << " #";
  else
    out << "  ";
  out << fp(result.estimate.standard_error(), 9)   << " ";

  if( get_analysis()->get_model().get_output_options().wide_out )
    out << fp(result.estimate.tvalue(), 9)         << " ";

  size_t pw = 12;

  if( result.estimate.raw_pvalue() == QNAN )
    out << "QNAN           ";
  else
  {
    if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
      out << pval_scientific(result.estimate.raw_pvalue(), 10, 3, 3) << " ";
    else
      out << pval(result.estimate.raw_pvalue(), pw, 10, 3) << " ";
  }

//  if( print_adjusted_pvalues() )
//  {
//    if( result.estimate.adjusted_pvalue() == QNAN )
//      out << "QNAN           ";
//    else
//      out << pval(result.estimate.adjusted_pvalue(), pw, 10, 3);
//  }
  
  // Print only when the term involves a marker (since that is all we permute)
  if(    get_analysis()->get_model().get_pvalue_options().is_on
      && param.markers.size() )
  {
    if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
      out << pval_scientific(result.estimate.empirical_pvalue(), 10, 3, 3) << " ";
    else
      out << pval(result.estimate.empirical_pvalue(), pw, 10, 3) << " ";

    out << "(" << setw(8) 
        << result.estimate.total_replicates()
        << ")";
  }

  out << endl;

  return;
}

void
regression_outfile::print_joint_F_test()
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  const F_result_type& f_re = get_analysis()->get_reg_results().get_F_result();

  size_t df1 = f_re.get_df1();
  size_t df2 = f_re.get_df2();

  double f_stat = f_re.get_F_statistic();
  double f_pval = f_re.get_F_pvalue();

#if 0
  cout << "df1 = " << df1 << endl
       << "df2 = " << df2 << endl
       << "f_stat = " << f_stat << endl
       << "f_pval = " << f_pval << endl;
#endif

  stringstream s;

  s << "      F-statistic : "
    << fp(f_stat, 6, 4)
    << " on " << df1 << " and " << df2 << " degrees of freedom ";

  out << left << setw(my_pre_width) << " ";
  out << "---------";

//  if( print_adjusted_pvalues() )
//    out << "                  ";
  
  if( get_analysis()->get_model().get_pvalue_options().is_on )
    out << "     ---------     ----------";
  out << endl;

  out << left << setw(my_pre_width) << s.str();

  if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
    out << pval_scientific(f_pval, 10, 3, 3) << " ";
  else
    out << pval(f_pval, 12, 10, 3) << " ";

  if( get_analysis()->get_model().get_pvalue_options().is_on )
  {
    if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
      out << pval_scientific(f_re.get_F_empirical_pvalue(), 10, 3, 3) << " ";
    else
      out << pval(f_re.get_F_empirical_pvalue(), 12, 10, 3) << " ";

    out << "(" << setw(8)
        << f_re.get_total_replicates()
        << ")";
  }

  out << endl;

  return;
}

//-----------------------------------
// Det file
//-----------------------------------

void
regression_outfile::print_det_results(size_t test_number)
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  print_double_line(out, true);
  out << "Test " << test_number << endl;
  print_double_line(out, true);

  print_model_info();
  print_sample_info();
  print_trait_info();
  print_dependent_info(test_number);
  print_regression_info();
  
  return;
}

void
regression_outfile::print_model_info()
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "-----" << endl
      << "Model" << endl
      << "-----" << endl
      << endl;

  const relative_pairs& pairs = get_analysis()->get_pairs();

  size_t t = get_analysis()->get_model().get_trait().trait_index;

  string line = pairs.trait_name(t) + " ~ Intercept";

  for( size_t i = 0; i < get_analysis()->get_model().get_parameter_count(); ++i )
  {
    const independent_variable& param = get_analysis()->get_model().get_parameter(i);
    
    string p_name = " + " + param.name(pairs);

    if( param.type == independent_variable::MARKER )
      p_name += "(" + param.short_effect_name() + ")";

    line += p_name;
  }

  line += " + e";

  out << line << endl << endl;

  return;
}

void
regression_outfile::print_sample_info()
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "------" << endl
      << "Sample" << endl
      << "------" << endl
      << endl;

  out <<  "Number of all sibs            = "
      << right << setw(4) << get_analysis()->get_model().get_trait().trait_all_sibs_info.count()
      << endl;

  if( get_analysis()->get_model().get_trait().trait_all_sibs_info.count()
      != get_analysis()->get_model().get_trait().trait_used_sibs_info.count() )
  {
    out <<  "Number of sibs used           = "
        << right << setw(4) << get_analysis()->get_model().get_trait().trait_used_sibs_info.count()
        << endl;
  }
  if( get_analysis()->get_use_pairs().first )
  {
    out << "Number of full sib pairs      = "
        << right << setw(4) << get_analysis()->fsib_pair_count() << endl;

    if( get_analysis()->get_model().is_x_linked() )
    {
      if( get_analysis()->get_use_x_pairs()[MM] && get_analysis()->mm_pair_count() )
        out << "  brother-brother pairs       = "
            << right << setw(4) << get_analysis()->mm_pair_count() << endl;

      if( get_analysis()->get_use_x_pairs()[MF] && get_analysis()->mf_pair_count() )
        out << "  brother-sister pairs        = "
            << right << setw(4) << get_analysis()->mf_pair_count() << endl;

      if( get_analysis()->get_use_x_pairs()[FF] && get_analysis()->ff_pair_count() )
        out << "  sister-sister pairs         = "
            << right << setw(4) << get_analysis()->ff_pair_count() << endl;
    }
  }
  if( get_analysis()->get_use_pairs().second )
  {
    out << "Number of half sib pairs      = "
        << right << setw(4) << get_analysis()->hsib_pair_count() << endl;
  }

  out << endl;

  return;
}

void
regression_outfile::print_trait_info()
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "-----" << endl
      << "Trait" << endl
      << "-----" << endl
      << endl;

  const dependent_variable& tparam = get_analysis()->get_model().get_trait();

  if(finite(tparam.fixed_mean))
    out << "Population mean               = "
        << fp(tparam.fixed_mean, 9) << endl;

  if(finite(tparam.trait_all_sibs_info.mean()))
    out << "Sample mean                   = "
        << fp(tparam.trait_all_sibs_info.mean(), 9) << endl;

  if(finite(tparam.trait_all_sibs_info.variance()))
    out << "Sample variance               = "
        << fp(tparam.trait_all_sibs_info.variance(), 9) << endl;

  if(finite(tparam.trait_all_sibs_info.skewness()))
    out << "Sample skewness               = "
        << fp(tparam.trait_all_sibs_info.skewness(), 9) << endl;

  if(finite(tparam.trait_all_sibs_info.kurtosis()))
    out << "Sample kurtosis               = "
        << fp(tparam.trait_all_sibs_info.kurtosis()-3.0, 9) << endl;

  if( get_analysis()->get_use_pairs().first )
  {
    out << "Pairwise full sib correlation = ";

    double rho = tparam.fisher_fsib_correlation;
    if(finite(rho))
      out << fp(rho,9) << endl;
    else
      out << "   ------" << endl;

    out << "Intra sibship correlation     = ";

    rho = tparam.correlation;
    if(finite(rho))
      out << fp(rho,9) << endl;
    else
      out << "   ------" << endl;
  }
  if( get_analysis()->get_use_pairs().second )
  {
    out << "Pairwise half sib correlation = ";

    double rho = tparam.fisher_hsib_correlation;
    if(finite(rho))
      out << fp(rho,9) << endl;
    else
      out << "   ------" << endl;
  }

  out << endl;

  return;
}

void
regression_outfile::print_dependent_info(size_t test_number)
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "-----------------" << endl
      << "Dependent variate" << endl
      << "-----------------" << endl
      << endl;

  if( get_analysis()->get_use_pairs().first )
  {
    if( get_analysis()->get_use_pairs().second )
      out << "- Full sib pairs" << endl;
  
    out << "Correlation between pairs with no sibs in common = ";
    double rho = get_analysis()->get_model().get_trait().p_fsib_correlation[0];
    if(finite(rho))
      out << fp(rho, 9) << endl;
    else
      out << "---------" << endl;

    out << "Correlation between pairs with one sib in common = ";
    rho = get_analysis()->get_model().get_trait().p_fsib_correlation[1];
    if(finite(rho))
      out << fp(rho,9) << endl;
    else
      out << "---------" << endl;

    if( toUpper(get_analysis()->get_model().get_regression_method_name()) == "W2" )
    {
      const map<size_t, pair<double, double> >& ws = get_analysis()->get_reg_results().get_full_w();
      double w1 = ws.begin()->second.first;
      double w2 = ws.begin()->second.second;

      double  w = 0.;
      double _w = 1.;

      if( fabs(w1) > 1.0e-10 && fabs(w2) > 1.0e-10 )
      {
        w = w2 / (w1+w2);
        _w = w1 / (w1+w2);
      }

      out << "w = " << w  << ", 1-w = " << _w << endl;
    }
    else if( toUpper(get_analysis()->get_model().get_regression_method_name().substr(0,1)) == "W" )
    {
      const map<size_t, pair<double, double> >& ws = get_analysis()->get_reg_results().get_full_w();
      map<size_t, pair<double, double> >::const_iterator wi = ws.begin();

      for( ; wi != ws.end(); ++wi )
      {
        double w1 = wi->second.first;
        double w2 = wi->second.second;

        double  w = 0.;
        double _w = 1.;

        if( fabs(w1) > 1.0e-10 && fabs(w2) > 1.0e-10 )
        {
          w = w2 / (w1+w2);
          _w = w1 / (w1+w2);
        }

        out << "sibship size " << wi->first
            << " : w = " << w
            << ", 1-w = " << _w << endl;
      }
    }
  }
  if( get_analysis()->get_use_pairs().second )
  {
    if( get_analysis()->get_use_pairs().first )
      out << endl
          << "- Half sib pairs" << endl;

    out << "Correlation between pairs with no sibs in common = ";
    double rho = get_analysis()->get_model().get_trait().p_hsib_correlation[0];
    if(finite(rho))
      out << fp(rho, 9) << endl;
    else
      out << "---------" << endl;

    out << "Correlation between pairs with one sib in common = ";
    rho = get_analysis()->get_model().get_trait().p_hsib_correlation[1];
    if(finite(rho))
      out << fp(rho,9) << endl;
    else
      out << "---------" << endl;

    if( toUpper(get_analysis()->get_model().get_regression_method_name()) == "W2" )
    {
      const map<size_t, pair<double, double> >& ws = get_analysis()->get_reg_results().get_half_w();
      double w1 = ws.begin()->second.first;
      double w2 = ws.begin()->second.second;

      double  w = 0.;
      double _w = 1.;

      if( fabs(w1) > 1.0e-10 && fabs(w2) > 1.0e-10 )
      {
        w = w2 / (w1+w2);
        _w = w1 / (w1+w2);
      }

      out << "w = " << w << ", 1-w = " << _w << endl;
    }
    else if( toUpper(get_analysis()->get_model().get_regression_method_name().substr(0,1)) == "W" )
    {
      const map<size_t, pair<double, double> >& ws = get_analysis()->get_reg_results().get_half_w();
      map<size_t, pair<double, double> >::const_iterator wi = ws.begin();

      for( ; wi != ws.end(); ++wi )
      {
        double w1 = wi->second.first;
        double w2 = wi->second.second;

        double  w = 0.;
        double _w = 1.;

        if( fabs(w1) > 1.0e-10 && fabs(w2) > 1.0e-10 )
        {
          w = w2 / (w1+w2);
          _w = w1 / (w1+w2);
        }

        out << "sibship size " << wi->first
            << " : w = " << w
            << ", 1-w = " << _w << endl;
      }
    }
  }

  out << endl;
  out << "Correlation between squared difference" << endl
      << "                and squared mean corrected sum   = "
      << fp(get_analysis()->get_model().get_trait().sum_diff_correlation, 9)
      << endl << endl;

  if( toUpper(get_analysis()->get_model().get_regression_method_name().substr(0,1)) == "W" )
  {
    if( get_analysis()->get_reg_intercept_results().get_sum_residual_variance() == 0.0 )
    {
      out << "--Note that residual variance of SUM regression = 0,"
          << endl
          << "  therefore the dependent variate is the same as for DIFF."
          << endl << endl;

      errors << priority(information)
                << "In Test " << test_number
                << ", note that residual variance of SUM regression = 0, "
                << "therefore the dependent variate is the same as for DIFF."
                << endl;
    }
    else if( get_analysis()->get_reg_intercept_results().get_diff_residual_variance() == 0.0 )
    {
      out << "--Note that residual variance of DIFF regression = 0,"
          << endl
          << "  therefore the dependent variate is for SUM."
          << endl << endl;

      errors << priority(information)
             << "In Test " << test_number
             << ", note that residual variance of DIFF regression = 0, "
             << "therefore the dependent variate is the same as for SUM."
             << endl;
    }
  }

  return;
}

void
regression_outfile::print_regression_info()
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "----------" << endl
      << "Regression" << endl
      << "----------" << endl
      << endl;

  size_t pw = 12;

  const relative_pairs& pairs = get_analysis()->get_pairs();

  size_t i_cnt_no_int = get_analysis()->get_reg_results().get_intercept_count();
  size_t i_cnt_int    = get_analysis()->get_reg_intercept_results().get_intercept_count();

  size_t r_cnt_int    = get_analysis()->get_reg_intercept_results().get_result_count();
#if 0
  size_t r_cnt_no_int = get_analysis()->get_reg_results().get_result_count();
  cout << "r_cnt_no_int = " << r_cnt_no_int << ", r_cnt_int = " << r_cnt_int << endl;
  cout << "i_cnt_no_int = " << i_cnt_no_int << ", i_cnt_int = " << i_cnt_int << endl;
#endif
  if( i_cnt_no_int != i_cnt_int )
  {
    out << "                  before intercept is set           # after intercept is set"
        << endl
        << "                  -----------------------------     -----------------------------"
        << endl;
  }

  out << "                   Estimate Std Error   P-value     ";
  if( i_cnt_no_int != i_cnt_int )
    out << " Estimate Std Error   P-value     ";
  out << endl;
  
  out << "----------------- --------- --------- ---------     ";
  if( i_cnt_no_int != i_cnt_int )
    out << "--------- --------- ---------     ";
  out << endl;

  // Print intercept
  for( size_t i = r_cnt_int-i_cnt_int, ii = r_cnt_int-i_cnt_int; i < r_cnt_int; ++i )
  {
    const reg_result& result = get_analysis()->get_reg_intercept_results().get_result(i);

    string name = "Intercept";
    double new_int = 0.0;

    if( my_x_type_count > 1 )
    {
      if( result.get_pair_type() == MM )
        name += "(BB)     ";
      else if( result.get_pair_type() == MF )
        name += "(BS)     ";
      else if( result.get_pair_type() == FF )
        name += "(SS)     ";

      if(    i_cnt_no_int != i_cnt_int
          && i_cnt_no_int
          && result.get_pair_type() == get_analysis()->get_reg_results().get_result(ii).get_pair_type() )
      {
        new_int = get_analysis()->get_reg_results().get_result(ii).estimate.value();
        ++ii;
      }
    }
    else
    {
      if( i_cnt_int > 1 && i == r_cnt_int-i_cnt_int )
        name += "(Full)   ";
      else if( i_cnt_int > 1 )
        name += "(Half)   ";
      else
        name += "         ";

      if(    i_cnt_no_int != i_cnt_int
          && i_cnt_no_int )
      {
        new_int = get_analysis()->get_reg_results().get_result(ii).estimate.value();

        if( name == "Intercept(Half)   " && result.estimate.value() < 0.0 )
          new_int = 0.0;
      }
    }

    out << name << fp(result.estimate.value(), 9) << " ";

    if( i_cnt_no_int != i_cnt_int )
    {
      out << "                        " << fp(new_int, 9);
    }

    out << endl;
  }

  // Print slope & covariate coeff.
  for( size_t i = 0; i < r_cnt_int-i_cnt_int; ++i )
  {
    const reg_result& result = get_analysis()->get_reg_intercept_results().get_result(i);
    const independent_variable& param = get_analysis()->get_model().get_parameter(result.index());

    string name = param.name(pairs);
    if(    param.type == independent_variable::MARKER
        && param.short_effect_name() != "A+D" )
      name += "(" + param.short_effect_name() + ")";

    if( my_x_type_count > 1 )
    {
      if( result.get_pair_type() == MM )
        name += "(BB)";
      else if( result.get_pair_type() == MF )
        name += "(BS)";
      else if( result.get_pair_type() == FF )
        name += "(SS)";
    }

    out << left << setw(17) << name << " ";
    out << fp(result.estimate.value(), 9) << " ";
    out << fp(result.estimate.standard_error(), 9) << " ";

    if( result.estimate.raw_pvalue() == QNAN )
      out << "QNAN           ";
    else
    {
      if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
        out << pval_scientific(result.estimate.raw_pvalue(), 10, 3, 3) << " ";
      else
        out << pval(result.estimate.raw_pvalue(), pw, 10, 3) << " ";
    }

    if( i_cnt_no_int != i_cnt_int )
    {
      const reg_result& re = get_analysis()->get_reg_results().get_result(i);

      out << fp(re.estimate.value(), 9) << " ";
      out << fp(re.estimate.standard_error(), 9) << " ";

      if( re.estimate.raw_pvalue() == QNAN )
        out << "QNAN           ";
      else
      {
        if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
          out << pval_scientific(re.estimate.raw_pvalue(), 10, 3, 3) << " ";
        else
          out << pval(re.estimate.raw_pvalue(), pw, 10, 3) << " ";
      }
    }

    out << endl;
  }

  out << endl;

  out << "Total variance    "
      << fp(get_analysis()->get_reg_intercept_results().get_total_variance(),9);
  if( i_cnt_no_int != i_cnt_int )
    out << "                         "
        << fp(get_analysis()->get_reg_results().get_total_variance(),9);
  out << endl;

  out << "Residual variance "
      << fp(get_analysis()->get_reg_intercept_results().get_residual_variance(),9);
  if( i_cnt_no_int != i_cnt_int )
    out << "                         "
        << fp(get_analysis()->get_reg_results().get_residual_variance(),9);
  out << endl;

  out << "Residual skewness "
      << fp(get_analysis()->get_reg_intercept_results().get_residual_info().skewness(),9);
  if( i_cnt_no_int != i_cnt_int )
    out << "                         "
        << fp(get_analysis()->get_reg_results().get_residual_info().skewness(),9);
  out << endl;

  out << "Residual kurtosis "
      << fp(get_analysis()->get_reg_intercept_results().get_residual_info().kurtosis()-3,9);
  if( i_cnt_no_int != i_cnt_int )
    out << "                         "
        << fp(get_analysis()->get_reg_results().get_residual_info().kurtosis()-3,9);
  out << endl;

  return;
}

//-----------------------------------
// Exp file
//-----------------------------------

void
regression_outfile::print_exp_results(size_t test_number)
{
  for( size_t i = 0; i < get_analysis()->get_reg_results().get_result_count(); ++i )
  {
    switch( get_analysis()->get_reg_results().get_result(i).type() )
    {
      case reg_result::none:
      case reg_result::intercept:
        break;

      case reg_result::param:
        print_exp_result(test_number, i);
        break;
    }
  }

  return;
}

void
regression_outfile::print_exp_result(size_t test_number, size_t i)
{
  ostream& out = exp_output_stream();
  if(!out)
    return;

  const relative_pairs& pairs = get_analysis()->get_pairs();

  reg_result result = get_analysis()->get_reg_results().get_result(i);

  const independent_variable& param = get_analysis()->get_model().get_parameter(result.index());

  if( !param.valid )
  {
    cerr << "Internal error: invalid parameter in marker results." << endl;
    return;
  }

  size_t t = get_analysis()->get_model().get_trait().trait_index;

  string tname = pairs.trait_name(t);

  string name      = param.name(pairs);
  string effect_name = param.effect_name();

  out << tname << "\t"
      << test_number << "\t"
      << name << "\t";

  if( pairs.valid_distance_exist() )
    out << fp(pairs.marker_distance(pairs.marker_find(name)), 5, 1) << "\t";
  else
    out << "--------" << "\t";
    
  if( my_x_type_count > 1 )
  {
    if( result.get_pair_type() == MM )
      out << "BB\t";
    else if( result.get_pair_type() == MF )
      out << "BS\t";
    else if( result.get_pair_type() == FF )
      out << "SS\t";
  }

  if( my_x_type_count > 1 )
  {
    if( result.get_pair_type() == MM )
      out << get_analysis()->mm_pair_count() << "\t";
    else if( result.get_pair_type() == MF )
      out << get_analysis()->mf_pair_count() << "\t";
    else if( result.get_pair_type() == FF )
      out << get_analysis()->ff_pair_count() << "\t";
  }
  else
    out << get_analysis()->pair_count() << "\t";

  out << result.estimate.value() << "\t";
  out << result.estimate.standard_error() << "\t";

  if( get_analysis()->get_model().get_output_options().wide_out )
    out << result.estimate.tvalue() << "\t";

  if( result.estimate.raw_pvalue() == QNAN )
    out << "QNAN\t";
  else
  {
    if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
      out << doub2str(result.estimate.raw_pvalue(), 10, 3, ios::scientific | ios::showpoint);
    else
      out << result.estimate.raw_pvalue();
  }

  // Print only when the term involves a marker (since that is all we permute)
  if(    get_analysis()->get_model().get_pvalue_options().is_on
      && param.markers.size() )
  {
    if( get_analysis()->get_model().get_output_options().pvalues_scientific_notation )
      out << "\t" << doub2str(result.estimate.empirical_pvalue(), 10, 3, ios::scientific | ios::showpoint);
    else
      out << "\t" << result.estimate.empirical_pvalue();

    out << "\t" << result.estimate.total_replicates();
  }

  out << "\t" << get_analysis()->get_reg_results().get_residual_square_info().sum();

  out << endl;

  return;
}

} // end of namespace SIBPAL
} // end of namespace SAGE
