//=============================================================================
// File:    output.cpp
//
// Author:  Yeunjoo Song
//
// History: Initial implementation
//
// Notes:   This file contains implementation for relpal_textfile.
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/output.h"

namespace SAGE   {
namespace RELPAL {

//--------------------------------------------------------------------
//
// relpal_outfile
//
//--------------------------------------------------------------------

relpal_outfile::relpal_outfile(ostream& sum, ostream& det, ostream& exp, const relpal_analysis& ra)
              : my_analysis(ra),
                detailed_out(my_analysis.get_output_options().detailed_out),
                export_out(my_analysis.get_output_options().export_out),
                my_sum_output(sum), my_det_output(det), my_exp_output(exp)
{
  my_var_opt.resize(0);
}

void
relpal_outfile::print_header(bool first_level)
{
  print_sum_header(first_level);

  if( detailed_out )
    print_det_header(first_level);

  if( export_out )
    print_exp_header();

  return;
}

void
relpal_outfile::print_results(size_t test_number, bool first_level, bool valid, bool reliable)
{
  if( first_level )
    print_first_level_results(test_number, valid);
  else
    print_second_level_results(test_number, valid, reliable);

  return;
}

void
relpal_outfile::print_footer(bool first_level)
{
  ostream& out = sum_output_stream();

  print_double_line(out);  

  out << endl << endl;

  if( detailed_out )
    det_output_stream() << endl;

  return;
}

//--------------------------------------------------------------------
//
// Header
//
//--------------------------------------------------------------------

void
relpal_outfile::print_sum_header(bool first_level)
{
  ostream& out = sum_output_stream();
  if( !out )
    return;

  print_double_line(out);
  print_analysis_title(out, first_level, "Summary");
  print_double_line(out);

  out << endl;

  print_trait_summary(out);
  print_subset_summary(out);
  print_legend(out);

  if( first_level )
  {
    print_double_line(out);
    print_reg_table_heading(out, true, false);
    print_reg_single_line(out, true, false);
  }
  else
  {
    print_pval_option(out);

    print_double_line(out);
    print_score_table_heading(out, true);
    print_score_single_line(out);
  }

  return;
}

void
relpal_outfile::print_det_header(bool first_level)
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  print_double_line(out, false);
  print_analysis_title(out, first_level, "Detailed");
  print_double_line(out, false);

  out << endl;

  print_trait_summary(out);
  print_subset_summary(out);
  print_legend(out);

  if( !first_level )
    print_pval_option(out);

  print_double_line(out, false);

  return;
}

void
relpal_outfile::print_exp_header()
{
  ostream& out = exp_output_stream();
  if( !out )
    return;

  out << "Variable,cM,P-value" << endl;

  return;
}

//--------------------------------------------------------------------

void
relpal_outfile::print_analysis_title(ostream& out, bool first_level, string out_type)
{
  if( !out )
    return;

  stringstream s;

  s << "  Two-Level Haseman-Elston Regression Analysis for General Pedigree "
    << endl;
    //<< "   - " << get_regression_type_to_string(get_analysis().get_regression_type())
    //<< " regression";

  s << endl;

  if( first_level )
    s << "  - First (Individual) Level Wald Tests ";
  else
    s << "  - Second (Pedigree) Level Score Tests ";

  s << out_type << " Output" << endl;

  out << s.str();

  return;
}

void
relpal_outfile::print_trait_summary(ostream& out)
{
  if( !out )
    return;

  const relative_pairs&   pairs = *get_analysis().get_pairs();
  const regression_model& model = get_analysis().get_current_model();

  if( model.get_trait_count() > 1 )
    out << "  Traits : ";
  else
    out << "  Trait :  ";

  for( size_t t = 0; t < model.get_trait_count(); ++t )
  {
    if( t )
      out << "           ";

    size_t ti = model.get_trait(t).trait_index;
    const RPED::RefTraitInfo &info = pairs.fped_info().trait_info(ti);

    if( model.get_trait_count() > 1 )
      out << "(" << t+1 << ") - ";

    out << pairs.trait_name(ti);

    if( info.type() == RPED::RefTraitInfo::continuous_trait )
    {
      out << ", Quantitative"
          << endl;
    }
    else if( info.type() == RPED::RefTraitInfo::binary_trait )
    {
      out << ", Binary"
          << ", affected = '" << info.string_affected_code()
          << "' , unaffected = '" << info.string_unaffected_code()
          << "'" << endl;
    }
  }
  out << endl;

  return;
}

void
relpal_outfile::print_subset_summary(ostream& out)
{
  if( !out )
    return;

  const data_options& data  = get_analysis().get_data_options();

  if( !data.subsets.size() )
    return;

  out << "  Subset used : ";

  for( size_t i = 0; i < data.subsets.size(); ++i )
  {
    if( i )
      out << ", ";
    out << data.subsets[i].name();
  }
  out << endl << endl;

  return;
}

void
relpal_outfile::print_pval_option(ostream& out)
{
  if( !out )
    return;

  const regression_model& model = get_analysis().get_current_model();

  if( model.get_analysis_options().naive_variance )
  {
    my_var_opt.push_back("nai");

    out << "    nai - naive variance" << endl;
  }

  if( model.get_analysis_options().sandwich_variance )
  {
    my_var_opt.push_back("sdw");

    out << "    sdw - robust sandwich variance" << endl;
  }

  if( model.get_analysis_options().alternative_variance )
  {
    my_var_opt.push_back("alt");

    out << "    alt - alternative variance" << endl;
  }
  if( model.get_analysis_options().IBD_variance )
  {
    my_var_opt.push_back("ibd");

    out << "    ibd - allele sharing variance" << endl;
  }
  out << endl;

  out << "    #   - The number of pedigrees is too small for this result to be" << endl
      << "          reliable when analyzing this number of traits." << endl
      << endl;

  out << "  Empirical p-value options used :" << endl
      << "    seed = " << model.get_pvalue_options().seed << endl
      << "    min replicates = " << model.get_pvalue_options().min_replicates << endl
      << "    max replicates = " << model.get_pvalue_options().max_replicates << endl;

  out << endl;

  return;
}

void
relpal_outfile::print_legend(ostream& out)
{
  if( !out )
    return;

  out << "  Legend :" << endl
      << "    *   - significance  .05 level" << endl
      << "    **  - significance  .01 level" << endl
      << "    *** - significance .001 level" << endl;

  out << endl;

  return;
}

void
relpal_outfile::print_table_heading_line1(ostream& out, bool print_cM, bool print_count)
{
  if( !out )
    return;

  const relative_pairs& pairs = *get_analysis().get_pairs();
  size_t max_name_size = get_analysis().get_output_options().param_name_max;

  out << left << setw(max_name_size) << "Test" ;

  if( print_cM && pairs.valid_distance_exist() )
    out << "      ";

  if( print_count )
    out << "      ";

  return;
}

void
relpal_outfile::print_table_heading_line2(ostream& out, bool sum, bool print_cM, bool print_count)
{
  if( !out )
    return;

  const relative_pairs& pairs = *get_analysis().get_pairs();
  size_t max_name_size = get_analysis().get_output_options().param_name_max;

  if( sum )
    out << left << setw(max_name_size) << "No   Variable";
  else
    out << left << setw(max_name_size) << "Variable";

  if( print_cM && pairs.valid_distance_exist() )
    out << "    cM";

  if( print_count )
    out << " Count";

  return;
}

void
relpal_outfile::print_double_line(ostream& out, bool sum)
{
  if( !out )
    return;

  //const relative_pairs& pairs = *get_analysis().get_pairs();
  size_t max_name_size = get_analysis().get_output_options().param_name_max;

  stringstream s;

  s << "=========================";

  if( max_name_size > DEFAULT_PARAM_NAME_MAX )
  {
    for( size_t i = DEFAULT_PARAM_NAME_MAX; i < max_name_size; ++i )
    {
      s << "=";
    }
  }

  s << "===========================================================";

  if( !sum )
  {
    size_t t_count = get_analysis().get_current_model().get_trait_count();

    for( size_t t = 2; t < t_count; ++t )
      s << "==========";
  }

  out << s.str() << endl;

  return;
}

void
relpal_outfile::print_single_line(ostream& out, bool print_cM, bool print_count)
{
  if( !out )
    return;

  const relative_pairs& pairs = *get_analysis().get_pairs();
  size_t max_name_size = get_analysis().get_output_options().param_name_max;

  stringstream s;

  s << "-------------------------";

  if( max_name_size > DEFAULT_PARAM_NAME_MAX )
  {
    for( size_t i = DEFAULT_PARAM_NAME_MAX; i < max_name_size; ++i )
      s << "-";
  }

  if( print_cM && pairs.valid_distance_exist() == true )
    s << " -----";

  if( print_count )
    s << " -----";

  out << s.str();

  return;
}

void
relpal_outfile::print_reg_table_heading(ostream& out, bool sum, bool print_cM, bool print_count)
{
  if( !out )
    return;

  size_t t_count = get_analysis().get_current_model().get_trait_count();

  print_table_heading_line1(out, print_cM, print_count);
  //out << "                                    Nominal    ";

  if( sum )
    out << "             Nominal    ";
  else
  {
    out << "             ";
    for( size_t t = 0; t < t_count; ++t )
      out << "          ";

    out << "             Nominal    ";
  }

  out << endl;

  print_table_heading_line2(out, sum, print_cM, print_count);

  if( sum )
    out << "   Chi-sq.   P-value    ";
  else
  {
    //out << "     Estimate Std Error   Z-value   P-value    ";

    out << "     Estimate";
    if( t_count == 1 )
      out << "  Variance";
    else
    {
      out << " Variance-Covariance";
      for( size_t t = 2; t < t_count; ++t )
        out << "          ";
    }

    out << "   Chi-sq.   P-value    ";
  }

  out << endl;

  return;
}

void
relpal_outfile::print_score_table_heading(ostream& out, bool sum, bool print_cM, bool print_count)
{
  if( !out )
    return;

  print_table_heading_line1(out, print_cM, print_count);
  out << "     Unadjusted   Adjusted Empirical     Number of ";

  out << endl;

  print_table_heading_line2(out, sum, print_cM, print_count);
  out << " Var    T-value    T-value   P-value     Replicates";

  out << endl;

  return;
}

void
relpal_outfile::print_reg_single_line(ostream& out, bool sum, bool print_cM, bool print_count)
{
  if( !out )
    return;

  size_t t_count = get_analysis().get_current_model().get_trait_count();

  print_single_line(out, print_cM, print_count);
  //out << " ------------ --------- --------- ---------    ";

  if( sum )
    out << " --------- ---------    ";
  else
  {
    out << " ------------";

    if( t_count == 1 )
      out << " ---------";
    else
    {
      out << " -------------------";
      for( size_t t = 2; t < t_count; ++t )
        out << "----------";
    }

    out << " --------- ---------    ";
  }

  out << endl;

  return;
}

void
relpal_outfile::print_score_single_line(ostream& out, bool print_cM, bool print_count)
{
  if( !out )
    return;

  print_single_line(out, print_cM, print_count);
  out << " --- ---------- ---------- ---------     ----------";

  out << endl;

  return;
}

//--------------------------------------------------------------------
//
// Results
//
//--------------------------------------------------------------------

void
relpal_outfile::print_first_level_results(size_t test_number, bool valid)
{
  const relative_pairs&   pairs = *get_analysis().get_pairs();
  const regression_model& model = get_analysis().get_current_model();

  size_t max_name_size = get_analysis().get_output_options().param_name_max;
  string mem_count_str = long2str(get_analysis().get_member_count(), 6);

  stringstream s;
  s << left << setw(5) << long2str(test_number, 1);

  string t_str = s.str();

  if( detailed_out )
  {
    print_det_first_level_results(test_number);
  }

  size_t t_count     = model.get_trait_count();
  bool   t_printed   = false;
  string last_name   = "";
  double last_pvalue = QNAN;

  for( size_t i = 0; i < model.get_ind_parameter_count(); )
  {
    const independent_variable& param = model.get_ind_parameter(i);

    i += t_count;

    if( param.type == independent_variable::INTERCEPT )
      continue;

    if( t_printed )
      t_str = "     ";

    string name = param.name(pairs);

    if( valid )
    {
      matrix beta, var;

      get_reg_param_results(i-t_count, t_count,
                            get_analysis().get_current_result().H0_result.ind_beta,
                            get_analysis().get_current_result().H0_result.ind_variance,
                            beta, var);

      pair<double, double> T_P = compute_multivariate_T_value(beta, var);

      print_reg_result_line(sum_output_stream(), t_str+name, max_name_size, "", mem_count_str, T_P.first, T_P.second);

      if( detailed_out )
        print_reg_result_block(det_output_stream(), name, max_name_size, "", "", beta, var, T_P.first, T_P.second);

      last_pvalue = T_P.second;
    }
    else
    {
      print_fail_line(sum_output_stream(), t_str+name, max_name_size, "", mem_count_str);

      if( detailed_out )
        print_fail_line(det_output_stream(), name, max_name_size, "", "");
    }

    t_printed = true;
    last_name = name;
  }

  if( detailed_out )
  {
    if( t_count == 1 || !valid )
      det_output_stream() << endl;
    print_double_line(det_output_stream(), false);
  }

  if( export_out )
    if( valid )
      exp_output_stream() << last_name << ",," << last_pvalue << endl;
    else
      exp_output_stream() << last_name << ",,FAIL" << endl;

  return;
}

void
relpal_outfile::print_second_level_results(size_t test_number, bool valid, bool reliable)
{
  const relative_pairs&   pairs = *get_analysis().get_pairs();
  const regression_model& model = get_analysis().get_current_model();

  size_t max_name_size  = get_analysis().get_output_options().param_name_max;
  string pair_count_str = long2str(get_analysis().get_pair_count(), 6);

  string test_no = long2str(test_number, 1);
  if( !reliable )
    test_no = test_no + " #";

  stringstream s;
  s << left << setw(5) << test_no;

  string t_str = s.str();
  string name = "";
  string cM = "";

  for( size_t i = 0; i < model.get_ped_parameter_count(); ++i )
  {
    const independent_variable& param = model.get_ped_parameter(i);

    if( !param.test_variable )
      continue;

    name = param.name(pairs);

    if( param.type == independent_variable::MARKER )
    {
      if( pairs.valid_distance_exist() )
        cM = " " + fp(pairs.marker_distance(pairs.marker_find(name)), 5, 1);
    }
    else
    {
      if( pairs.valid_distance_exist() )
        cM = " -----";
    }      

    break;
  }

  string cM_empty(cM.size(), ' ');
  string pc_empty(pair_count_str.size(), ' ');

  const vector<double>& T_unadj    = get_analysis().get_current_result().score_result.T;
  const vector<double>& correction = get_analysis().get_current_result().score_result.correction;
  const vector<double>& emp_pvalue = get_analysis().get_current_result().score_result.emp_pvalue;
  const vector<size_t>& rep_count  = get_analysis().get_current_result().score_result.rep_count;

  if( valid )
  {
    print_score_result_line(sum_output_stream(), t_str+name, max_name_size, cM, pair_count_str+" "+my_var_opt[0],
                            T_unadj[0], correction[0], emp_pvalue[0], rep_count[0]);

    for( size_t v = 1; v < T_unadj.size(); ++v )
      print_score_result_line(sum_output_stream(), "", max_name_size, cM_empty, pc_empty+" "+my_var_opt[v],
                              T_unadj[v], correction[v], emp_pvalue[v], rep_count[v]);
  }
  else
  {
    print_fail_line(sum_output_stream(), t_str+name, max_name_size, cM, pair_count_str);
  }

  if( detailed_out )
  {
    print_det_second_level_results(test_number, valid, reliable);

    if( valid )
    {
      print_score_result_line(det_output_stream(), name, max_name_size, cM, " "+my_var_opt[0],
                              T_unadj[0], correction[0], emp_pvalue[0], rep_count[0]);

      for( size_t v = 1; v < T_unadj.size(); ++v )
        print_score_result_line(det_output_stream(), "", max_name_size, cM_empty, " "+my_var_opt[v],
                                T_unadj[v], correction[v], emp_pvalue[v], rep_count[v]);
    }
    else
      print_fail_line(det_output_stream(), name, max_name_size, cM, "");

    det_output_stream() << endl;

    print_double_line(det_output_stream(), false);
  }

  if( export_out )
  {
    if( valid )
    {
      exp_output_stream() << name << "," << cM << "," << fp(emp_pvalue[0],9,7);
      for( size_t v = 1; v < T_unadj.size(); ++v )
        exp_output_stream() << "," << fp(emp_pvalue[v],9,7);
      exp_output_stream() << endl;
    }
    else
      exp_output_stream() << name << "," << cM << ",FAIL" << endl; 
  }

  return;
}

//--------------------------------------------------------------------

void
relpal_outfile::print_reg_result_line(ostream& out,
                                      string name,  size_t max_name_size,
                                      string cM,    string count,
                                      double T_val, double p_val)
{
  if( !out )
    return;

  out << left << setw(max_name_size) << name;
  out << cM << count;

  //out << " " << fp(beta, 12)
  //    << " " << fp(ss, 9);

  out << " " << fp(T_val, 9)
      << " " << pval(p_val, 12, 10, 3);

  out << endl;

  return;
}

void
relpal_outfile::print_reg_result_block(ostream& out,
                                       string name,         size_t max_name_size,
                                       string cM,           string count,
                                       const  matrix& beta, const matrix& var,
                                       double T_val,        double p_val)
{
  if( !out )
    return;

  size_t t_count = get_analysis().get_current_model().get_trait_count();

  for( size_t i = 0; i < t_count; ++i )
  {
    if( t_count > 1 )
      out << left << setw(max_name_size) << name + "(" + long2str(i+1,1) + ")";
    else
      out << left << setw(max_name_size) << name;

    if( i == 0 )
      out << cM << count;

    out << " " << fp(beta(i, 0), 12);

    for( size_t j = 0; j < t_count; ++j )
      out << " " << fp(var(i, j), 9);

    if( i == 0 )
      out << " " << fp(T_val, 9)
          << " " << pval(p_val, 12, 10, 3);

    out << endl;
  }

  if( t_count > 1 )
    out << endl;

  return;
}

void
relpal_outfile::print_score_result_line(ostream& out,
                                        string name,       size_t max_name_size,
                                        string cM,         string count,
                                        double T_unadj,    double correction,
                                        double emp_pvalue, size_t rep_count)
{
  if( !out )
    return;

  out << left << setw(max_name_size) << name;
  out << cM << count;

  out << " " << fp(T_unadj, 10)
      << " " << fp(T_unadj + correction, 10)
      << " " << pval(emp_pvalue, 12, 10, 3)
      << " " << right << setw(10) << rep_count;

  out << endl;

  return;
}

void
relpal_outfile::print_fail_line(ostream& out,
                                string name,  size_t max_name_size,
                                string cM,    string count)
{
  if( !out )
    return;

  out << left << setw(max_name_size) << name;
  out << cM << count;

  out << "       FAIL";

  out << endl;

  return;
}

void
relpal_outfile::get_reg_param_results(size_t i, size_t t_count,
                                      const matrix& beta, const matrix& var,
                                      matrix& param_beta, matrix& param_var)
{
  param_beta.resize(t_count, 1);
  param_var.resize(t_count, t_count);

  for( size_t t1 = 0; t1 < t_count; ++t1 )
  {
    param_beta(t1, 0) = beta(i+t1, 0);

    for( size_t t2 = 0; t2 < t_count; ++t2 )
    {
      param_var(t1, t2) = var(i+t1, i+t2);
    }
  }

  return;
}

pair<double, double>
relpal_outfile::compute_multivariate_T_value(const matrix& beta, const matrix& var)
{
  size_t df = get_analysis().get_current_model().get_trait_count();

  double T_val = QNAN;
  double P_val = QNAN;

  //T = beta^T(var^1)beta
  matrix vari, temp, T;
  SAGE::SVD svd;
  svd.compute(var);
  svd.inverse(vari);

  XTZ(beta, vari, temp);
  multiply(temp, beta, T);

#if 0
  print_matrix(beta, cout, "beta");
  print_matrix(var, cout, "var");
  print_matrix(vari, cout, "var^i");
  print_matrix(temp, cout, "temp");
  print_matrix(T, cout, "T");
#endif

  if( T.rows() )
  {
    T_val = T(0, 0);
    P_val = 1.0 - chi_square_cdf(df, T_val);
  }

  return make_pair(T_val, P_val);
}

pair<double, double>
relpal_outfile::compute_one_sided_T_value(const matrix& beta, const matrix& var)
{
  double T_val = QNAN;
  double P_val = QNAN;

  //T = beta^T(var^1)beta
  double g_beta = beta(0,0);
  double g_var  = var(0,0);

  if( g_beta < 0.0 )
    T_val = 0.0;
  else
    T_val = (g_beta * g_beta) / g_var;

  double c_val0 = 0.0;
  if( T_val > 0.0 )
    c_val0 = 0.5; 
  double c_val1 = 0.5 * chi_square_cdf(1, T_val);
  P_val = 1.0 - (c_val0 + c_val1);

#if 0
  cout << "g_beta = " << g_beta
       << ", g_var = " << g_var
       << ", T_val = " << T_val
       << ", p-val = " << P_val << endl;
#endif

  return make_pair(T_val, P_val);
}

//--------------------------------------------------------------------

void
relpal_outfile::print_det_first_level_results(size_t test_number)
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "Test " << test_number << endl;
  print_double_line(out, false);

  print_first_model_info();
  print_sample_info();

  print_reg_result_heading("Model");

  return;
}

void
relpal_outfile::print_det_second_level_results(size_t test_number, bool valid, bool reliable)
{
  ostream& out = det_output_stream();
  if( !out )
    return;

  out << "Test " << test_number;
  if( !reliable )
    out << " #";
  out << endl;
  print_double_line(out, false);

  print_second_model_info();
  print_sample_info();

  if( valid )
  {
    // Print null model regression results.
    //
    const regression_model& null_model = get_analysis().get_current_null_model();

    print_reg_result_heading("H0");
    print_reg_estimates(null_model, get_analysis().get_current_result().H0_result);

    // Print test model regression results if univariate.
    //
    if( null_model.get_trait_count() == 1 )
    {
      const regression_model& test_model = get_analysis().get_current_model();

      print_reg_result_heading("H1");
      print_reg_estimates(test_model, get_analysis().get_current_result().H1_result);
      out << endl;
    }
  }
  else
    out << endl;

  print_score_result_heading();

  return;
}

void
relpal_outfile::print_model_heading()
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  out << "-----" << endl
      << "Model" << endl
      << "-----" << endl << endl;

  return;
}

void
relpal_outfile::print_first_model_info()
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  print_model_heading();

  const regression_model& null_model = get_analysis().get_current_null_model();

  out << get_trait_names() << " ~ " << get_independent_variable_names(null_model) << endl;

  return;
}

void
relpal_outfile::print_second_model_info()
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  print_model_heading();

  const regression_model& null_model = get_analysis().get_current_null_model();
  const regression_model& test_model = get_analysis().get_current_model();

  out << "H0: " << get_trait_names() << " ~ " << get_independent_variable_names(null_model) << endl;
  out << "H1: " << get_trait_names() << " ~ " << get_independent_variable_names(test_model) << endl;

  return;
}

void
relpal_outfile::print_sample_info()
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  stringstream s;

  s << endl
    << "------------" << endl
    << "Sample Count" << endl
    << "------------" << endl << endl;

  size_t ind_count  = get_analysis().get_member_count();
  size_t pair_count = get_analysis().get_pair_count();

  s << "Number of individuals used at first level     = " << ind_count << endl;
  s << "Number of relative pairs used at second level = " << pair_count << endl;

  out << s.str();

  return;
}

void
relpal_outfile::print_reg_result_heading(string h)
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  stringstream d;
  stringstream s;

  d << "---------------";
  for( size_t i = 0; i < h.size(); ++i )
    d << "-";
  d << endl;

  s << "Estimates from " << h << endl;

  out << endl << d.str() << s.str() << d.str() << endl;

  print_reg_table_heading(out, false, false, false);
  print_reg_single_line(out, false, false, false);

  return;
}

void
relpal_outfile::print_score_result_heading()
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  stringstream s;

  s << "------------------" << endl
    << "Score Test Results" << endl
    << "------------------" << endl << endl;

  out << s.str();

  print_score_table_heading(out, false, true, false);
  print_score_single_line(out, true, false);

  return;
}

void
relpal_outfile::print_reg_estimates(const regression_model&  model,
                                    const regression_result& reg_result)
{
  ostream& out = det_output_stream();

  if( !out )
    return;

  const relative_pairs& pairs = *get_analysis().get_pairs();

  size_t max_name_size = get_analysis().get_output_options().param_name_max;

  size_t t_count = model.get_trait_count();

  for( size_t i = 0; i < model.get_ind_parameter_count(); )
  {
    const independent_variable& param = model.get_ind_parameter(i);

    string name = param.name(pairs);

    i += t_count;

    if( name == "INTERCEPT" )
      continue;

    matrix beta, var;

    get_reg_param_results(i-t_count, t_count,
                          reg_result.ind_beta,
                          reg_result.ind_variance,
                          beta, var);

    pair<double, double> T_P = compute_multivariate_T_value(beta, var);

    print_reg_result_block(out, name, max_name_size, "", "", beta, var, T_P.first, T_P.second);
  }

  if( t_count == 1 )
  {
    for( size_t i = 0; i < model.get_ped_parameter_count(); ++i )
    {
      const independent_variable& param = model.get_ped_parameter(i);

      string name = param.name(pairs);

      if( name == "POLYGENIC_EFF" || name == "RANDOM_EFF" )
        continue;

      matrix beta, var;

      get_reg_param_results(i, t_count,
                            reg_result.ped_beta,
                            reg_result.ped_variance,
                            beta, var);

      pair<double, double> T_P = compute_one_sided_T_value(beta, var);

      print_reg_result_block(out, name, max_name_size, "", "", beta, var, T_P.first, T_P.second);

      i += (((t_count+1)*t_count)/2);
    }
  }

  return;
}

string
relpal_outfile::get_trait_names() const
{
  const relative_pairs&   pairs = *get_analysis().get_pairs();
  const regression_model& model = get_analysis().get_current_model();

  string t_names = pairs.trait_name(model.get_trait(0).trait_index);

  for( size_t t = 1; t < model.get_trait_count(); ++t )
  {
    size_t ti = model.get_trait(t).trait_index;

    t_names += (", " + pairs.trait_name(ti));
  }

  return t_names;
}

string
relpal_outfile::get_independent_variable_names(const regression_model& model) const
{
  const relative_pairs&   pairs = *get_analysis().get_pairs();

  size_t t_count =  model.get_trait_count();

  string var_name = "[Intercept";

  for( size_t i = t_count; i < model.get_ind_parameter_count(); )
  {
    const independent_variable& param = model.get_ind_parameter(i);

    var_name += (" + " + param.name(pairs));

    i += t_count;
  }

  var_name += "]";

  size_t g_count = ((t_count+1)*t_count)/2;

  for( size_t i = 0; i < model.get_ped_parameter_count(); )
  {
    const independent_variable& param = model.get_ped_parameter(i);

    var_name += (" + " + param.name(pairs));

    i += g_count;
  }

  return var_name;
}

} // end of namespace RELPAL
} // end of namespace SAGE

