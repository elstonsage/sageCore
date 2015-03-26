//****************************************************************************
//* File:      lodpal_out_text.cpp                                           *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation.               kbj         *
//*                    1.0 Modification on output look.          yjs Jun. 01 *
//*                    2.0 Modification for X-linkage            yjs May. 02 *
//*                    3.0 Parent-of-Origin added.               yjs Feb. 03 *
//*                                                                          *
//* Notes:     This file implemnets lodpal_test_textfile class.              *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_out.h"

namespace SAGE   {
namespace LODPAL {

void
lodpal_test_viewer::print_title_heading(const ARP_base_analysis& test, bool split)
{
  ostream& out = output_stream();

  string pair_type = "Relative";
  if( test.pairs_info().parameters().sib_pairs_only() )
    pair_type = "Sib";

  out << "  Conditional Logistic Analysis of";
  if( split && !wide_output() )
  {
    if( test.parameters().covariate_count() < 2 )
    {
      if( !test.parameters().autosomal_model().parent_of_origin )
        out << endl << "    Affected " << pair_type << " Pairs";
      else if( test.parameters().autosomal_model().fixed != autosomal_model_type::none )
        out << endl << "    Affected " << pair_type << " Pairs";
      else if( test.parameters().covariate_count() < 1 )
        out << endl << "    Affected " << pair_type << " Pairs";
      else
        out << " Affected " << pair_type << " Pairs";
    }
    else
      out << " Affected " << pair_type << " Pairs";
  }
  else
    out << " Affected " << pair_type << " Pairs";

  if( test.pairs_info().parameters().multipoint() )
    out << " - multipoint";
  else
    out << " - singlepoint";

  out << endl;
}

void
lodpal_test_viewer::print_model_summary(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  string alignment = "                ";

  if( test.pairs_info().parameters().trait_count() )
  {
    string pair_type = "relative";
    if( test.pairs_info().parameters().sib_pairs_only() )
      pair_type = "sib";

    out << "    Trait    : ";

    for(size_t i=0; i < test.pairs_info().parameters().trait_count(); ++i)
    {
      if( i )
        out << alignment;
      out << test.pairs_info().parameters().trait_parameters(i).name(test.relative_pairs()) << endl;
      out << alignment;

      switch( test.pairs_info().parameters().trait_parameters(i).pair_select )
      {
        case trait_parameter::conaff:
             out << "concordantly affected " << pair_type << " pairs" << endl;
             break;

        case trait_parameter::condisc:
             out << "concordantly affected " << pair_type << " pairs" << endl
                 << alignment
                 << " with discordantly affected sib pairs" << endl
                 << alignment
                 << " and concordantly unaffected sib pairs" << endl;
             break;

        case trait_parameter::noconunaff:
             out << "concordantly affected " << pair_type << " pairs" << endl
                 << alignment
                 << " with discordantly affected sib pairs" << endl;
             break;

        case trait_parameter::contrast:
             out << "concordantly affected " << pair_type << " pairs" << endl
                 << alignment
                 << " contrast to discordantly affected " << pair_type << " pairs" << endl;
      }
      const RPED::RefTraitInfo& info = test.relative_pairs().fped_info().trait_info(
                                 test.pairs_info().parameters().trait_parameters(i).trait);
      if( info.type() != RPED::RefTraitInfo::binary_trait )
      {
        out << alignment;
        out << "cutpoint = " << test.pairs_info().parameters().trait_parameters(i).cutpoint << endl;
      }
    }
    out << endl;
  }

  if( test.pairs_info().parameters().covariate_count() )
  {
    out << "    Covariate: ";

    for(size_t i=0; i < test.pairs_info().parameters().covariate_count(); ++i)
    {
      if( test.pairs_info().parameters().covariate_count() > 1 )
      {
        if( i )
          out << endl << alignment;
        out << i+1 << ". " << test.pairs_info().parameters().covariate_parameters(i).name(test.relative_pairs());
      }
      else
        out << test.pairs_info().parameters().covariate_parameters(i).name(test.relative_pairs());

      switch( test.pairs_info().parameters().covariate_parameters(i).operation )
      {
        case covariate_type::none:
             out << " discordance status" << endl;
             break;

        case covariate_type::sum:
             out << endl << alignment
                 << "sum of two individual covariate values" << endl;
             break;

        case covariate_type::diff:
             out << endl << alignment
                 << "difference of two individual covariate values" << endl;
             break;

        case covariate_type::prod:
             out << endl << alignment
                 << "product of two individual covariate values" << endl;
             break;

        case covariate_type::single:
             out << endl << alignment
                 << "single covariate value" << endl;
             break;

        case covariate_type::pair:
             out << endl << alignment
                 << "pair covariate value" << endl;
             break;

        case covariate_type::avg:
             out << endl << alignment
                 << "average of two individual covariate values" << endl;
             break;

        default: break;
      }

      double ad_val      = test.pairs_info().parameters().covariate_parameters(i).adjust_value;
      double weight_sum  = test.pairs_info().parameters().weight_parameter().info.sum();
      double info_mean   = test.pairs_info().parameters().covariate_parameters(i).info.sum()/
                           weight_sum;
      double info_y_mean = test.pairs_info().parameters().covariate_parameters(i).info_y.sum()/
                           weight_sum;
//      double y_variance  = test.pairs_info().parameters().covariate_parameters(i).info_y.variance();
      double y_variance  = test.pairs_info().parameters().covariate_parameters(i).variance_y;

      out << alignment;

      string adjust_method = "centering";
      if( test.pairs_info().parameters().covariate_parameters(i).adjust == covariate_type::minimum )
        adjust_method = "adjusting";

      stringstream weighted;
      weighted  << alignment
                << "weighted mean before " << adjust_method << " = " << fp(info_mean, 8, 7) << "\n"
                << alignment
                << "weighted mean after " << adjust_method<< "  = " << fp(info_y_mean, 8, 7) << "\n"
                << alignment
                << "weighted std. deviation        = " << fp(sqrt(y_variance), 8, 7) << "\n";

      stringstream unweighted;
      unweighted << alignment
                 << "mean before " << adjust_method << " = " << fp(info_mean, 8, 7) << "\n"
                 << alignment
                 << "mean after " << adjust_method << "  = " << fp(info_y_mean, 8, 7) << "\n"
                 << alignment
                 << "std. deviation        = " << fp(sqrt(y_variance), 8, 7) << "\n";

      switch( test.pairs_info().parameters().covariate_parameters(i).adjust )
      {
        case covariate_type::mean:
             out << "mean centered " << endl;

             if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
               out << weighted.str();
             else
               out << unweighted.str();

             break;

        case covariate_type::minimum:
             out << "minimum adjusted " << endl;

             if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
               out << weighted.str();
             else
               out << unweighted.str();

             break;

        case covariate_type::prop:
             out << "proportion = " << ad_val << endl;

             if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
               out << weighted.str();
             else
               out << unweighted.str();

             break;

        case covariate_type::dsp:
             out << "mean           = " << info_mean << endl
                 << "                "
                 << "std. deviation = "
                 << fp(sqrt(test.pairs_info().parameters().covariate_parameters(i).info.variance()), 8, 7) << endl;
             break;
      }
    }
    out << endl;
  }

  if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
  {
    double weight_sum  = test.pairs_info().parameters().weight_parameter().info.sum();
    out << "    Weight   : ";
    out << test.pairs_info().parameters().weight_parameter().name(test.relative_pairs()) << endl;
    out << "                "
        << "sum of weights = " << weight_sum << endl;
    out << endl;
  }

  if( test.pairs_info().parameters().subset_count() )
  {
    out << "    Subset   : ";

    for(size_t i=0; i < test.pairs_info().parameters().subset_count(); ++i)
    {
      if( test.pairs_info().parameters().subset_count() > 1 )
      {
        if( i )
          out << "               ";
        out << i+1 << ". " << test.pairs_info().parameters().subset_parameters(i).name(test.relative_pairs()) << endl;
      }
      else
        out << test.pairs_info().parameters().subset_parameters(i).name(test.relative_pairs()) << endl;

    }
    out << endl;
  }

  if( test.pairs_info().parameters().turn_off_default() )
    out << "    Method   : non default analysis method";
  else
    out << "    Method   : default analysis method";
  out << endl << endl;
}

void
lodpal_test_viewer::print_skip_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  lodpal_parameters::marker_const_iterator pi = test.pairs_info().parameters().marker_begin();
  for( ; pi != test.pairs_info().parameters().marker_end(); ++pi)
  {
    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << test.relative_pairs().marker_name(pi->marker)
          << fp(test.relative_pairs().marker_distance(pi->marker), 6, 1) << " ";
    else
      out << setw(my_max_marker_name)
          << test.relative_pairs().marker_name(pi->marker) << " ";

    out << "  SKIPPED";
    out << endl;
  }
}

void
lodpal_test_viewer::set_format_params(const ARP_base_analysis& test)
{
  if( test.relative_pairs().valid_distance_exist() )
    my_max_marker_name = std::max(test.pairs_info().parameters().max_marker_name_size() + 6, my_max_marker_name);
  else
    my_max_marker_name = std::max(test.pairs_info().parameters().max_marker_name_size(), my_max_marker_name);

  my_max_ped_name = std::max(test.pairs_info().max_ped_name_size(), my_max_ped_name);
  my_max_ind_name = std::max(test.pairs_info().max_ind_name_size(), my_max_ind_name);
}

//---------------------------------------------------------------------------------

void
lodpal_test_textfile::print_header(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  set_format_params(test);

  print_double_line(test);

  print_title_heading(test, true);

  print_double_line(test);

  out << endl;

  print_model_summary(test);

  out << "    Model    : " << test.pairs_info().parameters().autosomal_model().name();
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
  {
    out << test.pairs_info().parameters().autosomal_model().poo_name();

    if( !test.pairs_info().parameters().sib_pairs_only() )
      out << endl << endl
          << "             ! WARNING: MODELS OF PARENT-OF-ORIGIN EFFECTS USING NON-SIB PAIRS"
          << endl
          << "                        MAY BE INVALID UNLESS ASCERTAINMENT CONDITIONS ARE MET.";
  }
  out << endl << endl;

  if( is_valid_for_emp_p_value(test) )
    out << "    Legend   :" << endl
        << "      *   - significance  .05 level;" << endl
        << "      **  - significance  .01 level;" << endl
        << "      *** - significance .001 level;" << endl << endl;

  print_double_line(test);

  print_table_heading(test);

  print_single_line(test);
}

double
get_asym_p_value(double lod_score, size_t cov_count)
{
  double p_k = 0.0;

  if( cov_count != 0 )
    p_k = 1.0 - chi_square_cdf(cov_count, 2.0*log(10.)*lod_score);  

  double p_k1 = 1.0 - chi_square_cdf(cov_count+1, 2.0*log(10.)*lod_score);

  return 0.5*p_k + 0.5*p_k1;
}

double
get_emp_p_value(double alpha, double ss, double cM, size_t cov_count, bool multi)
{
  double _log_alpha = -log(alpha);

#if 0
  cout << alpha
       << ", " << _log_alpha
       << ", " << ss
       << ", " << cM
       << ", " << cov_count;
#endif

  if( multi )
  {
    double _log_p = -1.4232 + 0.7273    * _log_alpha
                            + 0.00663   * ss
                            - 0.0058    * cM
                            - 0.000007  * (ss*ss)
                            - 0.000028  * (ss*_log_alpha)
                            + 0.00111   * (cM*_log_alpha)
                            + 0.0000212 * (ss*cM);

    if( cov_count == 4 )
      return _log_p;
    else if( cov_count == 3 )
      return _log_p + 0.1989;
    else if( cov_count == 2 )
      return _log_p + 0.4837;
    else if( cov_count == 1 )
      return _log_p + 0.9228;
    else if( cov_count == 0 )
      return _log_p + 1.605;
    else
      return numeric_limits<double>::quiet_NaN();
  }

  // Singlepoint
  //
  double _log_p = -1.4429 + 0.7522     * _log_alpha
                          + 0.00602    * ss
                          - 0.00000549 * (ss*ss);

  if( cov_count == 4 )
    return _log_p;
  else if( cov_count == 3 )
    return _log_p + 0.1968;
  else if( cov_count == 2 )
    return _log_p + 0.4747;
  else if( cov_count == 1 )
    return _log_p + 0.88688;
  else if( cov_count == 0 )
    return _log_p + 1.5434;
  else
    return numeric_limits<double>::quiet_NaN();
}

void
lodpal_test_textfile::print_results(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  lodpal_parameters::marker_const_iterator pi = test.parameters().marker_begin();
  for( ; pi != test.parameters().marker_end(); ++pi)
  {
    out << " ";
    out.setf(ios_base::left, ios_base::adjustfield);

    if( test.relative_pairs().valid_distance_exist() )
      out << setw(my_max_marker_name - 6)
          << test.relative_pairs().marker_name(pi->marker)
          << fp(test.relative_pairs().marker_distance(pi->marker), 6, 1) << " ";
    else
      out << setw(my_max_marker_name)
          << test.relative_pairs().marker_name(pi->marker) << " ";

    out << fp(test.get_lodpal_result().lod_score(),9,6)               << " ";

    out.setf(ios_base::right, ios_base::adjustfield);
    out << setw(6)   << test.pairs_info().fsib_pair_count()           << " ";
    if( wide_output() )
      out << setw(6) << test.pairs_info().hsib_pair_count()           << " "
          << setw(6) << test.pairs_info().other_pair_count()          << " ";
    
    out << setw(6)   << test.pairs_info().pair_count()                << " ";

    bool lambda = test.pairs_info().parameters().print_lambda();

    if( lambda )
    {
      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          out << " " << fp(pi->lambda1m.value(), 9, 6);
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          out << " " << fp(pi->lambda1p.value(), 9, 6);
      }
      else
        out << " " << fp(pi->lambda1.value(), 9, 6);

      if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
        out << " " << fp(pi->lambda2.value(), 9, 6);
    }
    else
    {
      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          out << " " << fp(pi->beta1m.value(), 9, 6);
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          out << " " << fp(pi->beta1p.value(), 9, 6);
      }
      else
        out << " " << fp(pi->beta1.value(), 9, 6);

      if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
        out << " " << fp(pi->beta2.value(), 9, 6);

      lodpal_parameters::covariate_const_iterator ci = test.parameters().covariate_begin();
      for( ; ci != test.parameters().covariate_end(); ++ci )
      {
        if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
        {
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
            out << " " << fp(ci->delta1m.value(), 9, 6);
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
            out << " " <<fp(ci->delta1p.value(), 9, 6);
        }
        else
          out << " " << fp(ci->delta1.value(), 9, 6);

        if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
          out << " " << fp(ci->delta2.value(), 9, 6);
      }
    }

    if( wide_output() )
    {
      if( lambda )
      {
        out << " ";
        if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
        {
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
            out << " " << fp(pi->lambda1m.first_derivative(), 9, 6);
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
            out << " " << fp(pi->lambda1p.first_derivative(), 9, 6);
        }
        else
          out << " " << fp(pi->lambda1.first_derivative(), 9, 6);

        if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
          out << " " << fp(pi->lambda2.first_derivative(), 9, 6);
      }
      else
      {
        out << " ";
        if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
        {
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
            out << " " << fp(pi->beta1m.first_derivative(), 9, 6);
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
            out << " " << fp(pi->beta1p.first_derivative(), 9, 6);
        }
        else
          out << " " << fp(pi->beta1.first_derivative(), 9, 6);

        if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
          out << " " << fp(pi->beta2.first_derivative(), 9, 6);

        lodpal_parameters::covariate_const_iterator ci = test.parameters().covariate_begin();
        for( ; ci != test.parameters().covariate_end(); ++ci )
        {
          if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
          {
            if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
              out << " " << fp(ci->delta1m.first_derivative(), 9, 6);
            if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
              out << " " << fp(ci->delta1p.first_derivative(), 9, 6);
          }
          else
            out << " " << fp(ci->delta1.first_derivative(), 9, 6);
        }
      }
    }

    // Condition to compute Emp p-value
    //  one-parameter model
    //  sample size : 20-320
    //  covariate count : 0 - 4
    //  average cM (multipoint) :1 - 20

#if 0
cout << "is_valid_for_emp_p_value? "
     << "  " << test.relative_pairs().is_x_linked(test.parameters().marker_parameters(0).marker)
     << ", " << test.pairs_info().parameters().autosomal_model().model
     << ", " << test.pairs_info().parameters().autosomal_model().parent_of_origin
     << ", " << test.pairs_info().parameters().trait_parameters(0).pair_select
     << ", " << test.pairs_info().pair_count()
     << ", " << test.pairs_info().parameters().covariate_count()
     << ", " << test.pairs_info().parameters().multipoint()
     << ", " << test.relative_pairs().get_average_marker_distance()
     << endl;
#endif

    if( is_valid_for_emp_p_value(test) )
    {
      double lod = test.get_lodpal_result().lod_score();

      if( SAGE::isnan(lod) )
        out << "  ---------      ---------";
      
      else
      {
        size_t cov_cnt = test.pairs_info().parameters().covariate_count();
        double alpha   = get_asym_p_value(lod, cov_cnt);
        double ss      = (double)test.pairs_info().pair_count();
        double cM      = test.relative_pairs().get_average_marker_distance();
        bool   multi   = test.pairs_info().parameters().multipoint();

        double emp_pval  = get_emp_p_value(alpha, ss, cM, cov_cnt, multi);

        if( emp_pval < 0. )
          emp_pval = 0.;

        if( get_pval_scientific_notation() )
          out << "  " << pval_scientific(alpha, 10, 3, 3) << "  " << pval_scientific(exp(-emp_pval), 10, 3, 3);
        else
          out << "  " << pval(alpha, 12, 10, 3) << "  " << pval(exp(-emp_pval), 12, 10, 3);
      }
    }

#if 0
    if( wide_output() )
    {
      out << " " << setw(4) << test.get_lodpal_result().function_evaluations()
          << " " << setw(3) << test.get_lodpal_result().last_error();

      out << " " << setw(4) << test.get_lodpal_result().ivfl()
          << " " << setw(4) << test.get_lodpal_result().igage()
          << " " <<            test.get_lodpal_result().method();
    }
#endif

    out << endl;
  }
}

void
lodpal_test_textfile::print_footer(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << endl << endl;
}

void
lodpal_test_textfile::print_table_heading_single_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << "----------";

  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "----------";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "----------";

  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "----------";

  return;
}

void
lodpal_test_textfile::print_param_names(const ARP_base_analysis& test, bool lambda)
{
  ostream& out = output_stream();

  if( lambda )
  {
    if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    {
      if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        out << " Lambda1m ";
      if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        out << " Lambda1p ";
    }
    else
      out << " Lambda1  ";

    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << " Lambda2  ";

    return;
  }

  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
  {
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      out << " Beta1m   ";
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      out << " Beta1p   ";
  }
  else
    out << " Beta1    ";

  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << " Beta2    ";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
  {
    if( test.pairs_info().parameters().covariate_parameters(ci).operation == covariate_type::none )
    {
      string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());

      string temp = c_name.substr(0,2);
      if( c_name.size() < 2 )
      {
        temp = c_name;
      }

      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          string temp_m = temp + " disc_m";
          out << " " << setw(9) << temp_m;
        }
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          string temp_p = temp + " disc_p";
          out << " " << setw(9) << temp_p;
        }
      }
      else
      {
        temp += " discor";
        out << " " << setw(9) << temp;
      }
    }
    else
    {
      string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
      out.setf(ios_base::left, ios_base::adjustfield);

      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        string temp = c_name.substr(0,7);
        if( c_name.size() < 7 )
        {
          temp = c_name;
        }

        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          string temp_m = temp + "_m";
          out << " " << setw(9) << temp_m;
        }
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          string temp_p = temp + "_p";
          out << " " << setw(9) << temp_p;
        }
      }
      else
        out << " " << setw(9) << c_name.substr(0,9);
    }
  }
  
  return;
}

void
lodpal_test_textfile::print_table_heading(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << " ";
  out << "           Full   ";
  if( wide_output() )
    out << "Half          ";
  if( test.pairs_info().parameters().covariate_count() )
    out << "        Parameter Estimates";
  else
  {
    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << "        Parameter Estimates";
    else if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
             && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      out << "        Parameter Estimates";
    else
      out << "        Param Est";
  }

  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
  {
    if( test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << "                    ";
    else
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << "          ";
  }
  else
    for( size_t ci = 1; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "          ";

  if( wide_output() )
  {
    if( test.pairs_info().parameters().covariate_count() )
      out << "  First Derivatives";
    else
    {
      if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
        out << "  First Derivatives";
      else if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
               && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
        out << "  First Derivatives";
      else
        out << "  1st Deriv";
    }
  }

  out << endl;

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << " ";
  out << " LOD       Sib    ";
  if( wide_output() )
    out << "Sib    Other  ";
  out << "All     ";
  out << "---------";

  print_table_heading_single_line(test);

/*
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << "----------";

  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "----------";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "----------";
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "----------";
*/
  if( wide_output() )
  {
    out << "  ---------";

    print_table_heading_single_line(test);

/*
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      out << "----------";
    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << "----------";
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "----------";
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << "----------";
*/
  }

  if( is_valid_for_emp_p_value(test) )
      out << "  Aymp.          Emp.";

  out << endl;

  out.setf(ios_base::left, ios_base::adjustfield);

  if( test.relative_pairs().valid_distance_exist() )
    out << setw(my_max_marker_name-4) << "MARKER" << setw(5) << "cM";
  else
    out << setw(my_max_marker_name+1) << "MARKER";
  out << " SCORE     Pairs  ";
  if( wide_output() )
    out << "Pairs  Pairs  ";
  out << "Pairs  ";

  bool lambda = test.pairs_info().parameters().print_lambda();

  print_param_names(test, lambda);
/*
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
  {
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      out << " Beta1m   ";
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      out << " Beta1p   ";
  }
  else
    out << " Beta1    ";

  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << " Beta2    ";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
  {
    if( test.pairs_info().parameters().covariate_parameters(ci).operation == covariate_type::none )
    {
      string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());

      string temp = c_name.substr(0,2);
      if( c_name.size() < 2 )
      {
        temp = c_name;
      }

      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          string temp_m = temp + " disc_m";
          out << " " << setw(9) << temp_m;
        }
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          string temp_p = temp + " disc_p";
          out << " " << setw(9) << temp_p;
        }
      }
      else
      {
        temp += " discor";
        out << " " << setw(9) << temp;
      }
    }
    else
    {
      string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
      out.setf(ios_base::left, ios_base::adjustfield);

      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        string temp = c_name.substr(0,7);
        if( c_name.size() < 7 )
        {
          temp = c_name;
        }

        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          string temp_m = temp + "_m";
          out << " " << setw(9) << temp_m;
        }
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          string temp_p = temp + "_p";
          out << " " << setw(9) << temp_p;
        }
      }
      else
        out << " " << setw(9) << c_name.substr(0,9);
    }
  }
*/

  if( wide_output() )
  {
    out << " ";
    print_param_names(test, lambda);
/*
    if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    {
      if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        out << " Beta1m   ";
      if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        out << " Beta1p   ";
    }
    else
      out << " Beta1    ";

    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << " Beta2    ";

    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    {
      if( test.pairs_info().parameters().covariate_parameters(ci).operation == covariate_type::none )
      {
        string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());

        string temp = c_name.substr(0,2);
        if( c_name.size() < 2 )
        {
          temp = c_name;
        }

        if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
        {
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          {
            string temp_m = temp + " disc_m";
            out << " " << setw(9) << temp_m;
          }
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          {
            string temp_p = temp + " disc_p";
            out << " " << setw(9) << temp_p;
          }
        }
        else
        {
          temp += " discor";
          out << " " << setw(9) << temp;
        }
      }
      else
      {
        string c_name = test.pairs_info().parameters().covariate_parameters(ci).name(test.relative_pairs());
        out.setf(ios_base::left, ios_base::adjustfield);

        if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
        {
          string temp = c_name.substr(0,7);
          if( c_name.size() < 7 )
          {
            temp = c_name;
          }

          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
          {
            string temp_m = temp + "_m";
            out << " " << setw(9) << temp_m;
          }
          if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
          {
            string temp_p = temp + "_p";
            out << " " << setw(9) << temp_p;
          }
        }
        else
          out << " " << setw(9) << c_name.substr(0,9);
      }
    }
*/
  }

  if( is_valid_for_emp_p_value(test) )
      out << "  P-value        P-value";

#if 0
  if( wide_output() )
  {
    out << " NFE  LFL ";
    out << " VFL  GAGE M";
  }
#endif

  out << endl;
}

void
lodpal_test_textfile::print_param_double_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << "==========";
  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "==========";
  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "==========";
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "==========";

  return;
}

void
lodpal_test_textfile::print_double_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << "=";
  out << "==================";
  if( wide_output() )
    out << "==============";
  out << "=================";

  print_param_double_line(test);
/*
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << "==========";
  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "==========";
  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "==========";
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "==========";
*/
  if( wide_output() )
  {
    out << "===========";

    print_param_double_line(test);
/*
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      out << "==========";
    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << "==========";
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << "==========";
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << "==========";
*/
#if 0
    out << "=========";
    out << "============";
#endif
  }

  if( is_valid_for_emp_p_value(test) )
      out << "==========================";

  out << endl;
}

void
lodpal_test_textfile::print_param_single_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << " ---------";
  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << " ---------";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << " ---------";
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << " ---------";

  return;
}

void
lodpal_test_textfile::print_single_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i <= my_max_marker_name; ++i )
    out << "-";
  out << " --------- ------ ";
  if( wide_output() )
    out << "------ ------ ";
  out << "------  ---------";

  print_param_single_line(test);
/*
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    out << " ---------";
  if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << " ---------";

  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << " ---------";
  if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
      && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << " ---------";
*/
  if( wide_output() )
  {
    out << "  ---------";

    print_param_single_line(test);
/*
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      out << " ---------";
    if( !test.pairs_info().parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      out << " ---------";
    for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
      out << " ---------";
    if(    test.pairs_info().parameters().autosomal_model().parent_of_origin
        && test.pairs_info().parameters().autosomal_model().fixed == autosomal_model_type::none )
      for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
        out << " ---------";
*/
#if 0
    out << " ---- ---";
    out << " ---- ---- -";
#endif
  }

  if( is_valid_for_emp_p_value(test) )
      out << "  ---------      ---------";

  out << endl;
}

bool
lodpal_test_textfile::is_valid_for_emp_p_value(const ARP_base_analysis& test) const
{
  if( test.relative_pairs().is_x_linked(test.parameters().marker_parameters(0).marker) )
    return false;

  if(    test.pairs_info().parameters().autosomal_model().model != autosomal_model_type::one_parameter
      || test.pairs_info().parameters().autosomal_model().parent_of_origin )
    return false;

  if( test.pairs_info().parameters().trait_parameters(0).pair_select != trait_parameter::conaff )
    return false;

  if( test.pairs_info().pair_count() < 20 || test.pairs_info().pair_count() > 350 )
    return false;

  if( test.pairs_info().parameters().covariate_count() > 4 )
    return false;

  if( test.pairs_info().parameters().multipoint() &&
      (test.relative_pairs().get_average_marker_distance() < 1.0 ||
       test.relative_pairs().get_average_marker_distance() > 20.0) )
    return false;

  return true;
}

} // end of namespace LODPAL
} // end of namespace SAGE
