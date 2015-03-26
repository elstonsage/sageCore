//****************************************************************************
//* File:      lodpal_out_diag.cpp                                           *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Jun. 01 *
//*                    0.1 Histogram added                       yjs Aug. 01 *
//*                    1.0 Modified for X-linkage.               yjs May. 02 *
//*                                                                          *
//* Notes:     This file implemnets lodpal_test_diagfile class.              *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "numerics/histogram.h"
#include "lodpal/lodpal_out.h"

namespace SAGE   {
namespace LODPAL {

//---------------------------------------------------------------------------------

void
lodpal_test_diagfile::print_header(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  set_format_params(test);

  print_double_line(test);

  print_title_heading(test, false);

  print_double_line(test);

  out << endl;

  print_model_summary(test);

  out << "    Model    : ";

  if( test.relative_pairs().is_x_linked(test.parameters().diagnostic_marker()) )
    out << test.parameters().x_linkage_model().name();
  else
  {
    out << test.parameters().autosomal_model().name();

    if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    {
      out << test.pairs_info().parameters().autosomal_model().poo_name();

      if( !test.pairs_info().parameters().sib_pairs_only() )
        out << endl << endl
            << "             ! WARNING: MODELS OF PARENT-OF-ORIGIN EFFECTS USING NON-SIB PAIRS"
            << endl
            << "                        MAY BE INVALID UNLESS ASCERTAINMENT CONDITIONS ARE MET.";
    }
  }
  out << endl << endl;

  out << "    Location : ";
  out << test.relative_pairs().marker_name(test.parameters().diagnostic_marker());
  if( test.relative_pairs().is_x_linked(test.parameters().diagnostic_marker()) )
    out << ", x_linked";
  out << endl << endl;

  print_double_line(test);

  out << endl;
}

void
lodpal_test_diagfile::print_results(const ARP_base_analysis& test)
{
  print_final_result_summary(test);

  print_histogram(test);

  print_individual_lod_table(test);
}

void
lodpal_test_diagfile::print_footer(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  print_double_line(test);

  out << endl << endl;
}

void
lodpal_test_diagfile::print_final_result_summary(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  out << "# Final Result Summary" << endl << endl;

  print_param_estimates(test);
  print_var_cov_matrix(test);
}

void
lodpal_test_diagfile::print_param_estimates(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  out << "    Parameter Estimates:" << endl;
  size_t p = 0;
  lodpal_parameters::marker_const_iterator    pi = test.parameters().marker_begin();
  for( ; pi != test.parameters().marker_end(); ++pi )
  {
    if( test.relative_pairs().is_x_linked(pi->marker) )
    {
      if( !test.parameters().x_linkage_model().lambda1_equal )
      {
        if( test.parameters().use_mm_pair() )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << "beta 1 male-male     = " << fp(pi->beta1mm.value(), 9, 6) << endl;
        }

        if( test.parameters().use_mf_pair() )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << "beta 1 male-female   = " << fp(pi->beta1mf.value(), 9, 6) << endl;
        }

        if( test.parameters().use_ff_pair() )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << "beta 1 female-female = " << fp(pi->beta1ff.value(), 9, 6) << endl;
        }
      }
      else
      {
        out << "      " << p+1 << ". ";
        p++;
        out << "beta 1 = " << fp(pi->beta1mm.value(), 9, 6) << endl;
      }

      if(     test.parameters().use_ff_pair()
          && !test.parameters().x_linkage_model().lambda2_fixed )
      {
        out << "      " << p+1 << ". ";
        p++;
        out << "beta 2 female-female = " << fp(pi->beta2ff.value(), 9, 6) << endl;
      }
    }
    else
    {
      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {
        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << "beta1m = " << fp(pi->beta1m.value(), 9, 6) << endl;
        }

        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << "beta1p = " << fp(pi->beta1p.value(), 9, 6) << endl;
        }
      }
      else
      {
        out << "      " << p+1 << ". ";
        p++;
        out << "beta 1 = " << fp(pi->beta1.value(), 9, 6) << endl;
      }

      if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      {
        out << "      " << p+1 << ". ";
        p++;
        out << "beta 2 = " << fp(pi->beta2.value(), 9, 6) << endl;
      }
    }
  }

  for(size_t i = 0; i < test.parameters().covariate_count(); ++i)
  {
    if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    {
      if( !test.parameters().x_linkage_model().lambda1_equal )
      {
        if( test.parameters().use_mm_pair() )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
          if(    test.parameters().covariate_parameters(i).operation == covariate_type::none )
            out << " discordance status";
          out << "(delta 1 male-male)     = " << fp(test.parameters().covariate_parameters(i).delta1mm.value(), 9, 6) << endl;
        }

        if( test.parameters().use_mf_pair() )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
          if(    test.parameters().covariate_parameters(i).operation == covariate_type::none )
            out << " discordance status";
          out << "(delta 1 male-female)   = " << fp(test.parameters().covariate_parameters(i).delta1mf.value(), 9, 6) << endl;
        }

        if( test.parameters().use_ff_pair() )
        {        out << "      " << p+1 << ". ";
          p++;
          out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
          if(    test.parameters().covariate_parameters(i).operation == covariate_type::none )
            out << " discordance status";
          out << "(delta 1 female-female) = " << fp(test.parameters().covariate_parameters(i).delta1ff.value(), 9, 6) << endl;
        }
      }
      else
      {
        out << "      " << p+1 << ". ";
        p++;
        out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
        if(    test.parameters().covariate_parameters(i).operation == covariate_type::none )
          out << " discordance status";
        out << "(delta 1) = " << fp(test.parameters().covariate_parameters(i).delta1mm.value(), 9, 6) << endl;
      }

      if(     test.parameters().use_ff_pair()
          && !test.parameters().x_linkage_model().lambda2_fixed )
      {
        out << "      " << p+1 << ". ";
        p++;
        out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
        out << "(delta 2 female-female) = " << fp(test.parameters().covariate_parameters(i).delta2ff.value(), 9, 6) << endl;
      }
    }
    else
    {
      if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      {

        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
          out << "(delta1m) = " << fp(test.parameters().covariate_parameters(i).delta1m.value(), 9, 6) << endl;
        }

        if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
        {
          out << "      " << p+1 << ". ";
          p++;
          out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
          out << "(delta1p) = " << fp(test.parameters().covariate_parameters(i).delta1p.value(), 9, 6) << endl;
        }
      }
      else
      {
        out << "      " << p+1 << ". ";
        p++;
        out << test.parameters().covariate_parameters(i).name(test.relative_pairs());

        if(    test.parameters().covariate_parameters(i).operation == covariate_type::none )
          out << " discordance status";
        out << "(delta 1) = " << fp(test.parameters().covariate_parameters(i).delta1.value(), 9, 6) << endl;
      }
      if( !test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
      {
        out << "      " << p+1 << ". ";
        p++;
        out << test.parameters().covariate_parameters(i).name(test.relative_pairs());
        out << "(delta 2) = " << fp(test.parameters().covariate_parameters(i).delta2.value(), 9, 6) << endl;
      }
    }
  }
  out << endl;
}

void
lodpal_test_diagfile::print_var_cov_matrix(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  out << "    Variance-Covariance Matrix(assuming independent pairs):" << endl;

  if( test.get_lodpal_result().ivfl() > 2 )
  {
    out << "      No valid variance-covariance matrix exist." << endl << endl;
    return;
  }

  size_t col_limit = test.get_lodpal_result().var_cov_matrix().cols();
  if( !test.parameters().covariate_count() )
    col_limit /= 2;

  size_t row_limit = test.get_lodpal_result().var_cov_matrix().rows();
  if( !test.parameters().covariate_count() )
    row_limit /= 2;

  out << "      ------------";
  for( size_t c = 0; c < col_limit; ++c )
    out << "------------";
  out << endl;

  out << "      |    \\     |";
  for( size_t c = 0; c < col_limit; ++c )
    out << "     " << c+1 << "     |";
  out << endl;

  out << "      ------------";
  for( size_t c = 0; c < col_limit; ++c )
    out << "------------";
  out << endl;

  for( size_t l = 0; l < row_limit; ++l )
  {
    out << "      |     " << l+1 << "    |";
    for( size_t c = 0; c < col_limit; ++c )
      out << fp(test.get_lodpal_result().var_cov_matrix()(l,c), 10, 6) << " |";
    out << endl;

    out << "      ------------";
    for( size_t c = 0; c < col_limit; ++c )
      out << "------------";
    out << endl;
  }
  out << endl;
}

void
lodpal_test_diagfile::print_histogram(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  if( test.get_lodpal_result().last_error() == -1 )
  {
    out << "# Didn't maximized at this point!!" << endl << endl;
    return;
  }

  double lod_max = -numeric_limits<double>::infinity();
  double lod_min =  numeric_limits<double>::infinity();

  for( size_t i = 0; i < test.pairs_info().pairs_info().size(); ++i )
  {
    const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

    if( a_pair.removed )
      continue;

    if(    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
        && a_pair.removed_x )
      continue;

    lod_max = max(lod_max, a_pair.likelihood);
    lod_min = min(lod_min, a_pair.likelihood);
  }

  if( lod_max == lod_min )
  {
    out << "# Not a good data to draw histogram!!" << endl << endl;
    return;
  }

  out << "# Histogram of Individual Lod Score Contribution" << endl << endl;
  out << "    Maximum Lod Score = " << fp(lod_max, 7, 4) << endl
      << "    Minimum Lod Score = " << fp(lod_min, 7, 4) << endl;

  Histogram histo_lod(10, lod_min - 0.0001, lod_max + 0.0001);

  for( size_t i = 0; i < test.pairs_info().pairs_info().size(); ++i )
  {
    const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

    if( a_pair.removed )
      continue;

    if(    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
        && a_pair.removed_x )
      continue;

    histo_lod.add(a_pair.likelihood);
  }

  double bin_size  = histo_lod.binsize();
  //size_t max_bin   = histo_lod.max_bins();
  size_t max_count = histo_lod.maxFrequency().second;

  double points_per_star = 1;

  //double base_line = 48.;
  double base_line = (11 + my_max_ped_name + 2*my_max_ind_name) + 19;
  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    base_line += 11.;
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    base_line += 11.;

  if( max_count > base_line )
    points_per_star = ceil(max_count / base_line);

  out << "    Bin Size          = " << fp(bin_size, 7, 4) << endl << endl;

  out << "    Interval            Count (one * is equal up to "
      << points_per_star
      << " pair(s).)"
      << endl;

  print_single_line(test);

  double prev_bin = histo_lod.lower_bound();

  for( size_t i = 0; i < histo_lod.max_bins(); ++i )
  {
    out << fp(prev_bin, 7, 4) << "  to  " ;

    prev_bin += bin_size;

    out << fp(prev_bin, 7, 4) << "\t" << setw(4);
    out.setf(ios_base::right, ios_base::adjustfield);
    out << histo_lod.binCount(i);

    if( !histo_lod[i] )
      out << endl;
    else
    {
      out << " ";

      size_t count = (size_t) (ceil(histo_lod.binCount(i) / points_per_star));

      if( !count )
        out << "*";
      else
        for( size_t j = 0; j < count; ++j )
          out << "*";

      out << endl;
    }
  }

  print_single_line(test);

  out << "                Total : " << setw(4);
  out.setf(ios_base::right, ios_base::adjustfield);
  out << histo_lod.point_count()
      << endl << endl;
}

double
get_beta1(const ARP_base_analysis& test, size_t i, size_t j)
{
  const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

  lodpal_parameters::marker_const_iterator pi = test.parameters().marker_begin();

  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
  {
    if( !test.parameters().x_linkage_model().lambda1_equal )
    {
      if( test.pairs_info().is_mm_pair(a_pair.prior_x_ibd_index) )
        return pi->beta1mm.value();
      else if( test.pairs_info().is_mf_pair(a_pair.prior_x_ibd_index) )
        return pi->beta1mf.value();
      else
        return pi->beta1ff.value();
    }
    else
      return pi->beta1mm.value();
  }
  else
  {
    if( j == 0 )
      return pi->beta1.value();
    else if( j == 1 )
      return pi->beta1m.value();
    else if( j == 2 )
      return pi->beta1p.value();
    else
      return pi->beta1.value();
  }
}

double
get_delta1y(const ARP_base_analysis& test, size_t i, size_t j)
{
  const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

  double delta1y = 0.;

  for( size_t c = 0; c < a_pair.lodpal_cov.size(); ++c )
  {
    double delta1 = std::numeric_limits<double>::quiet_NaN();

    if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    {
      if( !test.parameters().x_linkage_model().lambda1_equal )
      {
        if( test.pairs_info().is_mm_pair(a_pair.prior_x_ibd_index) )
          delta1 = test.parameters().covariate_parameters(c).delta1mm.value();

        else if( test.pairs_info().is_mf_pair(a_pair.prior_x_ibd_index) )
          delta1 = test.parameters().covariate_parameters(c).delta1mf.value();

        else
          delta1 = test.parameters().covariate_parameters(c).delta1ff.value();
      }
      else
        delta1 = test.parameters().covariate_parameters(c).delta1mm.value();
    }
    else
    {
      if( j == 0 )
        delta1 = test.parameters().covariate_parameters(c).delta1.value();
      else if( j == 1 )
        delta1 = test.parameters().covariate_parameters(c).delta1m.value();
      else if( j == 2 )
        delta1 = test.parameters().covariate_parameters(c).delta1p.value();
    }

    double y = a_pair.lodpal_cov[c].ad_pair_value;

    if(    test.pairs_info().re_built_pairs_info()
        && c == test.pairs_info().re_built_covariate() )
      y = a_pair.lodpal_cov[c].re_ad_pair_value;

    delta1y += (delta1 * y);
  }

  return delta1y;
}

double
get_linkage_probability(const ARP_base_analysis& test, size_t i)
{
  const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

  const PALBASE::rel_pair& p = a_pair.lodpal_pair;

  double alpha = test.pairs_info().parameters().autosomal_model().alpha;

  lodpal_parameters::marker_const_iterator pi = test.parameters().marker_begin();

  double phiy_hat = 0.;

  if( !test.pairs_info().parameters().autosomal_model().parent_of_origin )
  {
    double f1 = p.prob_share(pi->marker, 1);
    double f2 = p.prob_share(pi->marker, 2);

    double beta1   = get_beta1(test, i, 0);
    double delta1y = get_delta1y(test, i, 0);

    double lambda1y = exp(beta1 + delta1y);

    if( lambda1y < 1. )
      lambda1y = 1.;

    double phiy = ((lambda1y - 1.) * (alpha + 3.)) / ((lambda1y - 1.) * (alpha + 3.) + 4.);
    double a    = 2. / (alpha + 3.);

    double z1y = phiy * a        + (1. - phiy)/2.;
    double z2y = phiy * (1. - a) + (1. - phiy)/4.;

    phiy_hat = phiy * ( (a * f1) / z1y + ((1. - a) * f2) / z2y );
  }
  else
  {
    double f1   = p.prob_share(pi->marker, 1);
    double f2   = p.prob_share(pi->marker, 2);
    double f1mp = p.prob_share(pi->marker, 3);
    double f1m  = (f1 + f1mp) / 2.0;
    double f1p  = (f1 - f1mp) / 2.0;

    double beta1m = 0.;
    double beta1p = 0.;

    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      beta1m = get_beta1(test, i, 1);
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      beta1p = get_beta1(test, i, 2);

    double delta1my = 0.;
    double delta1py = 0.;

    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::maternal )
      delta1my = get_delta1y(test, i, 1);
    if( test.pairs_info().parameters().autosomal_model().fixed != autosomal_model_type::paternal )
      delta1py = get_delta1y(test, i, 2);


    double lambda1my = exp(beta1m + delta1my);
    double lambda1py = exp(beta1p + delta1py);

    if( lambda1my < 1. )
      lambda1my = 1.;

    if( lambda1py < 1. )
      lambda1py = 1.;

    double lambda1y = (lambda1my + lambda1py) / 2.0;

    double phiy = ((lambda1y - 1.) * (alpha + 3.)) / ((lambda1y - 1.) * (alpha + 3.) + 4.);

    double b    = (lambda1my - 1.) / ((alpha + 3.) * (lambda1y - 1.));
    double c    = (lambda1py - 1.) / ((alpha + 3.) * (lambda1y - 1.));

    if( lambda1y == 1. )
      b = c = 1. / (alpha + 3.);

    double bc   = 2. / (alpha + 3.);

    double z1my = phiy * b         + (1. - phiy)/4.;
    double z1py = phiy * c         + (1. - phiy)/4.;
    double z2y  = phiy * (1. - bc) + (1. - phiy)/4.;

    phiy_hat = phiy * ( (b * f1m) / z1my + (c * f1p) / z1py + ((1. - bc) * f2) / z2y );

    if( phiy_hat <= -0. )
      phiy_hat = 0.;

#if 0
  cout << "b1m  = " << beta1m   << ", b1p  = " << beta1p   << "  ";
  cout << "d1my = " << delta1my << ", d1py = " << delta1py << endl;
  cout << "lambda1y = " << lambda1y << "  "
       << "phiy = " << phiy << "  "
       << "b    = " << b    << "  "
       << "c    = " << c    << "  "
       << "bc   = " << bc   << "  "
       << "phiy_hat = " << phiy_hat << endl;
#endif
  }

  return phiy_hat;
}

void
lodpal_test_diagfile::print_individual_lod_table(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  out << "# Individual Lod Score Contribution" << endl;

  if( test.pairs_info().removed_biggest_pair() )
    out << endl << "    *R* indicates the removed pair from the analysis at this location.";
  if(    test.pairs_info().mm_invalid_pair_count() > 0
      || test.pairs_info().mf_invalid_pair_count() > 0
      || test.pairs_info().ff_invalid_pair_count() > 0 )
    out << endl << "    *I* indicates the invalid or unused pair for the analysis at this location.";

  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    print_type_regend(test);
  else
    out << endl;

  out << endl;
  out.setf(ios_base::left, ios_base::adjustfield);
  out << setw(my_max_ped_name) << " " << " "
      << setw(my_max_ind_name) << " " << " "
      << setw(my_max_ind_name) << " " << " ";
  out << "PAIR  ";
  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    out << "PAIR  ";
//  out << "                          ";
  out << "              ";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << "           ";
  out << "            ";

  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    out << " Covariate ";

  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    out << " Weight    ";
  
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << "  LOD SCORE       LINKAGE" << endl;
  else
    out << "    LOD SCORE" << endl;

  out << setw(my_max_ped_name) << "FAMID" << " "
      << setw(my_max_ind_name) << "ID1"   << " "
      << setw(my_max_ind_name) << "ID2"   << " ";
  out << "TYPE  ";
  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    out << "SEX   ";

//  out << "    F0         F2         ";
  out << "    F0        ";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << " F1m-p     ";
  out << " F2         ";

  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
  {
    string c_name = test.parameters().covariate_parameters(ci).name(test.relative_pairs());
    out << " " << setw(10) << c_name.substr(0,10);
  }
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
  {
    string w_name = test.parameters().weight_parameter().name(test.relative_pairs());
    out << " " << setw(10) << w_name.substr(0,10);
  }
  
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << "  CONTRIBUTION    PROBABILITY" << endl;
  else
    out << "    CONTRIBUTION" << endl;

  for( size_t i = 0; i < my_max_ped_name; ++i )
    out << "-";
  out << " ";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << " ";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << " ";
  out << "----- ";
  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    out << "------";
//  out << "    ---------- ---------- ";
  out << "    ----------";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << " ----------";
  out << " ---------- ";

  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    out << " ----------";
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    out << " ----------";
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << "  -------------   -----------";
  else
  {
    out << "  ----------------";
    if( !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
      out << "------";
  }
  out << endl;

  double lod_sum  = 0.;
  double prob_sum = 0.;

  for( size_t i = 0; i < test.pairs_info().pairs_info().size(); ++i )
  {
    const lodpal_pairs::lodpal_pair_info& a_pair = test.pairs_info().pairs_info()[i];

    const PALBASE::rel_pair& p = a_pair.lodpal_pair;

    string rel1 = p.rels().pair.first->name();
    string rel2 = p.rels().pair.second->name();
    string ped1 = p.rels().pair.first->pedigree()->name(); 
    string ped2 = p.rels().pair.second->pedigree()->name();

    assert(ped1 == ped2);

    if(    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
        && a_pair.removed_x )
      continue;

    out << setw(my_max_ped_name) << ped1 << " " 
        << setw(my_max_ind_name) << rel1 << " " 
        << setw(my_max_ind_name) << rel2 << " ";

    string type = "other";
    if( p.is_fsib_pair() )
      type = "sib";
    else if( p.is_hsib_pair() )
      type = "h.sib";

    out << setw(5) << type << " ";

    if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
      if( a_pair.prior_x_ibd_index < test.pairs_info().prior_x_ibds().size() )
        out << setw(6) << test.pairs_info().prior_x_ibds()[a_pair.prior_x_ibd_index].subtype;
      else
        out << "      ";

    if( a_pair.removed )
      out << "*R* ";
    else if(    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
             && a_pair.removed_x )
      out << "*I* ";
    else
    {
      out << "    ";
      lod_sum += a_pair.likelihood;
    }

//    out << fp(a_pair.f0, 10, 7) << " "
//        << fp(a_pair.f2, 10, 7) << " ";

    out << fp(a_pair.f0, 10, 7) << " ";
    if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
      out << fp(a_pair.f1mp, 10, 7) << " ";
    out << fp(a_pair.f2, 10, 7) << " ";

    for( size_t c = 0; c < a_pair.lodpal_cov.size(); ++c )
      if(    test.parameters().covariate_parameters(c).adjust != covariate_type::dsp
          && test.pairs_info().re_built_pairs_info()
          && c == test.pairs_info().re_built_covariate() )
        out << " " << fp(a_pair.lodpal_cov[c].re_ad_pair_value, 10, 7);
      else
        out << " " << fp(a_pair.lodpal_cov[c].ad_pair_value, 10, 7);
    if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
      out << " " << fp(a_pair.lodpal_weight, 10, 7);

    if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
        && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
                 && test.parameters().x_linkage_model().lambda2_fixed )
            || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
                 && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    {
      if( test.get_lodpal_result().last_error() != -1 )
        out << "  " << fp(a_pair.likelihood, 13, 10);
      else
        out << "  " << fp(std::numeric_limits<double>::quiet_NaN(), 13, 10);

      if( type != "sib" )
      {
        out << "   " << fp(std::numeric_limits<double>::quiet_NaN(), 10, 7);
      }
      else
      {
        double linkage_prob = get_linkage_probability(test, i);
        prob_sum += linkage_prob;

        out << "   " << fp(linkage_prob, 10, 7);
      }
    }
    else
    {
      if( test.get_lodpal_result().last_error() != -1 )
        out << "    " << fp(a_pair.likelihood, 13, 10);
      else
        out << "    " << fp(std::numeric_limits<double>::quiet_NaN(), 13, 10);
    }

    out << endl;
  }

  for( size_t i = 0; i < my_max_ped_name; ++i )
    out << "-";
  out << " ";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << " ";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << " ";
  out << "----- ";
  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    out << "------";
//  out << "    ---------- ---------- ";
  out << "    ----------";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << " ----------";
  out << " ---------- ";

  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    out << " ----------";
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    out << " ----------";
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << "  -------------   -----------";
  else
  {
    out << "  ----------------";
    if( !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
      out << "------";
  }

  out << endl;

  out.setf(ios_base::left,ios_base::adjustfield);
  out << "Total Pair Count = " << setw(5) << test.pairs_info().pairs_info().size();

  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << setw(my_max_ped_name+my_max_ind_name*2-5) << " ";
  else
    out << setw(my_max_ped_name+my_max_ind_name*2-3) << " ";
  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    out << "           ";
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    out << "           ";

  if( test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker()) )
    out << setw(6) << " ";

  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << setw(11) << " ";

  out << "Total Lod Score = ";
  if( test.get_lodpal_result().last_error() != -1 )
    out << fp(lod_sum, 13, 10);
  else
    out << fp(std::numeric_limits<double>::quiet_NaN(), 13, 10);

  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && (   (    test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().x_linkage_model().lambda2_fixed )
          || (   !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
               && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter ) ) )
    out << "   " << fp(prob_sum, 10, 7);

  out << endl;
}

void
lodpal_test_diagfile::print_type_regend(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  out << endl;
  out << "    ------------------------------------------------------------------------------------" << endl;
  out << "     Pair Type      Pair Sex  Definition" << endl;
  out << "    ------------------------------------------------------------------------------------" << endl;
  out << "     sib" << endl;
  if( test.parameters().use_mm_pair() )
    out << "                    M,M       brother:brother" << endl;
  if( test.parameters().use_mf_pair() )
    out << "                    M,F       brother:sister" << endl;
  if( test.parameters().use_ff_pair() )
    out << "                    F,F       sister:sister" << endl;
  out << "     h.sib" << endl;
  if( test.parameters().use_mm_pair() )
    out << "                    M,F,M     maternal half-brother:half-brother" << endl;
  if( test.parameters().use_mf_pair() )
    out << "                    M,F,F     maternal half-brother:half-sister" << endl;
  if( test.parameters().use_ff_pair() )
    out << "                    F,F,F     maternal half-sister:half-sister" << endl;
  out << "     grandparental" << endl;
  if( test.parameters().use_mm_pair() )
    out << "                    MFM       grandfather-through-mother:grandson" << endl;
  if( test.parameters().use_mf_pair() )
  {
    out << "                    MFF       grandfather-through-mother:granddaughter" << endl;
    out << "                    FFM       grandmother-through-mother:grandson" << endl;
  }
  if( test.parameters().use_ff_pair() )
    out << "                    FFF       grandmother-through-mother:granddaughter" << endl;
  out << "     avuncular" << endl;
  if( test.parameters().use_mm_pair() )
    out << "                    M,FM      uncle-through-mother:nephew" << endl;
  if( test.parameters().use_mf_pair() )
  {
    out << "                    M,FF      uncle-through-mother:niece" << endl;
    out << "                    M,MF      uncle-through-father:niece" << endl;
    out << "                    F,FM      aunt-through-mother:nephew" << endl;
  }
  if( test.parameters().use_ff_pair() )
  {
    out << "                    F,FF      aunt-through-mother:niece" << endl;
    out << "                    F,MF      aunt-through-father:niece" << endl;
  }
  out << "     cousins" << endl;
  if( test.parameters().use_mm_pair() )
    out << "                    MF,FM     male-cousin-through-mother:male-cousin-through-mother" << endl;
  if( test.parameters().use_mf_pair() )
  {
    out << "                    MF,FF     male-cousin-through-mother:female-cousin-through-mother" << endl;
    out << "                    FM,FM     female-cousin-through-father:male-cousin-through-mother" << endl;
  }
  if( test.parameters().use_ff_pair() )
  {
    out << "                    FM,MF     female-cousin-through-father:female-cousin-through-father" << endl;
    out << "                    FM,FF     female-cousin-through-father:female-cousin-through-mother" << endl;
    out << "                    FF,FF     female-cousin-through-mother:female-cousin-through-mother" << endl;
  }
  out << "    ------------------------------------------------------------------------------------" << endl;
}

void
lodpal_test_diagfile::print_double_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i < my_max_ped_name; ++i )
    out << "=";
  out << "=";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "=";
  out << "=";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "=";
  out << "=";
  out << "================================";
  for( size_t ci = 0; ci < test.pairs_info().parameters().covariate_count(); ++ci )
    out << "===========";
  if( test.pairs_info().parameters().weight_parameter().weight != (size_t)-1 )
    out << "===========";
  out << "========================";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << "===========";
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
      && test.parameters().x_linkage_model().lambda2_fixed )
    out << "===========";
  else if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
           && !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
           && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "=====";
  out << endl;
}

void
lodpal_test_diagfile::print_single_line(const ARP_base_analysis& test)
{
  ostream& out = output_stream();

  for( size_t i = 0; i < my_max_ped_name; ++i )
    out << "-";
  out << "-";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << "-";
  for( size_t i = 0; i < my_max_ind_name; ++i )
    out << "-";
  out << "-";
  out << "--------------------------------";
  if( test.pairs_info().parameters().autosomal_model().parent_of_origin )
    out << "-----------";

  for( size_t ci = 0; ci < test.parameters().covariate_count(); ++ci )
    out << "-----------";
  if( test.parameters().weight_parameter().weight != (size_t)-1 )
    out << "-----------";
  out << "------------------------";
  if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
      && test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
      && test.parameters().x_linkage_model().lambda2_fixed )
    out << "-----------";
  else if(    test.pairs_info().parameters().trait_parameters(0).pair_select == trait_parameter::conaff
           && !test.relative_pairs().is_x_linked(test.pairs_info().parameters().diagnostic_marker())
           && test.parameters().autosomal_model().model == autosomal_model_type::one_parameter )
    out << "-----";
  out << endl;
}

/*
void
lodpal_test_diagfile::print_bad_sib_pair(const ARP_base_analysis&& test)
{
  ostream& out = output_stream();
  bad_sib_pair_type::type_const_iterator ti;
  bad_sib_pair_type::value_const_iterator vi;

  out << "===============================================" << endl;
  out << "Location/Marker  FamID     ID/Sib1     ID/Sib2   " << endl;
  for( ti = test.bad_sib_pair().type_begin(); ti != test.bad_sib_pair().type_end(); ++ti )
  {
    out << "-----------------------------------------------" << endl;
    out << setw(16) << ti->first.name(test.relative_pairs()) << "  ";
    for( vi = ti->second.begin(); vi != ti->second.end(); ++vi )
    {
      const SAGE::rel_pair& p = *vi;
      string rel1 = p.rels().pair.first->name();
      string rel2 = p.rels().pair.second->name();
      string ped1  = p.rels().pair.first->pedigree()->name(); 
      string ped2  = p.rels().pair.second->pedigree()->name();

      assert(ped1 == ped2);

      if( vi != ti->second.begin() )
        out << setw(16) << " " << "  ";

      out << setw(10) << ped1 << "  " 
          << setw(10) << rel1 << "  " 
          << setw(10) << rel2 << "  "
          << endl;
    }
  }
  out << "===============================================" << endl;
  out << endl;
}
*/

} // end of namespace LODPAL
} // end of namespace SAGE

