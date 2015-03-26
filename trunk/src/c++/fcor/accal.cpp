//****************************************************************************
//* File:      accal.cpp                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Oct 99 *
//*            1.0 Maintype computation added                     yjs May 01 *
//*            1.1 Changed comp. method for std.err calculation.  yjs Jul 01 *
//*                                                                          *
//* Notes:     This source file calculates the asymptotic covariance of      *
//*            multipedigree correlations.                                   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/accal.h"

namespace SAGE {
namespace FCOR {

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of ACCal
// ---------------------------------------------------------------------------

ACCal::ACCal(const CorrelationCal* c)
{
  my_correlation = c;
}

void
ACCal::compute_covariance(const pairset_by_pedigree_type&  rel_type1,
                          const pairset_by_pedigree_type&  rel_type2,
                          const vector<CorrelationInfo>&   cor1,
                          const vector<CorrelationInfo>&   cor2,
                          const weight_matrix_by_pedigree& optimal_w1,
                          const weight_matrix_by_pedigree& optimal_w2,
                          Matrix2D< vector<double> >&      ac_matrix,
                          Matrix2D< vector<bool> >&        re_matrix,
                          Matrix2D< vector<size_t> >&      lp_matrix) const
{
  if( my_correlation == NULL )
    return;

  size_t trait_count = my_correlation->get_parser()->get_trait_count();
  size_t wcnt = WEIGHT_COUNT;

  if( optimal_w1.size() )
    wcnt = WEIGHT_COUNT + 1;

  size_t matrix_size = trait_count * trait_count;

  vector<double> ac(wcnt, std::numeric_limits<double>::quiet_NaN());
  vector<bool>   re(wcnt, false);
  vector<size_t> lp(wcnt, size_t(-1));

  ac_matrix.resize(matrix_size, matrix_size, ac);
  re_matrix.resize(matrix_size, matrix_size, re);
  lp_matrix.resize(matrix_size, matrix_size, lp);

  if( rel_type1 != rel_type2 )
  {
    Matrix2D< vector<double> >  ac_matrix1(matrix_size, matrix_size, ac);
    Matrix2D< vector<double> >  ac_matrix2(matrix_size, matrix_size, ac);

    Matrix2D< vector<bool> >    re_matrix1(matrix_size, matrix_size, re);
    Matrix2D< vector<bool> >    re_matrix2(matrix_size, matrix_size, re);

    Matrix2D< vector<size_t> >  lp_matrix1(matrix_size, matrix_size, lp);
    Matrix2D< vector<size_t> >  lp_matrix2(matrix_size, matrix_size, lp);

    SubACCalM ac_m(my_correlation);

    ac_m.compute_asymptotic_covariance(rel_type1,  rel_type2,
                                       cor1,       cor2,
                                       optimal_w1, optimal_w2,
                                       ac_matrix1, re_matrix1, lp_matrix1);
    ac_m.compute_asymptotic_covariance(rel_type2,  rel_type1,
                                       cor2,       cor1,
                                       optimal_w2, optimal_w1,
                                       ac_matrix2, re_matrix2, lp_matrix2);

    for( size_t i = 0; i < ac_matrix.rows(); ++i )
      for( size_t j = 0; j < ac_matrix.cols(); ++j )
        for( size_t w = 0; w < ac_matrix(i,j).size() ; ++w )
        {
          double m = std::numeric_limits<double>::quiet_NaN();

          if( !(isnan(ac_matrix1(i,j)[w]) || isnan(ac_matrix2(i,j)[w])) )
          {
            m = (ac_matrix1(i,j)[w] + ac_matrix2(j,i)[w]) / 2.;

            if( fabs(m) < 1.0e-11 )
              m = 0.0;
          }

          ac_matrix(i,j)[w] = m;

          if( re_matrix1(i,j)[w] || re_matrix2(j,i)[w] )
            re_matrix(i,j)[w] = true;

          lp_matrix(i,j)[w] = std::min(lp_matrix1(i,j)[w], lp_matrix2(j,i)[w]);
        }
  }
  else
  {
    SubACCalM ac_m(my_correlation);

    ac_m.compute_asymptotic_covariance(rel_type1,  rel_type1,
                                       cor1,       cor1,
                                       optimal_w1, optimal_w1,
                                       ac_matrix,  re_matrix,  lp_matrix);

    for( size_t i = 0; i < ac_matrix.rows(); ++i )
      for( size_t j = 0; j <= i; ++j )
        for( size_t w = 0; w < ac_matrix(i,j).size() ; ++w )
        {
          double m = std::numeric_limits<double>::quiet_NaN();

          if( !(isnan(ac_matrix(i,j)[w]) || isnan(ac_matrix(j,i)[w])) )
          {
            m = (ac_matrix(i,j)[w] + ac_matrix(j,i)[w]) / 2.0;

            if( fabs(m) < 1.0e-11 )
              m = 0.0;
          }

          ac_matrix(i,j)[w] = ac_matrix(j,i)[w] = m;

          if( re_matrix(i,j)[w] || re_matrix(j,i)[w] )
            re_matrix(i,j)[w] = re_matrix(j,i)[w] = true;
        }
  }

  return;
}

void
ACCal::compute_variance(const pairset_by_pedigree_type&  rel_type1,
                        const vector<CorrelationInfo>&   cor1,
                        const weight_matrix_by_pedigree& optimal_w1,
                        Matrix2D< vector<double> >&      av_matrix,
                        Matrix2D< vector<bool> >&        re_matrix,
                        Matrix2D< vector<size_t> >&      lp_matrix) const
{
  if( my_correlation == NULL )
    return;

  size_t trait_count = my_correlation->get_parser()->get_trait_count();
  size_t wcnt = WEIGHT_COUNT;

  if( optimal_w1.size() )
    wcnt = WEIGHT_COUNT + 1;

  size_t matrix_size = trait_count * trait_count;

  vector<double> ac(wcnt, std::numeric_limits<double>::quiet_NaN());
  vector<bool>   re(wcnt, false);
  vector<size_t> lp(wcnt, size_t(-1));

  av_matrix.resize(matrix_size, matrix_size, ac);
  re_matrix.resize(matrix_size, matrix_size, re);
  lp_matrix.resize(matrix_size, matrix_size, lp);

  SubAVCalM av_m(my_correlation);
  av_m.compute_asymptotic_variance(rel_type1, cor1, optimal_w1, av_matrix, re_matrix, lp_matrix);

  return;
}

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of StdErrCal
// ---------------------------------------------------------------------------

StdErrCal::StdErrCal(const CorrelationCal& c)
{
  my_correlation = &c;
}

void
StdErrCal::compute_standard_errors(pairset_result_vector& results)
{
  if( my_correlation == NULL )
    return;

  const pairset_vector&      pairsets = *my_correlation->get_pairset();
  const pairset_info_vector& pinfos   = *my_correlation->get_pairset_info();
  const corinfo_vector&      corinfos = my_correlation->get_corinfo();

  weight_type w           = my_correlation->get_parser()->get_analysis_options().class_weight;
  size_t generation_limit = my_correlation->get_parser()->get_analysis_options().generation_limit;
  size_t trait_count      = my_correlation->get_parser()->get_trait_count();
  size_t matrix_size      = trait_count * trait_count;

  for( size_t r = 0; r < pairsets.size(); ++r )
  {
    if( is_invalid_pair_type(pinfos[r], generation_limit, trait_count) )
      continue;

    // new asymptotic covariance method
    //
    Matrix2D< vector<double> >  ac_matrix;
    Matrix2D< vector<bool> >    re_matrix;
    Matrix2D< vector<size_t> >  lp_matrix;

    weight_matrix_by_pedigree weight;

    ACCal cal(my_correlation);

    cal.compute_covariance(pairsets[r], pairsets[r],
                           corinfos[r], corinfos[r],
                           weight,      weight,
                           ac_matrix,   re_matrix,   lp_matrix);

    set_standard_error(results[r], ac_matrix, re_matrix, lp_matrix);

#if 0
    print_vec_matrix(ac_matrix, cout, pinfos[r].gname);
#endif

    if( w == WEIGHT_COUNT )
    {
      optimal_weight_finder owf;
      owf.estimate_variance_by_quadradic(results[r]);
    }

    if( pinfos[r].type == OTHER && my_correlation->get_parser()->get_trait_count() > 1 )
    {
      Matrix2D<double> new_ac_matrix;
      new_ac_matrix.resize(matrix_size, matrix_size, std::numeric_limits<double>::quiet_NaN());

      if( w == WEIGHT_COUNT )
      {
        optimal_weight_finder owf;
        owf.estimate_weighted_covariance(trait_count, results[r], results[r], ac_matrix, new_ac_matrix);
      }
      else
      {
        for( size_t t1 = 0; t1 < ac_matrix.rows(); ++t1 )
          for( size_t t2 = 0; t2 < ac_matrix.cols(); ++t2 )
            if( !isnan(ac_matrix(t1, t2)[w]) )
            {
              if( fabs(ac_matrix(t1, t2)[w]) < 1.0e-11 )
                new_ac_matrix(t1, t2) = 0.;
              else
                new_ac_matrix(t1, t2) = ac_matrix(t1, t2)[w];
            }
      }
      
#if 0
      print_matrix(new_ac_matrix, cout, pinfos[r].gname);
#endif

      compute_pooled_cross_correlation(results[r], new_ac_matrix);
    }
  }

  return;
}

void
StdErrCal::set_standard_error(pairset_result& pr, const Matrix2D< vector<double> >  ac_matrix,
                                                  const Matrix2D< vector<bool> >    re_matrix,
                                                  const Matrix2D< vector<size_t> >  lp_matrix) const
{
  size_t trait_count = my_correlation->get_parser()->get_trait_count();

  size_t t = 0;
  for( size_t t1 = 0; t1 < trait_count; ++t1 )
    for( size_t t2 = 0; t2 < trait_count; ++t2, ++t )
      for( size_t w = 0; w < ac_matrix(t, t).size(); ++w )
      {
        if( !isnan(ac_matrix(t, t)[w]) )
        {
          if( fabs(ac_matrix(t, t)[w]) < 1.0e-11 )
            pr.std_err(t1, t2).standard_error[w] = 0.;
          else if( ac_matrix(t, t)[w] >= 1.0e-11 )
            pr.std_err(t1, t2).standard_error[w] = ac_matrix(t, t)[w];

          pr.std_err(t1, t2).replaced         = re_matrix(t, t)[w];
          pr.std_err(t1, t2).least_pair_count = lp_matrix(t, t)[w];
        }
      }

  return;
}

void
StdErrCal::compute_pooled_cross_correlation(pairset_result&           pr, 
                                            const Matrix2D< double >& ac_matrix) const
{
  size_t trait_count = my_correlation->get_parser()->get_trait_count();
  size_t df          = (trait_count*(trait_count-1)) / 2;
  weight_type w      = my_correlation->get_parser()->get_analysis_options().class_weight;

  for( size_t t1 = 0; t1 < trait_count; ++t1 )
  {
    for( size_t t2 = 0; t2 < t1; ++t2)
    {
      double x = pr.corr(t1, t2).correlation[w];
      double y = pr.corr(t2, t1).correlation[w];

      double Vx   = ac_matrix(t1*trait_count + t2, t1*trait_count + t2);
      double Vy   = ac_matrix(t2*trait_count + t1, t2*trait_count + t1);
      double CVxy = ac_matrix(t1*trait_count + t2, t2*trait_count + t1);
      double ow   = (Vy - CVxy) / (Vx + Vy - 2.0*CVxy);

      double pooled_cor = ow*x + (1.0-ow)*y;
      double var        = ow*ow*Vx + (1.0-ow)*(1.0-ow)*Vy + 2.0*ow*(1.0-ow)*CVxy;
      double chiq       = ((x-y)*(x-y)) / (Vx + Vy - 2.0*CVxy);
#if 0
  cout << "t1 = " << t1 << ", t2 = " << t2 << ", x = " << x << ", y = " << y
       << ", Vx = " << Vx << ", Vy = " << Vy << ", CVxy = " << CVxy << ", w = " << ow
       << ", pooled_cor = " << pooled_cor << ", var = " << var << ", chisq = " << chiq << endl;
#endif
      pr.pooled_cross_corr(t1, t2).correlation = pooled_cor;
      pr.pooled_cross_corr(t1, t2).variance    = var;
      pr.pooled_cross_corr(t1, t2).chi_square  = chiq;
    }
  }

  if( trait_count > 2 )
  {
    Matrix2D<double> x_y(df, 1, std::numeric_limits<double>::quiet_NaN());
    Matrix2D<double> Vx_y(df, df, std::numeric_limits<double>::quiet_NaN());

    for( size_t t1_1 = 0, r = 0; t1_1 < trait_count; ++t1_1 )
    {
      for( size_t t2_1 = 0; t2_1 < t1_1; ++t2_1, ++r)
      {
        double x1 = pr.corr(t1_1, t2_1).correlation[w];
        double y1 = pr.corr(t2_1, t1_1).correlation[w];

        x_y(r, 0) = x1 - y1;

        for( size_t t1_2 = 0, c = 0; t1_2 < trait_count; ++t1_2 )
        {
          for( size_t t2_2 = 0; t2_2 < t1_2; ++t2_2, ++c)
          {
            double CVx1x2 = ac_matrix(t1_1*trait_count + t2_1, t1_2*trait_count + t2_2);
            double CVx1y2 = ac_matrix(t1_1*trait_count + t2_1, t2_2*trait_count + t1_2);
            double CVy1x2 = ac_matrix(t2_1*trait_count + t1_1, t1_2*trait_count + t2_2);
            double CVy1y2 = ac_matrix(t2_1*trait_count + t1_1, t2_2*trait_count + t1_2);

            Vx_y(r, c) = CVx1x2 - CVx1y2 - CVy1x2 + CVy1y2;
          }
        }
      }
    }

#if 0
    print_matrix(x_y, cout, "x-y");
    print_matrix(Vx_y, cout, "Vx-y");
#endif

    Matrix2D<double> x_y_transpose;
    Transpose(x_y, x_y_transpose);

    Matrix2D<double> Vx_y_inverse;
    Inverse(Vx_y, Vx_y_inverse);

#if 0
    print_matrix(Vx_y_inverse, cout, "(Vx-y) inverse");
#endif

    Matrix2D<double> temp1, chi_val;
    Multiply(x_y_transpose, Vx_y_inverse, temp1);
    Multiply(temp1, x_y, chi_val);

    if( chi_val.rows() && chi_val.cols() )
    {
      double chi_square = chi_val(0,0);

#if 0
  cout << "chisq = " << chi_square << endl;
#endif
      if( !SAGE::isnan(chi_square) && chi_square > std::numeric_limits<double>::epsilon() )
      {
        pr.pooled_cross_corr_chi_square = chi_square;
        pr.pooled_cross_corr_p_value    = 1.0 - chi_square_cdf(df, chi_square);
      }
    }
  }
  else
  {
    double chi_square = pr.pooled_cross_corr(1, 0).chi_square;
    pr.pooled_cross_corr_chi_square = chi_square;
    pr.pooled_cross_corr_p_value    = 1.0 - chi_square_cdf(df, chi_square);

#if 0
  cout << "chisq = " << chi_square << endl;
#endif
  }

  return;
}

void
StdErrCal::compute_standard_errors(const weight_matrix_vector&  weights,
                                         pairset_result_vector& results)
{
  if( my_correlation == NULL )
    return;

  const pairset_vector&      pairsets = *my_correlation->get_pairset();
  const pairset_info_vector& pinfos   = *my_correlation->get_pairset_info();
  const corinfo_vector&      corinfos = my_correlation->get_corinfo();

  size_t generation_limit = my_correlation->get_parser()->get_analysis_options().generation_limit;
  size_t trait_count      = my_correlation->get_parser()->get_trait_count();

  for( size_t r = 0; r < pairsets.size(); ++r )
  {
    if( is_invalid_pair_type(pinfos[r], generation_limit, trait_count) )
      continue;

    Matrix2D< vector<double> >  ac_matrix;
    Matrix2D< vector<bool> >    re_matrix;
    Matrix2D< vector<size_t> >  lp_matrix;

    ACCal cal(my_correlation);

    cal.compute_covariance(pairsets[r], pairsets[r],
                           corinfos[r], corinfos[r],
                           weights[r],  weights[r],
                           ac_matrix,   re_matrix,   lp_matrix);

    //set_KE_standard_error(results[r], ac_matrix, re_matrix, lp_matrix);
    size_t t = 0;
    for( size_t t1 = 0; t1 < trait_count; ++t1 )
      for( size_t t2 = 0; t2 < trait_count; ++t2, ++t )
      {
        if( !isnan(ac_matrix(t, t)[0]) )
        {
          if( fabs(ac_matrix(t, t)[0]) < 1.0e-11 )
            results[r].std_err(t1, t2).standard_error[WEIGHT_COUNT + 1] = 0.;
          else if( ac_matrix(t, t)[0] >= 1.0e-11 )
            results[r].std_err(t1, t2).standard_error[WEIGHT_COUNT + 1] = ac_matrix(t, t)[0];

          results[r].std_err(t1, t2).replaced         = re_matrix(t, t)[0];
          results[r].std_err(t1, t2).least_pair_count = lp_matrix(t, t)[0];
        }
      }
#if 0
    print_vec_matrix(ac_matrix, cout, pinfos[r].gname);
#endif
  }

  return;
}

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of VarCovCal
// ---------------------------------------------------------------------------

VarCovCal::VarCovCal(const CorrelationCal& sc, const CorrelationCal& mc)
{
  my_sub_correlation  = &sc;
  my_main_correlation = &mc;
}

void
VarCovCal::compute_variance_covariances(const pairset_result_vector& sub_results,
                                        const pairset_result_vector& main_results,
                                              var_cov_result_vector& vc_results)
{
  if( my_sub_correlation == NULL || my_main_correlation == NULL )
    return;

  size_t trait_count = my_sub_correlation->get_parser()->get_trait_count();
  size_t matrix_size = trait_count * trait_count;

  const vector<var_cov_param>& vc = my_sub_correlation->get_parser()->get_analysis_options().var_covs;

  const pairset_info_vector& sinfos = *my_sub_correlation->get_pairset_info();
  const pairset_info_vector& minfos = *my_main_correlation->get_pairset_info();

  my_var_covs.resize(vc.size());
  vc_results.resize(vc.size());

  for( size_t i = 0; i < vc.size(); ++i )
  {
    my_var_covs[i].resize(matrix_size, matrix_size, var_cov_result());
    vc_results[i].resize(matrix_size, matrix_size, var_cov_result());

    string reltype1 = vc[i].correlations.first;
    string reltype2 = vc[i].correlations.second;

    size_t m1 = find_reltype(minfos, reltype1);
    size_t m2 = m1;
    if( reltype1 != reltype2 )
      m2 = find_reltype(minfos, reltype2);

    if( m1 >= minfos.size() )
    {
      size_t s1 = find_reltype(sinfos, reltype1);
      size_t s2 = s1;
      if( reltype1 != reltype2 )
        s2 = find_reltype(sinfos, reltype2);

      if( s1 < sinfos.size() )
      {
        compute_variance_covariance(i, s1, s2, sub_results, my_sub_correlation);
      }
      else
      {
        string err_msg = "Relative pair type : " + reltype1;
        if( reltype1 != reltype2 )
          err_msg += ("/" + reltype2);

        sage_cerr << priority(warning) << err_msg           
                  << " does not exist! \n                     Skipping..." << endl; 
      }
    }
    else
    {
      compute_variance_covariance(i, m1, m2, main_results, my_main_correlation);
    }

    vc_results[i] = my_var_covs[i];
  }

  return;
}

void
VarCovCal::compute_variance_covariance(size_t i, size_t p1, size_t p2,
                                       const pairset_result_vector& results,
                                       const CorrelationCal*        corcal)
{
  const pairset_vector& pset = *corcal->get_pairset();
  const corinfo_vector& cors = corcal->get_corinfo();

  const pairset_by_pedigree_type& pairset1 = pset[p1];
  const pairset_by_pedigree_type& pairset2 = pset[p2];

  const corinfo_by_weight_type&   cor_type1 = cors[p1];
  const corinfo_by_weight_type&   cor_type2 = cors[p2];

  Matrix2D< vector<double> >  ac_matrix;
  Matrix2D< vector<bool> >    re_matrix;
  Matrix2D< vector<size_t> >  lp_matrix;

  weight_matrix_by_pedigree weight;

  ACCal cal(corcal);

  cal.compute_covariance(pairset1,  pairset2,
                         cor_type1, cor_type2,
                         weight,    weight,
                         ac_matrix, re_matrix, lp_matrix);

#if 0
  print_vec_matrix(ac_matrix, cout, "ac_matrix");
#endif

  if( my_sub_correlation->get_parser()->get_analysis_options().class_weight != WEIGHT_COUNT )
  {
    weight_type w = my_sub_correlation->get_parser()->get_analysis_options().class_weight;

    for( size_t r = 0; r < ac_matrix.rows(); ++r )
      for( size_t c = 0; c < ac_matrix.cols(); ++c )
        my_var_covs[i](r, c).var_cov = ac_matrix(r, c)[w];
  }
  else
  {
    const pairset_result& result1 = results[p1];
    const pairset_result& result2 = results[p2];

    size_t trait_count = my_sub_correlation->get_parser()->get_trait_count();
    size_t matrix_size = trait_count * trait_count;

    Matrix2D<double> optimal_ac_matrix;
    optimal_ac_matrix.resize(matrix_size, matrix_size, std::numeric_limits<double>::quiet_NaN());

    optimal_weight_finder owf;
    owf.estimate_weighted_covariance(trait_count, result1, result2, ac_matrix, optimal_ac_matrix);

#if 0
  print_matrix(optimal_ac_matrix, cout, "optimal_ac_matrix");
#endif

    for( size_t r = 0; r < optimal_ac_matrix.rows(); ++r )
      for( size_t c = 0; c < optimal_ac_matrix.cols(); ++c )
        my_var_covs[i](r, c).var_cov = optimal_ac_matrix(r, c);
  }

  for( size_t r = 0; r < re_matrix.rows(); ++r )
    for( size_t c = 0; c < re_matrix.cols(); ++c )
      for( size_t w = 0; w < re_matrix(r, c).size(); ++w )
      {
        my_var_covs[i](r, c).replaced         = re_matrix(r, c)[w];
        my_var_covs[i](r, c).least_pair_count = lp_matrix(r, c)[w];
      }

  return;
}

} // end of namespace FCOR
} // end of namespace SAGE
