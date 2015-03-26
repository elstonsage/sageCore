//****************************************************************************
//* File:      htest.cpp                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                    0.1 Revised to print var-cov matrix         Feb 01 01 * 
//*                                                                          *
//* Notes:     This source file calculates the chi-square & p-value of       *
//*            multipedigree correlations.                                   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/htest.h"

namespace SAGE {
namespace FCOR {


//---------------------------------------------------------------------------
// Out-of-Line Implementation of Htest
//---------------------------------------------------------------------------

Htest::Htest(const PairSetData& g, const CorrelationCal& c)
{
  my_pairsetdata = &g;
  my_correlation = &c;
}

Htest::~Htest()
{}

void
Htest::compute_homogeneity_test(const pairset_result_vector& sub_results,
                                      Htest_result_vector&   H_results)
{
  if( my_pairsetdata == NULL || my_correlation == NULL )
    return;

  // local declaration to simplify the code
  const pairset_info_vector& minfos = my_pairsetdata->get_maintype_info();
  const pairset_info_vector& sinfos = my_pairsetdata->get_subtype_info();

  const main_to_sub_map& g = my_pairsetdata->get_pairset_group();

  main_to_sub_type_const_iterator rel_group     = g.type_begin();
  main_to_sub_type_const_iterator rel_group_end = g.type_end();

  bool   individual_homog = my_correlation->get_parser()->get_analysis_options().individual_homog;
  size_t generation_limit = my_correlation->get_parser()->get_analysis_options().generation_limit;
  size_t trait_count      = my_correlation->get_parser()->get_trait_count();

  for( size_t r = 0; rel_group != rel_group_end; ++rel_group, ++r )
  {
    H_results.push_back(Htest_result());

    if( is_invalid_pair_type(minfos[r], generation_limit, trait_count) )
      continue;

    size_t subtype_count = rel_group->second.size();

    if( subtype_count < 2 )
      continue;

    bool valid_sub = true;
    for( size_t s = 0; s < rel_group->second.size(); ++s )
    {
      size_t s_type = rel_group->second[s];

      if( sinfos[s_type].total_pair_count < 3 )
      {
        valid_sub = false;
        break;
      }
    }

    if( !valid_sub )
      continue;

    // Get asymptotic covariance matrix.
    //
    vector< vector< Matrix2D<double> > > ac_matrix;
    bool replaced = compute_asymptotic_covariance(rel_group, sub_results, ac_matrix);

    if( !individual_homog )
    {
      size_t df = trait_count * trait_count * (subtype_count - 1);

      Matrix2D<internal_real_type> omega(df, df, std::numeric_limits<double>::quiet_NaN());

      compute_omega(ac_matrix, omega, trait_count);

      if( is_intraclass(rel_group->first) )
      {
        df = (trait_count * (trait_count+1) * (subtype_count-1)) / 2;

        Matrix2D<internal_real_type> new_omega(df, df, std::numeric_limits<double>::quiet_NaN());
        compute_omega_intra(omega, new_omega, df, trait_count);

        omega = new_omega;
      }

      Matrix2D<internal_real_type> y(df, 1, 0.0);

      compute_linear_correlation(rel_group, sub_results, y, (size_t)-1);

      double chi_square = compute_chi_square(y, omega);

      if(    SAGE::isnan(chi_square)
          || chi_square < std::numeric_limits<double>::epsilon() )
      {
        H_results[r].chi_square.push_back(std::numeric_limits<double>::quiet_NaN());
        H_results[r].p_value.push_back(std::numeric_limits<double>::quiet_NaN());
      }
      else
      {
        H_results[r].chi_square.push_back(chi_square);
        H_results[r].p_value.push_back(1.0 - chi_square_cdf(df, chi_square));
        H_results[r].replaced = replaced;
      }
    }
    else
    {
      for( size_t t = 0; t < trait_count; ++t )
      {
        vector< vector< Matrix2D<double> > > ac_matrix_ind;
        
        ac_matrix_ind.resize(0);
        ac_matrix_ind.resize(subtype_count);

        for( size_t i = 0; i < ac_matrix_ind.size(); ++i )
        {
          ac_matrix_ind[i].resize(0);
          ac_matrix_ind[i].resize(rel_group->second.size());

          for( size_t j = 0; j < ac_matrix_ind[i].size(); ++j )
          {
            ac_matrix_ind[i][j].resize(0, 0);
            ac_matrix_ind[i][j].resize(1, 1, std::numeric_limits<double>::quiet_NaN());

            ac_matrix_ind[i][j](0, 0) = ac_matrix[i][j]((trait_count+1)*t, (trait_count+1)*t);
#if 0
           cout << "ac_matrix_ind[" << i << "," << j << "] : " << endl;
           print_matrix(ac_matrix_ind[i][j], cout, "ac_matrix_ind");
#endif
          }
        }

        size_t df_ind = subtype_count - 1;

        Matrix2D<internal_real_type> y(df_ind, 1, 0.0);

        compute_linear_correlation(rel_group, sub_results, y, t);

        Matrix2D<internal_real_type> omega(df_ind, df_ind, std::numeric_limits<double>::quiet_NaN());

        compute_omega(ac_matrix_ind, omega, 1);

        double chi_square = compute_chi_square(y, omega);

        if(    SAGE::isnan(chi_square)
            || chi_square < std::numeric_limits<double>::epsilon() )
        {
          H_results[r].chi_square.push_back(std::numeric_limits<double>::quiet_NaN());
          H_results[r].p_value.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
          H_results[r].chi_square.push_back(chi_square);
          H_results[r].p_value.push_back(1.0 - chi_square_cdf(df_ind, chi_square));
          H_results[r].replaced = replaced;
        }
      }
    }
  }
}

bool
Htest::compute_asymptotic_covariance(main_to_sub_type_const_iterator       rel_group,
                                     const pairset_result_vector&          results,
                                     vector< vector< Matrix2D<double> > >& ac_matrix) const
{
  bool replaced = false;

  const pairset_vector&  pairsets = *my_correlation->get_pairset();
  const corinfo_vector&  corinfos = my_correlation->get_corinfo();

  size_t trait_count = my_correlation->get_parser()->get_trait_count();
  size_t matrix_size = trait_count * trait_count;

  ac_matrix.resize(rel_group->second.size());

#if 0
  cout << "rel group size = " << rel_group->second.size() << endl;
  cout << "pairsets size = " << pairsets.size() << endl;
  cout << "corinfos size = " << corinfos.size() << endl;
  cout << "results size = " << results.size() << endl;
#endif

  for( size_t i = 0; i < ac_matrix.size(); ++i )
  {
    ac_matrix[i].resize(rel_group->second.size());

    for( size_t j = 0; j < ac_matrix[i].size(); ++j )
    {
      ac_matrix[i][j].resize(matrix_size, matrix_size, std::numeric_limits<double>::quiet_NaN());

      size_t sub_type_i = rel_group->second[i];
      size_t sub_type_j = rel_group->second[j];

      //cout << "subtype i = " << sub_type_i << " "
      //     << "subtype j = " << sub_type_j << endl;

      // 1.
      const pairset_by_pedigree_type& pairset1 = pairsets[sub_type_i];
      const pairset_by_pedigree_type& pairset2 = pairsets[sub_type_j];

      const corinfo_by_weight_type&  cor_type1 = corinfos[sub_type_i];
      const corinfo_by_weight_type&  cor_type2 = corinfos[sub_type_j];

      Matrix2D< vector<double> >  ac_matrix_1;
      Matrix2D< vector<bool> >    re_matrix;
      Matrix2D< vector<size_t> >  lp_matrix;

      weight_matrix_by_pedigree weight;

      ACCal cal(my_correlation);

      cal.compute_covariance(pairset1,    pairset2,
                             cor_type1,   cor_type2,
                             weight,      weight,
                             ac_matrix_1, re_matrix, lp_matrix);

#if 0
      print_vec_matrix(ac_matrix_1, cout, "ac_matrix[i,j]");
#endif

      for( size_t r = 0; r < re_matrix.rows(); ++r )
        for( size_t c = 0; c < re_matrix.cols(); ++c )
          for( size_t w = 0; w < re_matrix(r, c).size(); ++w )
            if( re_matrix(r, c)[w] )
              replaced = true;

      // 2.
      weight_type w = my_correlation->get_parser()->get_analysis_options().class_weight;
      if( w != WEIGHT_COUNT )
      {
        for( size_t r = 0; r < ac_matrix_1.rows(); ++r )
          for( size_t c = 0; c < ac_matrix_1.cols(); ++c )
            ac_matrix[i][j](r, c) = ac_matrix_1(r, c)[w];
      }
      else
      {
        const pairset_result& result1 = results[sub_type_i];
        const pairset_result& result2 = results[sub_type_j];

        optimal_weight_finder owf;
        owf.estimate_weighted_covariance(trait_count, result1, result2, ac_matrix_1, ac_matrix[i][j]);
      }
    }
  }

  return replaced;
}

void
Htest::compute_linear_correlation(main_to_sub_type_const_iterator rel_group,
                                  const pairset_result_vector&    results,
                                  Matrix2D<internal_real_type>&   y,
                                  size_t                          t) const
{
  bool   conservative = my_correlation->get_parser()->get_analysis_options().conservative;
  size_t trait_count  = my_correlation->get_parser()->get_trait_count();
  size_t rel_type1    = rel_group->second[0];

  weight_type w = my_correlation->get_parser()->get_analysis_options().class_weight;

  const pairset_result& result1 = results[rel_type1];

  for( size_t s = 1, i = 0; s < rel_group->second.size(); ++s )
  {
    size_t rel_type2 = rel_group->second[s];

    const pairset_result& result2 = results[rel_type2];

    if( t < trait_count )
    {
      size_t t1 = t; //trait_count + t;
      size_t t2 = t;

      double cor_type1 = result1.corr(t1, t2).correlation[w];
      double cor_type2 = result2.corr(t1, t2).correlation[w];

      double type1 = 0.0;

      if( !SAGE::isnan(cor_type1) )
        type1 = cor_type1;
      else if( conservative )
        type1 = std::numeric_limits<double>::quiet_NaN();
      
      double type2 = 0.0;

      if( !SAGE::isnan(cor_type2) )
        type2 = cor_type2;
      else if( conservative )
        type2 = std::numeric_limits<double>::quiet_NaN();

      if( !SAGE::isnan(type1) && !SAGE::isnan(type2) )
        y(i, 0) = type1 - type2;
      else if( conservative )
        y(i, 0) = std::numeric_limits<double>::quiet_NaN();

      ++i;
    }
    else
    {
      for( size_t t1 = 0; t1 < trait_count; ++t1 )
      {

        size_t t2_start = 0;

        if(    is_intraclass(rel_group->first)
            || is_selfclass(rel_group->first) )
          t2_start = t1;

        for( size_t t2 = t2_start; t2 < trait_count; ++t2, ++i )
        {
          double cor_type1 = result1.corr(t1, t2).correlation[w];
          double cor_type2 = result2.corr(t1, t2).correlation[w];

          double type1 = 0.0;

          if( !SAGE::isnan(cor_type1) )
            type1 = cor_type1;
          else if( conservative )
            type1 = std::numeric_limits<double>::quiet_NaN();
          
          double type2 = 0.0;

          if( !SAGE::isnan(cor_type2) )
            type2 = cor_type2;
          else if( conservative )
            type2 = std::numeric_limits<double>::quiet_NaN();

          if( !SAGE::isnan(type1) && !SAGE::isnan(type2) )
            y(i, 0) = type1 - type2;
          else if( conservative )
            y(i, 0) = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }

#if 0
  print_matrix(y, cout, "y");
#endif
}

void
Htest::compute_omega(const vector< vector< Matrix2D<double> > >& ac_matrix,
                     Matrix2D<internal_real_type>&               omega,
                     size_t                                      trait_count) const
{
  bool conservative = my_correlation->get_parser()->get_analysis_options().conservative;

  const Matrix2D<double>& cor_type1 = ac_matrix[0][0];

  for( size_t i = 0; i < omega.rows(); ++i )
  {
    size_t type2 = (i / (trait_count*trait_count) + 1);
    const Matrix2D<double>& cor_type2 = ac_matrix[type2][0];

    size_t t1 = i % (trait_count * trait_count);

    for( size_t j = 0; j < omega.cols(); ++j )
    {
      size_t type3 = (j / (trait_count*trait_count) + 1);

      const Matrix2D<double>& cor_type3 = ac_matrix[0][type3];
      const Matrix2D<double>& cor_type4 = ac_matrix[type2][type3];

      size_t t2 = j % (trait_count * trait_count);

      if(    !SAGE::isnan(cor_type1(t1, t2))
          && !SAGE::isnan(cor_type2(t1, t2))
          && !SAGE::isnan(cor_type3(t1, t2))
          && !SAGE::isnan(cor_type4(t1, t2)) )
      {
        omega(i, j) =   (cor_type1(t1, t2) - cor_type2(t1, t2))
                      - (cor_type3(t1, t2) - cor_type4(t1, t2));
      }
      else if( conservative )
        omega(i, j) = std::numeric_limits<double>::quiet_NaN();
    }
  }

#if 0
    print_matrix(omega, cout, "omega");
#endif
}

void
Htest::compute_omega_intra(const Matrix2D<internal_real_type>& old_omega,
                                 Matrix2D<internal_real_type>& new_omega,
                           size_t                              df,
                           size_t                              trait_count) const
{
  size_t subtype_count = ((2*df) / (trait_count * (trait_count+1))) + 1;

#if 0
  cout << "compute_omega_intra..." << endl;
  cout << "old_o size = " << old_omega.rows() << "," << old_omega.cols() << endl;
  cout << ", sub count = " << subtype_count << endl;
#endif

  for( size_t sub_i = 0; sub_i < subtype_count-1; ++sub_i )
  {
    size_t new_i = sub_i * ((trait_count*(trait_count+1)) / 2);

    size_t old_i = sub_i * (trait_count*trait_count);

    for( size_t i = 0; i < trait_count; ++i )
    {
      for( size_t t1 = 0; t1 < trait_count; ++t1, ++old_i )
      {
        if( t1 < i )
          continue;

        for( size_t sub_j = 0; sub_j < subtype_count-1; ++sub_j )
        {
          size_t new_j = sub_j * ((trait_count*(trait_count+1)) / 2);

          size_t old_j = sub_j * (trait_count*trait_count);

          for( size_t j = 0; j < trait_count; ++j )
          {
            for( size_t t2 = 0; t2 < trait_count; ++t2, ++old_j )
            {
              if( t2 < j )
                continue;

              //cout << "new_omega(" << new_i << "," << new_j << ") = old_omega("
              //     << old_i << "," << old_j << ")" << endl;

              new_omega(new_i, new_j) = old_omega(old_i, old_j);

              ++new_j;
            }
          }
        }
        ++new_i;
      }
    }
  }

#if 0
    print_matrix(old_omega, cout, "old_omega");
    print_matrix(new_omega, cout, "new_omega");
#endif
}

/* Another way to calculate omega by constructing ci matrix.
   - debug purpose
void
Htest::compute_omega(const vector< vector< Matrix2D<double> > >& ac_matrix,
                     Matrix2D<internal_real_type>&               omega,
                     size_type                                   df) const
{
  size_t subtype_count = (df / (trait_count() * trait_count())) + 1;
  size_t ac_size       = subtype_count * trait_count() * trait_count();

  Matrix2D< double > ac(ac_size, ac_size, 0.0);

  for( size_t i = 0; i < ac_matrix.size(); ++i )
    for( size_t j = 0; j < ac_matrix.size(); ++j )
      for( size_t t1 = 0; t1 < ac_matrix[i][j].rows(); ++t1 )
        for( size_t t2 = 0; t2 < ac_matrix[i][j].cols(); ++t2 )
          ac(i * trait_count() * trait_count() + t1,
             j * trait_count() * trait_count() + t2) = ac_matrix[i][j](t1, t2);

  for( size_t i = 0; i < ac.rows(); ++i )
    for( size_t j = i; j < ac.cols() ; ++j )
      ac(i,j) = ac(j,i) = (ac(i,j) + ac(j,i)) / 2.;

  Matrix2D< double > ci(df, ac_size, 0.0);

  for( size_t i = 0; i < ci.rows(); ++i )
  {
    size_t t1 = i % (trait_count() * trait_count());

    for( size_t j = 0; j < trait_count() * trait_count(); ++j )
    {
      if( j == t1 )
        ci(i, j) = 1.0;
    }
  }

  for( size_t i = 0; i < ci.rows(); ++i )
  {
    size_t t1 = i % (trait_count() * trait_count());

    for( size_t j = i + trait_count() * trait_count(); j < ci.cols(); ++j )
    {
      if( j == i + trait_count() * trait_count() )
        ci(i, j) = -1.0;
    }
  }

  print_matrix(ci, cout, "ci");

  Matrix2D< double >  ci_transpose;
  Matrix2D< double >& ci_re = Transpose(ci, ci_transpose);

  Matrix2D< double >  ciac;
  Matrix2D< double >& ciac_re1 = Multiply(ci, ac, ciac);

  Matrix2D< double >  ciacciT;
  Matrix2D< double >& ciac_re2 = Multiply(ciac, ci_transpose, ciacciT);

  omega.resize(df, df, std::numeric_limits<double>::quiet_NaN());
  assert( ciacciT.rows() == omega.rows() && ciacciT.cols() == omega.cols() );

  for( size_t i = 0; i < omega.rows(); ++i )
    for( size_t j = 0; j < omega.cols(); ++j )
      omega(i, j) = ciacciT(i, j);
}
*/
double
Htest::compute_chi_square(const Matrix2D<internal_real_type>& y,
                          const Matrix2D<internal_real_type>& omega) const
{
  Matrix2D<internal_real_type> y_transpose;
  Transpose(y, y_transpose);

  Matrix2D<internal_real_type> omega_inverse;

  /* For later use when svd.h is ready.
  Matrix2D<internal_real_type> u;
  Matrix2D<internal_real_type> v;
  Matrix2D<internal_real_type> w;

  singular_value_decomposition(omega, u, w, v);

  int singular_count = singularities(w, 1.0e-12);

  omega_inverse = inverse(u, w, v);
  */

  Inverse(omega, omega_inverse);

#if 0
  print_matrix(y_transpose, cout, "y_transpose");
  print_matrix(omega, cout, "omega");
  print_matrix(omega_inverse, cout, "omega_inverse");
  print_matrix(y, cout, "y");
#endif

  Matrix2D<internal_real_type> chi_val1;
  Multiply(y_transpose, omega_inverse, chi_val1);

#if 0
  print_matrix(chi_val1, cout, "chi_val1");
#endif

  if( !chi_val1.rows() || !chi_val1.cols() )
    return std::numeric_limits<double>::quiet_NaN();

  Matrix2D<internal_real_type> chi_val2(1, 1);
  Multiply(chi_val1, y, chi_val2);

#if 0
  print_matrix(chi_val2, cout, "chi_val2");
#endif

  if( !chi_val2.rows() || !chi_val2.cols() )
    return std::numeric_limits<double>::quiet_NaN();

  assert(chi_val2.rows() == 1 && chi_val2.cols() == 1);

  /* For later use when svd.h is ready
  if( singular_count > 0 )
    return -1.0;
  */

  return chi_val2(0, 0);
}

// end of Htest Out-of-Line Implementation

} // end of namespace FCOR
} // end of namespace SAGE
