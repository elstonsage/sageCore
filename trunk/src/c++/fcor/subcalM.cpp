//****************************************************************************
//* File:      subcalM.cpp                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file calculates the asymptotic (co)variance of    *
//*            a multipedigree correlation.                                  *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/subcalM.h"

namespace SAGE {
namespace FCOR {
    
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubCalBaseM                                                  ~
// ~                                                                         ~
// ~ Purpose:   Implements the base class for SubAVCalM & SubACCalM.         ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of SubCalBaseM
// ---------------------------------------------------------------------------

SubCalBaseM::SubCalBaseM()
{
  my_correlation = NULL;

  my_pedigree_count = 0;
  my_trait_count    = 0;
  my_weight_count   = 0;
}

SubCalBaseM::SubCalBaseM(const CorrelationCal* cc)
{
  my_correlation = cc;

  my_pedigree_count = cc->get_parser()->get_pedigree_count();
  my_trait_count    = cc->get_parser()->get_trait_count();
}

SubCalBaseM::~SubCalBaseM()
{}

void
SubCalBaseM::compute_te_we(Matrix2D< vector<internal_real_type> >& wete,
                           const vector<CorrelationInfo>*          type)
{
  vector<internal_real_type> wt(my_weight_count, 0.0);

  wete.resize(trait_count(), trait_count(), wt);

  // Compute We(or Te) = 1/(W(or T) - sum(w(or t)^2)/W(or T))
  //
  for( size_t t1 = 0; t1 < wete.rows(); ++t1 )
    for( size_t t2 = 0; t2 < wete.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        double cor              = (*type)[w].correlation      (t1 + trait_count(), t2);
        internal_real_type wt   = (*type)[w].sum_weight       (t1 + trait_count(), t2);
        internal_real_type wwtt = (*type)[w].sum_weight_square(t1 + trait_count(), t2);

        if( !SAGE::isnan(cor) && wt > std::numeric_limits<double>::epsilon() )
          wete(t1, t2)[w] = 1.0 / (wt - wwtt / wt);
      }
}

void
SubCalBaseM::build_weight_matrix(Matrix2D< vector<double> >&     weight,  const size_t& i,
                                 const pairset_by_pedigree_type* pairset, const size_t& ped,
                                 const weight_matrix&            opt_weight)
{
  const pedigree_member_type* member1;
  const pedigree_member_type* member2;

  member1 = (*pairset)[ped][i].member_pair.second;
  member2 = (*pairset)[ped][i].member_pair.first;

  if( opt_weight.first.rows() && opt_weight.second.rows() )
  {
    vector<double> init_weight(my_weight_count, 0.0);
    weight.resize(trait_count(), trait_count(), init_weight);

    for( size_t t1 = 0; t1 < weight.rows(); ++t1 )
      for( size_t t2 = 0; t2 < weight.cols(); ++t2 )
      {
        for( size_t w = 0; w < my_weight_count; ++w )
          weight(t1, t2)[w] = opt_weight.first(t1, t2);
      }

#if 0
  if( ped < 3 && i < 2 )
  {
    for( size_t t1 = 0; t1 < opt_weight.first.size(); ++t1 )
    {
      for( size_t t2 = 0; t2 <= t1; ++t2 )
        cout << opt_weight.first(t1, t2) << "  ";
      cout << endl;
    }
    cout << endl;

  //  cout << "weight matrix for pedigree " << ped << ", pair " << i << endl;
  //  print_vec_matrix(weight, cout);
  }
#endif
  }
  else
    weight.resize(trait_count(), trait_count(), (*pairset)[ped][i].pair_weight);

  for( size_t t1 = 0; t1 < trait_count(); ++t1 )
    for( size_t t2 = 0; t2 < trait_count(); ++t2 )
    {
      if(    SAGE::isnan( member1->pedigree()->info().trait(member1->index(), t1) ) 
          || SAGE::isnan( member2->pedigree()->info().trait(member2->index(), t2) ) )
      {
        for( size_t w = 0; w < my_weight_count; ++w )
          weight(t1, t2)[w] = 0.0;
      }
    }

#if 0
  //if( ped < 3 && i < 2 )
  if( ped == 0 && i == 0 )
  {
    cout << "weight matrix for pedigree " << ped << ", pair " << i << endl;
    print_vec_matrix(weight, cout);
  }
#endif
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubAVCalM                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic variance of a correlation(XY & XY).     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of SubAVCalM
// ---------------------------------------------------------------------------

SubAVCalM::SubAVCalM(const CorrelationCal* cc)
          :SubCalBaseM(cc)
{ }

SubAVCalM::~SubAVCalM()
{ }

void
SubAVCalM::compute_asymptotic_variance(const pairset_by_pedigree_type&  pairset,
                                       const vector<CorrelationInfo>&   xy,
                                       const weight_matrix_by_pedigree& optimal_w,
                                       Matrix2D< vector<double> >&      av,
                                       Matrix2D< vector<bool> >&        replaced,
                                       Matrix2D< vector<size_t> >&      lp_matrix)
{
  if( !pairset.size() )
    return;

  my_corinfo_xy = &xy;
  my_pairset_xy = &pairset;

  my_conservative = true;

  my_weight_count = WEIGHT_COUNT;

  if( optimal_w.size() )
    my_weight_count = 1;

  vector<size_t> lp(my_weight_count, (size_t)-1);

  my_least_pair_count.resize(lp_matrix.rows(), lp_matrix.cols(), lp);

  init(optimal_w);
  compute(av);

  bool parser_conservative = my_correlation->get_parser()->get_analysis_options().conservative;

  if( parser_conservative )
  {
    lp_matrix = my_least_pair_count;

    return;
  }

  my_conservative = false;

  my_least_pair_count.resize(lp_matrix.rows(), lp_matrix.cols(), lp);

  vector<double> av_vec(my_weight_count, std::numeric_limits<double>::quiet_NaN());

  Matrix2D< vector<double> > av_replaced(av.rows(), av.cols(), av_vec);

  init(optimal_w);
  compute(av_replaced);

  // Compare between robust & conservative, then set replaced true if they are not same.
  //
  for( size_t i = 0; i < av.rows(); ++i )
    for( size_t j = 0; j < av.cols(); ++j )
    {
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if( av(i,j)[w] != av_replaced(i,j)[w] )
        {
          av(i,j)[w] = av_replaced(i,j)[w];
          replaced(i,j)[w] = true;
        }

        lp_matrix(i,j)[w] = my_least_pair_count(i,j)[w];
      }
    }

  return;
}

void
SubAVCalM::init(const weight_matrix_by_pedigree& optimal_w)
{
  compute_te_we(my_we, my_corinfo_xy);

  my_weight_w.resize(0);
  my_weight_w.resize(pedigree_count());

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    my_weight_w[p].resize(0);
    my_weight_w[p].resize(pair_count(my_pairset_xy, p));

    weight_matrix ped_weight;

    for( size_t i = 0; i < pair_count(my_pairset_xy, p); ++i )
    {
      if( optimal_w.size() )
        build_weight_matrix(my_weight_w[p][i], i, my_pairset_xy, p, optimal_w[p]);
      else
        build_weight_matrix(my_weight_w[p][i], i, my_pairset_xy, p, ped_weight);
    }
  }

  my_sum_kl_xx.resize(0);
  my_sum_kl_xy.resize(0);
  my_sum_kl_yy.resize(0);

  my_sum_ikj_x.resize(0);
  my_sum_ik_xy.resize(0);
  my_sum_kj_xy.resize(0);
  my_sum_ikj_y.resize(0);

  my_sum_kl_xx.resize(pedigree_count());
  my_sum_kl_xy.resize(pedigree_count());
  my_sum_kl_yy.resize(pedigree_count());

  my_sum_ikj_x.resize(pedigree_count());
  my_sum_ik_xy.resize(pedigree_count());
  my_sum_kj_xy.resize(pedigree_count());
  my_sum_ikj_y.resize(pedigree_count());

  vector<internal_real_type> init_sum_p(my_weight_count, 0.0);

  Matrix2D< vector<internal_real_type> > sum_p;

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    sum_p.resize(0, 0);
    sum_p.resize(trait_count(), trait_count(), init_sum_p);

    pedigree_weight(my_weight_w[p], sum_p);

    my_sum_ikj_x[p].resize(pair_count(my_pairset_xy, p));
    my_sum_ik_xy[p].resize(pair_count(my_pairset_xy, p));
    my_sum_kj_xy[p].resize(pair_count(my_pairset_xy, p));
    my_sum_ikj_y[p].resize(pair_count(my_pairset_xy, p));

    compute_sub(sum_p, my_sum_ikj_x[p], my_sum_kl_xx[p], p, XX);
    compute_sub(sum_p, my_sum_ik_xy[p], my_sum_kl_xy[p], p, XY);
    compute_sub(sum_p, my_sum_kj_xy[p], my_sum_kl_xy[p], p, YX);
    compute_sub(sum_p, my_sum_ikj_y[p], my_sum_kl_yy[p], p, YY);
  }
}

void
SubAVCalM::compute_sub(const Matrix2D< vector<internal_real_type> >&   sum_p,
                     vector< Matrix2D< vector<internal_real_type> > >& ikj,
                             Matrix2D< vector<internal_real_type> >&   kl, 
                       size_t p, pair_type p_t)
{
  vector<internal_real_type> init(my_weight_count, 0.0);

  kl.resize(trait_count(), trait_count(), init);

  Matrix2D< vector<internal_real_type> > local_sum_w (trait_count(), trait_count(), init);
  Matrix2D< vector<internal_real_type> > local_sum_ww(trait_count(), trait_count(), init);
    
  for( size_t i = 0; i < pair_count(my_pairset_xy, p); ++i )
  {
    ikj[i].resize(trait_count(), trait_count(), init);

    local_sum_w.resize(0, 0);
    local_sum_w.resize(trait_count(), trait_count(), init);

    for( size_t j = 0; j < pair_count(my_pairset_xy, p); ++j )
    {
      for( size_t t1 = 0; t1 < ikj[i].rows(); ++t1 )
        for( size_t t2 = 0; t2 < ikj[i].cols(); ++t2 )
        {
          vector< pair<double,size_t> > cor = get_corinfo(my_pairset_xy, my_pairset_xy,
                                                          i, j, p, p_t,
                                                          t1 + trait_count(), t2);

          for( size_t w = 0; w < my_weight_count; ++w )
          {
            if( !SAGE::isnan(cor[w].first) && !SAGE::isnan(kl(t1, t2)[w]) && !SAGE::isnan(ikj[i](t1, t2)[w]) )
            {
              // Compute:
              //  ti  = my_weight_w[p][i]
              //  wj  = my_weight_w[p][j]
              //  kl += wi*wj*kl
              kl(t1, t2)[w] += my_weight_w[p][i](t1, t2)[w] * my_weight_w[p][j](t1, t2)[w] * cor[w].first;

              // Compute:
              //  ikj += my_weight_w[p][j] * ikj
              ikj[i](t1, t2)[w] += my_weight_w[p][j](t1, t2)[w] * cor[w].first;

              local_sum_w(t1, t2)[w]  += my_weight_w[p][j](t1, t2)[w];
              local_sum_ww(t1, t2)[w] += my_weight_w[p][i](t1, t2)[w] * my_weight_w[p][j](t1, t2)[w];

              my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor[w].second);
            }
            else if( my_conservative )
            {
              kl(t1, t2)[w]     = std::numeric_limits<double>::quiet_NaN();
              ikj[i](t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
            }
          }
        }
    }

    for( size_t t1 = 0; t1 < ikj[i].rows(); ++t1 )
      for( size_t t2 = 0; t2 < ikj[i].cols(); ++t2 )
        for( size_t w = 0; w < my_weight_count; ++w )
        {
          local_sum_w(t1, t2)[w] =   (*my_corinfo_xy)[w].sum_weight(t1 + trait_count(), t2)
                                   - sum_p(t1, t2)[w] + local_sum_w(t1, t2)[w];

          if( !SAGE::isnan(ikj[i](t1, t2)[w]) && local_sum_w(t1, t2)[w] > std::numeric_limits<double>::epsilon() )
            ikj[i](t1, t2)[w] /= local_sum_w(t1, t2)[w];
        }
  }

  for( size_t t1 = 0; t1 < kl.rows(); ++t1 )
    for( size_t t2 = 0; t2 < kl.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        internal_real_type sum_ww =   (*my_corinfo_xy)[w].sum_weight(t1 + trait_count(), t2)
                                    * (*my_corinfo_xy)[w].sum_weight(t1 + trait_count(), t2);

        internal_real_type sum_pp = sum_p(t1, t2)[w] * sum_p(t1, t2)[w];

        local_sum_ww(t1, t2)[w] = sum_ww - sum_pp + local_sum_ww(t1, t2)[w];

        if( !SAGE::isnan(kl(t1, t2)[w]) && local_sum_ww(t1, t2)[w] > std::numeric_limits<double>::epsilon() )
          kl(t1, t2)[w] /= local_sum_ww(t1, t2)[w];
      }
} 

void
SubAVCalM::compute(Matrix2D< vector<double> >& ac)
{
  vector<internal_real_type> init_result(my_weight_count, 0.0);

  Matrix2D< vector<internal_real_type> > m1 (trait_count(), trait_count(), init_result);
  Matrix2D< vector<internal_real_type> > m23(trait_count(), trait_count(), init_result);
  Matrix2D< vector<internal_real_type> > m4 (trait_count(), trait_count(), init_result);
  Matrix2D< vector<internal_real_type> > m57(trait_count(), trait_count(), init_result);
  Matrix2D< vector<internal_real_type> > m68(trait_count(), trait_count(), init_result);
  Matrix2D< vector<internal_real_type> > m9 (trait_count(), trait_count(), init_result);

  Matrix2D< vector<internal_real_type> > sub_xx;
  Matrix2D< vector<internal_real_type> > sub_xy;
  Matrix2D< vector<internal_real_type> > sub_yx;
  Matrix2D< vector<internal_real_type> > sub_yy;

  for( size_t p = 0; p < pedigree_count(); ++p )
    for( size_t i = 0; i < pair_count(my_pairset_xy, p); ++i )
      for( size_t j = 0; j < pair_count(my_pairset_xy, p); ++j )
      {
        sub_module(sub_xx, i, j, p, XX);
        sub_module(sub_xy, i, j, p, XY);
        sub_module(sub_yx, i, j, p, YX);
        sub_module(sub_yy, i, j, p, YY);

        module1_8(sub_xx, sub_xx, m1,  i, j, p);
        module1_8(sub_xy, sub_xy, m23, i, j, p);
        module1_8(sub_yy, sub_yy, m4,  i, j, p);

        module1_8(sub_xx, sub_xy, m57, i, j, p);
        module1_8(sub_yx, sub_yy, m68, i, j, p);

        module9(sub_xx, sub_yy, sub_xy, sub_yx, m9, i, j, p);
      }

  vector<double> init(my_weight_count, 0.0);

  Matrix2D< vector<double> > cor_xy  (trait_count(), trait_count(), init);
  Matrix2D< vector<double> > cov_xy  (trait_count(), trait_count(), init);
  Matrix2D< vector<double> > cov_xxyy(trait_count(), trait_count(), init);

  for( size_t t1 = 0; t1 < cov_xy.rows(); ++t1 )
    for( size_t t2 = 0; t2 < cov_xy.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        cor_xy(t1, t2)[w] = (*my_corinfo_xy)[w].correlation(t1 + trait_count(), t2);
        cov_xy(t1, t2)[w] = (*my_corinfo_xy)[w].covariance (t1 + trait_count(), t2);

        double xx = (*my_corinfo_xy)[w].covariance(t1 + trait_count(), t1 + trait_count());
        double yy = (*my_corinfo_xy)[w].covariance(t2, t2);

        if( !SAGE::isnan(xx) && !SAGE::isnan(yy) )
          cov_xxyy(t1, t2)[w] = sqrt(xx * yy);
        else if( my_conservative )
          cov_xxyy(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }

  for( size_t t1 = 0; t1 < m1.rows(); ++t1 )
    for( size_t t2 = 0; t2 < m1.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(m1 (t1, t2)[w]) && !SAGE::isnan(m23(t1, t2)[w]) && !SAGE::isnan(m4(t1, t2)[w])
            && !SAGE::isnan(m57(t1, t2)[w]) && !SAGE::isnan(m68(t1, t2)[w]) && !SAGE::isnan(m9(t1, t2)[w])
            && !SAGE::isnan(cov_xxyy(t1, t2)[w]) && !SAGE::isnan(cor_xy(t1, t2)[w])
            && fabs(cov_xy(t1, t2)[w]) > std::numeric_limits<double>::epsilon()  )
        {
          m1 (t1, t2)[w] *= 2.0 * my_we(t1, t2)[w] * my_we(t1, t2)[w];
          m23(t1, t2)[w] *= 2.0 * my_we(t1, t2)[w] * my_we(t1, t2)[w];
          m4 (t1, t2)[w] *= 2.0 * my_we(t1, t2)[w] * my_we(t1, t2)[w];
          m57(t1, t2)[w] *= 2.0 * my_we(t1, t2)[w] * my_we(t1, t2)[w] * cov_xxyy(t1, t2)[w];
          m68(t1, t2)[w] *= 2.0 * my_we(t1, t2)[w] * my_we(t1, t2)[w] * cov_xxyy(t1, t2)[w];
          m9 (t1, t2)[w] *=       my_we(t1, t2)[w] * my_we(t1, t2)[w] * cov_xxyy(t1, t2)[w] * cov_xxyy(t1, t2)[w];

          // Compute: m57, m68, m9 /= cov
          m57(t1, t2)[w] /= cov_xy(t1, t2)[w];
          m68(t1, t2)[w] /= cov_xy(t1, t2)[w];
          m9 (t1, t2)[w] /= ( cov_xy(t1, t2)[w] * cov_xy(t1, t2)[w] );

          double m0  = cor_xy(t1, t2)[w] * cor_xy(t1, t2)[w] / 4.0;

          // Compute ac = m0*(m1+m23+m23+m4 - 2*(m57+m68+m57+m68) + 4*m9)
          ac(t1, t2)[w] =   m0  * ( m1 (t1, t2)[w] + m23(t1, t2)[w] + m23(t1, t2)[w] + m4 (t1, t2)[w]
                          - 2.0 * ( m57(t1, t2)[w] + m68(t1, t2)[w] + m57(t1, t2)[w] + m68(t1, t2)[w] )
                          + 4.0 *   m9 (t1, t2)[w] );
        }
      }
}

void
SubAVCalM::module1_8(const Matrix2D< vector<internal_real_type> >& s1,
                     const Matrix2D< vector<internal_real_type> >& s2,
                           Matrix2D< vector<internal_real_type> >& result,
                     size_t i, size_t j, size_t p) const
{
  // Compute: result += s1 * s2 * my_weight_w[p][i] * my_weight_w[p][j]
  for( size_t t1 = 0; t1 < result.rows(); ++t1 )
    for( size_t t2 = 0; t2 < result.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(result(t1, t2)[w])
            && !SAGE::isnan(s1(t1, t2)[w]) && !SAGE::isnan(s2(t1, t2)[w]) )
        {
          result(t1, t2)[w] +=   s1(t1, t2)[w] * my_weight_w[p][i](t1, t2)[w]
                               * s2(t1, t2)[w] * my_weight_w[p][j](t1, t2)[w];
        }
        else if( my_conservative )
          result(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }
}

void
SubAVCalM::module9(const Matrix2D< vector<internal_real_type> >& s_xx,
                   const Matrix2D< vector<internal_real_type> >& s_yy, 
                   const Matrix2D< vector<internal_real_type> >& s_xy,
                   const Matrix2D< vector<internal_real_type> >& s_yx, 
                         Matrix2D< vector<internal_real_type> >& result,
                   size_t i, size_t j, size_t p) const
{
  // Compute: result += (s_xx*s_yy + s_xy*s_yx) * wi * wj
  for( size_t t1 = 0; t1 < result.rows(); ++t1 )
    for( size_t t2 = 0; t2 < result.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(result(t1, t2)[w])
            && !SAGE::isnan(s_xx(t1, t2)[w]) && !SAGE::isnan(s_yy(t1, t2)[w])
            && !SAGE::isnan(s_xy(t1, t2)[w]) && !SAGE::isnan(s_yx(t1, t2)[w]) )
        {
          result(t1, t2)[w] +=   ( s_xx(t1, t2)[w] * s_yy(t1, t2)[w] + s_xy(t1, t2)[w] * s_yx(t1, t2)[w] )
                               *   my_weight_w[p][i](t1, t2)[w]
                               *   my_weight_w[p][j](t1, t2)[w];
        }
        else if( my_conservative )
          result(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }
}

void
SubAVCalM::sub_module(Matrix2D< vector<internal_real_type> >& sub,
                      size_t i, size_t j, size_t p, pair_type p_t)
{
  vector<internal_real_type> init_result(my_weight_count, 0.0);

  sub.resize(0, 0);
  sub.resize(trait_count(), trait_count(), init_result);

  vector<double> init(my_weight_count, 0.0);

  Matrix2D< vector<double> > cor(trait_count(), trait_count(), init);

  for( size_t t1 = 0; t1 < cor.rows(); ++t1 )
  {
    for( size_t t2 = 0; t2 < cor.cols(); ++t2 )
    {
      vector< pair<double, size_t> > cor_pair = get_corinfo(my_pairset_xy, my_pairset_xy,
                                                            i, j, p, p_t,
                                                            t1 + trait_count(), t2);

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if( !SAGE::isnan(cor_pair[w].first) )
        {
          cor(t1, t2)[w] = cor_pair[w].first;

          my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor_pair[w].second);
        }
        else if(my_conservative)
          cor(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  const Matrix2D< vector<internal_real_type> >* ik = NULL;
  const Matrix2D< vector<internal_real_type> >* kj = NULL;
  const Matrix2D< vector<internal_real_type> >* kl = NULL;
  
  if( p_t == XX )
  {
    ik = &my_sum_ikj_x[p][i];
    kj = &my_sum_ikj_x[p][j];
    kl = &my_sum_kl_xx[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_xx_type(sub, cor, ik, kj, kl);

    return;
  }

  if( p_t == XY )
  {
    ik = &my_sum_ik_xy[p][i];
    kj = &my_sum_kj_xy[p][j];
    kl = &my_sum_kl_xy[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_xy_type(sub, cor, ik, kj, kl);

    return;
  }

  if( p_t == YX )
  {
    ik = &my_sum_kj_xy[p][i];
    kj = &my_sum_ik_xy[p][j];
    kl = &my_sum_kl_xy[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_yx_type(sub, cor, ik, kj, kl);

    return;
  } 

  if( p_t == YY )
  {
    ik = &my_sum_ikj_y[p][i];
    kj = &my_sum_ikj_y[p][j];
    kl = &my_sum_kl_yy[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_yy_type(sub, cor, ik, kj, kl);

    return;
  }

  assert( ik != NULL && kj != NULL && kl != NULL );

}

void
SubAVCalM::sub_xx_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1 = 0; t1 < sub.rows(); ++t1 )
    for( size_t t2 = 0; t2 < sub.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
        if(    !SAGE::isnan(cor(t1, t1)[w])   && !SAGE::isnan((*ik)(t1, t1)[w])
            && !SAGE::isnan((*kj)(t1, t1)[w]) && !SAGE::isnan((*kl)(t1, t1)[w]) )
        {
          sub(t1, t2)[w] = cor(t1, t1)[w] - (*ik)(t1, t1)[w] - (*kj)(t1, t1)[w] + (*kl)(t1, t1)[w];

          my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w],
                                                    my_least_pair_count(t1, t1)[w]);
        }
        else if( my_conservative )
          sub(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
}

void
SubAVCalM::sub_xy_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1 = 0; t1 < sub.rows(); ++t1 )
    for( size_t t2 = 0; t2 < sub.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
        if(    !SAGE::isnan(cor(t1, t2)[w])   && !SAGE::isnan((*ik)(t1, t2)[w])
            && !SAGE::isnan((*kj)(t1, t2)[w]) && !SAGE::isnan((*kl)(t1, t2)[w]) )
        {
          sub(t1, t2)[w] = cor(t1, t2)[w] - (*ik)(t1, t2)[w] - (*kj)(t1, t2)[w] + (*kl)(t1, t2)[w];
        }
        else if( my_conservative )
          sub(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
}

void
SubAVCalM::sub_yx_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1 = 0; t1 < sub.rows(); ++t1 )
    for( size_t t2 = 0; t2 < sub.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
        if(    !SAGE::isnan(cor(t2, t1)[w])   && !SAGE::isnan((*ik)(t2, t1)[w])
            && !SAGE::isnan((*kj)(t2, t1)[w]) && !SAGE::isnan((*kl)(t2, t1)[w]) )
        {
          sub(t1, t2)[w] = cor(t2, t1)[w] - (*ik)(t2, t1)[w] - (*kj)(t2, t1)[w] + (*kl)(t2, t1)[w];

          my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w],
                                                    my_least_pair_count(t2, t1)[w]);
        }
        else if( my_conservative )
          sub(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
}

void
SubAVCalM::sub_yy_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1 = 0; t1 < sub.rows(); ++t1 )
    for( size_t t2 = 0; t2 < sub.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
        if(    !SAGE::isnan(cor(t2, t2)[w])   && !SAGE::isnan((*ik)(t2, t2)[w])
            && !SAGE::isnan((*kj)(t2, t2)[w]) && !SAGE::isnan((*kl)(t2, t2)[w]) )
        {
          sub(t1, t2)[w] = cor(t2, t2)[w] - (*ik)(t2, t2)[w] - (*kj)(t2, t2)[w] + (*kl)(t2, t2)[w];

          my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w],
                                                    my_least_pair_count(t2, t2)[w]);
        }
        else if( my_conservative )
          sub(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
}

// end of SubAVCalM Implementation

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubACCalM                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic covariance of two correlations(XY & UV).~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of SubACCalM
// ---------------------------------------------------------------------------

SubACCalM::SubACCalM(const CorrelationCal* cc)
          :SubCalBaseM(cc)
{ }

SubACCalM::~SubACCalM()
{ }

void
SubACCalM::compute_asymptotic_covariance(const pairset_by_pedigree_type&  pairset1,
                                         const pairset_by_pedigree_type&  pairset2,
                                         const vector<CorrelationInfo>&   xy,
                                         const vector<CorrelationInfo>&   uv,
                                         const weight_matrix_by_pedigree& optimal_w1,
                                         const weight_matrix_by_pedigree& optimal_w2,
                                         Matrix2D< vector<double> >&      ac,
                                         Matrix2D< vector<bool> >&        replaced,
                                         Matrix2D< vector<size_t> >&      lp_matrix)
{
  if( !pairset1.size() || !pairset2.size() )
    return;

  my_corinfo_xy = &xy;
  my_corinfo_uv = &uv;
  
  my_pairset_xy = &pairset1;
  my_pairset_uv = &pairset2;

  my_conservative = true;

  my_weight_count = WEIGHT_COUNT;

  if( optimal_w1.size() )
    my_weight_count = 1;

  vector<size_t> lp(my_weight_count, (size_t)-1);

  my_least_pair_count.resize(trait_count(), trait_count(), lp);
  my_extended_least_pair_count.resize(lp_matrix.rows(), lp_matrix.cols(), lp);

  init(optimal_w1, optimal_w2);
  compute(ac);

  bool parser_conservative = my_correlation->get_parser()->get_analysis_options().conservative;

  if( parser_conservative )
  {
    lp_matrix = my_extended_least_pair_count;

    return;
  }

  my_conservative = false;

  my_least_pair_count.resize(0, 0);
  my_extended_least_pair_count.resize(0, 0);

  my_least_pair_count.resize(trait_count(), trait_count(), lp);
  my_extended_least_pair_count.resize(lp_matrix.rows(), lp_matrix.cols(), lp);

  vector<double> ac_vec(my_weight_count, std::numeric_limits<double>::quiet_NaN());

  Matrix2D< vector<double> > ac_replaced(ac.rows(), ac.cols(), ac_vec);

  init(optimal_w1, optimal_w2);
  compute(ac_replaced);

#if 0
  print_vec_matrix(my_least_pair_count, cout, "least_pair_count");
  print_vec_matrix(my_extended_least_pair_count, cout, "extended_least_pair_count");
#endif

  // Compare between robust & conservative, then set replaced true if they are not same.
  //
  for( size_t i = 0; i < ac.rows(); ++i )
    for( size_t j = 0; j < ac.cols(); ++j )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if( ac(i,j)[w] != ac_replaced(i,j)[w] )
        {
          ac(i,j)[w] = ac_replaced(i,j)[w];
          replaced(i,j)[w] = true;
        }

        lp_matrix(i,j)[w] = my_extended_least_pair_count(i,j)[w];
      }

  return;
}

void
SubACCalM::init(const weight_matrix_by_pedigree& w1, const weight_matrix_by_pedigree& w2)
{
  compute_te_we(my_we, my_corinfo_xy);
  compute_te_we(my_te, my_corinfo_uv);

  my_weight_w.resize(0);
  my_weight_t.resize(0);

  my_weight_w.resize(pedigree_count());
  my_weight_t.resize(pedigree_count());

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    my_weight_w[p].resize(0);
    my_weight_t[p].resize(0);

    my_weight_w[p].resize(pair_count(my_pairset_xy, p));
    my_weight_t[p].resize(pair_count(my_pairset_uv, p));

    weight_matrix ped_weight;

    for( size_t i = 0; i < pair_count(my_pairset_xy, p); ++i )
    {
      if( w1.size() )
        build_weight_matrix(my_weight_w[p][i], i, my_pairset_xy, p, w1[p]);
      else
        build_weight_matrix(my_weight_w[p][i], i, my_pairset_xy, p, ped_weight);
    }

    for( size_t i = 0; i < pair_count(my_pairset_uv, p); ++i )
    {
      if( w2.size() )
        build_weight_matrix(my_weight_t[p][i], i, my_pairset_uv, p, w2[p]);
      else
        build_weight_matrix(my_weight_t[p][i], i, my_pairset_uv, p, ped_weight);
    }
  }

  my_sum_kl_ux.resize(0);
  my_sum_kl_uy.resize(0);
  my_sum_kl_vx.resize(0);
  my_sum_kl_vy.resize(0);

  my_sum_ik_ux.resize(0);
  my_sum_ik_uy.resize(0);
  my_sum_ik_vx.resize(0);
  my_sum_ik_vy.resize(0);

  my_sum_kj_ux.resize(0);
  my_sum_kj_uy.resize(0);
  my_sum_kj_vx.resize(0);
  my_sum_kj_vy.resize(0);

  my_sum_kl_ux.resize(pedigree_count());
  my_sum_kl_uy.resize(pedigree_count());
  my_sum_kl_vx.resize(pedigree_count());
  my_sum_kl_vy.resize(pedigree_count());

  my_sum_ik_ux.resize(pedigree_count());
  my_sum_ik_uy.resize(pedigree_count());
  my_sum_ik_vx.resize(pedigree_count());
  my_sum_ik_vy.resize(pedigree_count());

  my_sum_kj_ux.resize(pedigree_count());
  my_sum_kj_uy.resize(pedigree_count());
  my_sum_kj_vx.resize(pedigree_count());
  my_sum_kj_vy.resize(pedigree_count());

  vector<internal_real_type> init_sum_p(my_weight_count, 0.0);

  Matrix2D< vector<internal_real_type> > sum_pw;
  Matrix2D< vector<internal_real_type> > sum_pt;

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    sum_pw.resize(0, 0);
    sum_pt.resize(0, 0);

    sum_pw.resize(trait_count(), trait_count(), init_sum_p);
    sum_pt.resize(trait_count(), trait_count(), init_sum_p);

    pedigree_weight(my_weight_w[p], sum_pw);
    pedigree_weight(my_weight_t[p], sum_pt);

    my_sum_ik_ux[p].resize(pair_count(my_pairset_uv, p));
    my_sum_ik_uy[p].resize(pair_count(my_pairset_uv, p));
    my_sum_ik_vx[p].resize(pair_count(my_pairset_uv, p));
    my_sum_ik_vy[p].resize(pair_count(my_pairset_uv, p));

    for( size_t i = 0; i < pair_count(my_pairset_uv, p); ++i )
    {
      compute_ik(sum_pw, my_sum_ik_ux[p][i], i, p, UX);
      compute_ik(sum_pw, my_sum_ik_uy[p][i], i, p, UY);
      compute_ik(sum_pw, my_sum_ik_vx[p][i], i, p, VX);
      compute_ik(sum_pw, my_sum_ik_vy[p][i], i, p, VY);
    }

    my_sum_kj_ux[p].resize(pair_count(my_pairset_xy, p));
    my_sum_kj_uy[p].resize(pair_count(my_pairset_xy, p));
    my_sum_kj_vx[p].resize(pair_count(my_pairset_xy, p));
    my_sum_kj_vy[p].resize(pair_count(my_pairset_xy, p));

    for( size_t j = 0; j < pair_count(my_pairset_xy, p); ++j )
    {
      compute_kj(sum_pt, my_sum_kj_ux[p][j], j, p, UX);
      compute_kj(sum_pt, my_sum_kj_uy[p][j], j, p, UY);
      compute_kj(sum_pt, my_sum_kj_vx[p][j], j, p, VX);
      compute_kj(sum_pt, my_sum_kj_vy[p][j], j, p, VY);
    }

    compute_kl(sum_pt, sum_pw, my_sum_kl_ux[p], p, UX);
    compute_kl(sum_pt, sum_pw, my_sum_kl_uy[p], p, UY);
    compute_kl(sum_pt, sum_pw, my_sum_kl_vx[p], p, VX);
    compute_kl(sum_pt, sum_pw, my_sum_kl_vy[p], p, VY);
  }
}

void
SubACCalM::compute_ik(const Matrix2D< vector<internal_real_type> >& sum_pw,
                            Matrix2D< vector<internal_real_type> >& ik,
                            size_t i, size_t p, pair_type p_t)
{
  vector<internal_real_type> init(my_weight_count, 0.0);

  ik.resize(trait_count(), trait_count(), init);

  Matrix2D< vector<internal_real_type> > local_sum_w(trait_count(), trait_count(), init);
  
  for( size_t k = 0; k < pair_count(my_pairset_xy, p); ++k )
  {
    for( size_t t1 = 0; t1 < ik.rows(); ++t1 )
    {
      for( size_t t2 = 0; t2 < ik.cols(); ++t2 )
      {
        vector< pair<double, size_t> > cor = get_corinfo(my_pairset_uv, my_pairset_xy,
                                                         i, k, p, p_t,
                                                         t1 + trait_count(), t2);
        for( size_t w = 0; w < my_weight_count; ++w )
        {
          // Compute: ik += my_weight_w[p][k] * ik    
          if( !SAGE::isnan(cor[w].first) && !SAGE::isnan(ik(t1, t2)[w]) )
          {
            ik(t1, t2)[w]          += my_weight_w[p][k](t1, t2)[w] * cor[w].first;
            local_sum_w(t1, t2)[w] += my_weight_w[p][k](t1, t2)[w];

            my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor[w].second);
          }
          else if( my_conservative )
            ik(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }

#if 0
  print_vec_matrix(my_least_pair_count, cout, "least_pair_count in ik()");
#endif

  // Compute : ik /= W
  // W has been replaced with local_W since we replaced nan_cor with 0.
  //
  for( size_t t1 = 0; t1 < ik.rows(); ++t1 )
  {
    for( size_t t2 = 0; t2 < ik.cols(); ++t2 )
    {
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        local_sum_w(t1, t2)[w] =   (*my_corinfo_xy)[w].sum_weight(t1 + trait_count(), t2)
                                 - sum_pw(t1, t2)[w] + local_sum_w(t1, t2)[w]; 

        if( !SAGE::isnan(ik(t1, t2)[w]) && local_sum_w(t1, t2)[w] > std::numeric_limits<double>::epsilon() )
          ik(t1, t2)[w] /= local_sum_w(t1, t2)[w];
      }
    }
  }
}

void
SubACCalM::compute_kj(const Matrix2D< vector<internal_real_type> >& sum_pt,
                            Matrix2D< vector<internal_real_type> >& kj,
                            size_t j, size_t p, pair_type p_t)
{
  vector<internal_real_type> init(my_weight_count, 0.0);

  kj.resize(trait_count(), trait_count(), init);
  
  Matrix2D< vector<internal_real_type> > local_sum_t(trait_count(), trait_count(), init);

  for( size_t k = 0; k < pair_count(my_pairset_uv, p); ++k )
  {
    for(size_t t1 = 0; t1 < kj.rows(); ++t1)
    {
      for(size_t t2 = 0; t2 < kj.cols(); ++t2)
      {
        vector< pair<double, size_t> > cor = get_corinfo(my_pairset_uv, my_pairset_xy,
                                                         k, j, p, p_t,
                                                         t1 + trait_count(), t2);

        for( size_t w = 0; w < my_weight_count; ++w )
        {
          // Compute: kj += my_weight_t[p][k] * kj    
          if( !SAGE::isnan(cor[w].first) && !SAGE::isnan(kj(t1,t2)[w]) )
          {
            kj(t1, t2)[w]          += my_weight_t[p][k](t1, t2)[w] * cor[w].first;
            local_sum_t(t1, t2)[w] += my_weight_t[p][k](t1, t2)[w];

            my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor[w].second);
          }
          else if( my_conservative )
            kj(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }
  }

#if 0
  print_vec_matrix(my_least_pair_count, cout, "least_pair_count in kj()");
#endif

  // Compute : kj /= T
  // T has been replaced with local_T since we replaced nan_cor with 0.
  //
  for( size_t t1 = 0; t1 < kj.rows(); ++t1 )
  {
    for( size_t t2 = 0; t2 < kj.cols(); ++t2 )
    {
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        local_sum_t(t1, t2)[w] =   (*my_corinfo_uv)[w].sum_weight(t1 + trait_count(), t2)
                                 - sum_pt(t1, t2)[w] + local_sum_t(t1, t2)[w];

        if( !SAGE::isnan(kj(t1, t2)[w]) && local_sum_t(t1, t2)[w] > std::numeric_limits<double>::epsilon() )
          kj(t1, t2)[w] /= local_sum_t(t1, t2)[w];
      }
    }
  }
} 

void
SubACCalM::compute_kl(const Matrix2D< vector<internal_real_type> >& sum_pt,
                      const Matrix2D< vector<internal_real_type> >& sum_pw,
                            Matrix2D< vector<internal_real_type> >& kl,
                            size_t p, pair_type p_t)
{
  vector<internal_real_type> init(my_weight_count, 0.0);

  kl.resize(trait_count(), trait_count(), init);

  Matrix2D< vector<internal_real_type> > local_sum_tw(trait_count(), trait_count(), init);
    
  for( size_t k = 0; k < pair_count(my_pairset_uv, p); ++k )
  {
    for( size_t l = 0; l < pair_count(my_pairset_xy, p); ++l )
    {
      for(size_t t1 = 0; t1 < kl.rows(); ++t1)
      {
        for(size_t t2 = 0; t2 < kl.cols(); ++t2)
        {
          vector< pair<double, size_t> > cor = get_corinfo(my_pairset_uv, my_pairset_xy,
                                                           k, l, p, p_t,
                                                           t1 + trait_count(), t2);

          for( size_t w = 0; w < my_weight_count; ++w )
          {
            if( !SAGE::isnan(cor[w].first) && !SAGE::isnan(kl(t1, t2)[w]) )
            {
              // Compute:
              //  tk  = my_weight_t[p][k]
              //  wl  = my_weight_w[p][l]
              //  kl += tk*wl*kl
              kl(t1, t2)[w]  += my_weight_t[p][k](t1, t2)[w] * my_weight_w[p][l](t1, t2)[w] * cor[w].first;

              local_sum_tw(t1, t2)[w] += my_weight_t[p][k](t1, t2)[w] * my_weight_w[p][l](t1, t2)[w];

              my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor[w].second);
            }
            else if( my_conservative )
              kl(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          }
        }
      }
    }
  }

#if 0
  print_vec_matrix(my_least_pair_count, cout, "least_pair_count in kl()");
#endif

  // Compute : kl /= TW
  // TW has been replaced with local_TW since we replaced nan_cor with 0.
  //
  for(size_t t1 = 0; t1 < kl.rows(); ++t1)
  {
    for(size_t t2 = 0; t2 < kl.cols(); ++t2)
    {
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        internal_real_type sum_tw   =   (*my_corinfo_xy)[w].sum_weight(t1 + trait_count(), t2)
                                      * (*my_corinfo_uv)[w].sum_weight(t1 + trait_count(), t2);                                             

        internal_real_type sum_ptpw = sum_pt(t1, t2)[w] * sum_pw(t1, t2)[w];                        

        local_sum_tw(t1, t2)[w] = sum_tw - sum_ptpw + local_sum_tw(t1, t2)[w];

        if( !SAGE::isnan(kl(t1, t2)[w]) && local_sum_tw(t1, t2)[w] > std::numeric_limits<double>::epsilon() )
          kl(t1, t2)[w] /= local_sum_tw(t1, t2)[w];
      }
    }
  }
} 

void
SubACCalM::compute(Matrix2D< vector<double> >& ac)
{
  size_t matrix_size = trait_count() * trait_count();

  vector<internal_real_type> init_result(my_weight_count, 0.0);

  Matrix2D< vector<internal_real_type> > m1(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m2(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m3(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m4(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m5(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m6(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m7(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m8(matrix_size, matrix_size, init_result);
  Matrix2D< vector<internal_real_type> > m9(matrix_size, matrix_size, init_result);

  Matrix2D< vector<internal_real_type> > sub_ux;
  Matrix2D< vector<internal_real_type> > sub_uy;
  Matrix2D< vector<internal_real_type> > sub_vx;
  Matrix2D< vector<internal_real_type> > sub_vy;

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    for( size_t i = 0; i < pair_count(my_pairset_uv, p); ++i )
    {
      for( size_t j = 0; j < pair_count(my_pairset_xy, p); ++j )
      {
        sub_module(sub_ux, i, j, p, UX);
        sub_module(sub_uy, i, j, p, UY);
        sub_module(sub_vx, i, j, p, VX);
        sub_module(sub_vy, i, j, p, VY);

        module1_8(sub_ux, sub_ux, m1, i, j, p);
        module1_8(sub_uy, sub_uy, m2, i, j, p);
        module1_8(sub_vx, sub_vx, m3, i, j, p);
        module1_8(sub_vy, sub_vy, m4, i, j, p);

        module1_8(sub_ux, sub_uy, m5, i, j, p);
        module1_8(sub_vx, sub_vy, m6, i, j, p);
        module1_8(sub_ux, sub_vx, m7, i, j, p);
        module1_8(sub_uy, sub_vy, m8, i, j, p);

        module9(sub_ux, sub_vy, sub_uy, sub_vx, m9, i, j, p);
      }
    }
  }

  vector<double> init(my_weight_count, 0.0);

  Matrix2D< vector<double> > cor_xy(trait_count(), trait_count(), init);
  Matrix2D< vector<double> > cor_uv(trait_count(), trait_count(), init);

  Matrix2D< vector<double> > cov_xy(trait_count(), trait_count(), init);
  Matrix2D< vector<double> > cov_uv(trait_count(), trait_count(), init);

  Matrix2D< vector<double> > cov_xxyy(trait_count(), trait_count(), init);
  Matrix2D< vector<double> > cov_uuvv(trait_count(), trait_count(), init);

  for( size_t t1 = 0; t1 < cov_xy.rows(); ++t1 )
    for( size_t t2 = 0; t2 < cov_xy.cols(); ++t2 )
      for( size_t w = 0; w < my_weight_count; ++w )
      {
        cor_xy(t1, t2)[w] = (*my_corinfo_xy)[w].correlation(t1 + trait_count(), t2);     
        cor_uv(t1, t2)[w] = (*my_corinfo_uv)[w].correlation(t1 + trait_count(), t2);

        cov_xy(t1, t2)[w] = (*my_corinfo_xy)[w].covariance(t1 + trait_count(), t2);
        cov_uv(t1, t2)[w] = (*my_corinfo_uv)[w].covariance(t1 + trait_count(), t2);

        double v_xx = (*my_corinfo_xy)[w].covariance(t1 + trait_count(), t1 + trait_count());     
        double v_uu = (*my_corinfo_uv)[w].covariance(t1 + trait_count(), t1 + trait_count());     
        double v_yy = (*my_corinfo_xy)[w].covariance(t2, t2);
        double v_vv = (*my_corinfo_uv)[w].covariance(t2, t2);

        if(    !SAGE::isnan(v_xx) && !SAGE::isnan(v_yy)  
            && !SAGE::isnan(v_uu) && !SAGE::isnan(v_vv) )
        {
          cov_xxyy(t1, t2)[w] = sqrt(v_xx * v_yy);
          cov_uuvv(t1, t2)[w] = sqrt(v_uu * v_vv);     
        }
        else if( my_conservative )
        {
          cov_xxyy(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          cov_uuvv(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
        }
      }

  for( size_t t1 = 0, t1_t = 0, t2_t = 0; t1 < m1.rows(); ++t1, ++t2_t )
  {
    if( t1 != 0 && t1 % trait_count() == 0 )
    {
      ++t1_t;
      t2_t = 0;
    }

    for( size_t t2 = 0, t1_w = 0, t2_w = 0; t2 < m1.cols(); ++t2, ++t2_w )
    {
      if( t2 != 0 && t2 % trait_count() == 0 )
      {
        ++t1_w;
        t2_w = 0;
      }

      for( size_t w = 0; w < my_weight_count; ++w )
      {
#if 0
  cout << "constituent cor & cov for t1 = " << t1 << ", t2 = " << t2 << " :" << endl;
  cout.setf(ios_base::scientific, ios_base::floatfield);
  cout << "m1 = " << m1(t1, t2)[w] << endl
       << "m2 = " << m2(t1, t2)[w] << endl
       << "m3 = " << m3(t1, t2)[w] << endl
       << "m4 = " << m4(t1, t2)[w] << endl
       << "m5 = " << m5(t1, t2)[w] << endl
       << "m6 = " << m6(t1, t2)[w] << endl
       << "m7 = " << m7(t1, t2)[w] << endl
       << "m8 = " << m8(t1, t2)[w] << endl
       << "m9 = " << m9(t1, t2)[w] << endl
       << "cov_xxyy = " << cov_xxyy(t1_w, t2_w)[w] << endl
       << "cor_xy   = " << cor_xy(t1_w, t2_w)[w] << endl
       << "cov_uuvv = " << cov_uuvv(t1_t, t2_t)[w] << endl
       << "cor_uv   = " << cor_uv(t1_t, t2_t)[w] << endl
       << "cov_xy   = " << cov_xy(t1_w, t2_w)[w] << endl
       << "cov_uv   = " << cov_uv(t1_t, t2_t) << endl
       << endl;
#endif

        if(    !SAGE::isnan(m1(t1, t2)[w]) && !SAGE::isnan(m2(t1, t2)[w]) && !SAGE::isnan(m3(t1, t2)[w])
            && !SAGE::isnan(m4(t1, t2)[w]) && !SAGE::isnan(m5(t1, t2)[w]) && !SAGE::isnan(m6(t1, t2)[w])
            && !SAGE::isnan(m7(t1, t2)[w]) && !SAGE::isnan(m8(t1, t2)[w]) && !SAGE::isnan(m9(t1, t2)[w]) )
        {
          m1(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m2(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m3(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m4(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m5(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m6(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m7(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m8(t1, t2)[w] *= 2.0 * my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
          m9(t1, t2)[w] *=       my_te(t1_t, t2_t)[w] * my_we(t1_w, t2_w)[w];
        }

        if(    !SAGE::isnan(m5(t1, t2)[w]) && !SAGE::isnan(m6(t1, t2)[w]) && !SAGE::isnan(m7(t1, t2)[w])
            && !SAGE::isnan(m8(t1, t2)[w]) && !SAGE::isnan(m9(t1, t2)[w])
            && !SAGE::isnan(cov_xxyy(t1_w, t2_w)[w])
            && !SAGE::isnan(cov_uuvv(t1_t, t2_t)[w]) )
        {
          m5(t1, t2)[w] *= cov_xxyy(t1_w, t2_w)[w];
          m6(t1, t2)[w] *= cov_xxyy(t1_w, t2_w)[w];
          m7(t1, t2)[w] *= cov_uuvv(t1_t, t2_t)[w];
          m8(t1, t2)[w] *= cov_uuvv(t1_t, t2_t)[w];
          m9(t1, t2)[w] *= ( cov_xxyy(t1_w, t2_w)[w] * cov_uuvv(t1_t, t2_t)[w] );
        }
        else if( my_conservative )
        {
          m5(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          m6(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          m7(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          m8(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
          m9(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
        }

        if(    !SAGE::isnan(m1(t1, t2)[w]) && !SAGE::isnan(m2(t1, t2)[w]) && !SAGE::isnan(m3(t1, t2)[w])
            && !SAGE::isnan(m4(t1, t2)[w]) && !SAGE::isnan(m5(t1, t2)[w]) && !SAGE::isnan(m6(t1, t2)[w])
            && !SAGE::isnan(m7(t1, t2)[w]) && !SAGE::isnan(m8(t1, t2)[w]) && !SAGE::isnan(m9(t1, t2)[w])
            && !SAGE::isnan(cor_xy(t1_w, t2_w)[w])
            && !SAGE::isnan(cor_uv(t1_t, t2_t)[w])
            && fabs(cov_xy(t1_w,t2_w)[w]) > std::numeric_limits<double>::epsilon()
            && fabs(cov_uv(t1_t,t2_t)[w]) > std::numeric_limits<double>::epsilon() )
        {
          m5(t1, t2)[w] /= cov_xy(t1_w, t2_w)[w];
          m6(t1, t2)[w] /= cov_xy(t1_w, t2_w)[w];
          m7(t1, t2)[w] /= cov_uv(t1_t, t2_t)[w];
          m8(t1, t2)[w] /= cov_uv(t1_t, t2_t)[w];
          m9(t1, t2)[w] /= ( cov_xy(t1_w, t2_w)[w] * cov_uv(t1_t, t2_t)[w] );

          double m0 = cor_xy(t1_w, t2_w)[w] * cor_uv(t1_t, t2_t)[w] / 4.0;

          // Compute ac = m0*(m1+m2+m3+m4 - 2*(m5+m6+m7+m8) + 4*m9)
          ac(t1, t2)[w] =   m0  * ( m1(t1, t2)[w] + m2(t1, t2)[w] + m3(t1, t2)[w] + m4(t1, t2)[w]
                          - 2.0 * ( m5(t1, t2)[w] + m6(t1, t2)[w] + m7(t1, t2)[w] + m8(t1, t2)[w] )
                          + 4.0 *   m9(t1, t2)[w] );

          size_t min_pair_count = std::min(my_least_pair_count(t1_w, t2_w)[w],
                                           my_least_pair_count(t1_t, t2_t)[w]);

          my_extended_least_pair_count(t1, t2)[w] = std::min(my_extended_least_pair_count(t1, t2)[w], min_pair_count);

        }
      }
    }
  }

#if 0
  print_vec_matrix(my_least_pair_count, cout, "least_pair_count");
  print_vec_matrix(my_extended_least_pair_count, cout, "extended_least_pair_count");
  print_vec_matrix(ac, cout, "ac_matrix");
#endif
}

void
SubACCalM::module1_8(const Matrix2D< vector<internal_real_type> >& s1,
                     const Matrix2D< vector<internal_real_type> >& s2,
                           Matrix2D< vector<internal_real_type> >& result,
                     size_t i, size_t j, size_t p) const
{
  // Compute: result += s1 * s2 * my_weight_t[p][i] * my_weight_w[p][j]
  for(size_t t1 = 0, t1_t = 0, t2_t = 0; t1 < result.rows(); ++t1, ++t2_t)
  {
    if( t1 != 0 && t1 % trait_count() == 0 )
    {
      ++t1_t;
      t2_t = 0;
    }

    for(size_t t2 = 0, t1_w = 0, t2_w = 0; t2 < result.cols(); ++t2, ++t2_w)
    {
      if( t2 != 0 && t2 % trait_count() == 0 )
      {
        ++t1_w;
        t2_w = 0;
      }

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(result(t1, t2)[w])
            && !SAGE::isnan(s1(t1, t2)[w]) && !SAGE::isnan(s2(t1, t2)[w]) )
        {
          result(t1, t2)[w] +=   s1(t1, t2)[w] * my_weight_t[p][i](t1_w, t2_w)[w]
                               * s2(t1, t2)[w] * my_weight_w[p][j](t1_t, t2_t)[w];
        }
        else if( my_conservative )
          result(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void
SubACCalM::module9(const Matrix2D< vector<internal_real_type> >& s_ux,
                   const Matrix2D< vector<internal_real_type> >& s_vy,
                   const Matrix2D< vector<internal_real_type> >& s_uy,
                   const Matrix2D< vector<internal_real_type> >& s_vx,
                         Matrix2D< vector<internal_real_type> >& result,
                   size_t i, size_t j, size_t p) const
{
  // Compute: result += (s_ux*s_vy + s_uy*s_vx) * ti * wj
  for(size_t t1 = 0, t1_t = 0, t2_t = 0; t1 < result.rows(); ++t1, ++t2_t)
  {
    if( t1 != 0 && t1 % trait_count() == 0 )
    {
      ++t1_t;
      t2_t = 0;
    }

    for(size_t t2 = 0, t1_w = 0, t2_w = 0; t2 < result.cols(); ++t2, ++t2_w)
    {
      if( t2 != 0 && t2 % trait_count() == 0 )
      {
        ++t1_w;
        t2_w = 0;
      }

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(result(t1, t2)[w])                                          
            && !SAGE::isnan(s_ux(t1, t2)[w]) && !SAGE::isnan(s_vy(t1, t2)[w])
            && !SAGE::isnan(s_uy(t1, t2)[w]) && !SAGE::isnan(s_vx(t1, t2)[w]) )                  
        {
          result(t1, t2)[w] +=   ( s_ux(t1, t2)[w] * s_vy(t1, t2)[w] + s_uy(t1, t2)[w] * s_vx(t1, t2)[w] )
                               *   my_weight_t[p][i](t1_t, t2_t)[w]
                               *   my_weight_w[p][j](t1_w, t2_w)[w];
        }
        else if( my_conservative )                                        
          result(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();                        
      }
    }
  }
}

void
SubACCalM::sub_module(Matrix2D< vector<internal_real_type> >& sub,
                      size_t i, size_t j, size_t p, pair_type p_t)
{
#if 0
  cout << "before least_pair_count in sub_module() :" << endl;
  print_vec_matrix(my_least_pair_count, cout);
#endif

  vector<internal_real_type> init_result(my_weight_count, 0.0);

  sub.resize(0, 0);
  sub.resize(trait_count() * trait_count(), trait_count() * trait_count(), init_result);

  vector<double> init(my_weight_count, 0.0);

  Matrix2D< vector<double> > cor(trait_count(), trait_count(), init);

  for( size_t t1 = 0; t1 < cor.rows(); ++t1 )
  {
    for( size_t t2 = 0; t2 < cor.cols(); ++t2 )
    {
      vector< pair<double, size_t> > cor_pair = get_corinfo(my_pairset_uv, my_pairset_xy,
                                                            i, j, p, p_t,
                                                            t1 + trait_count(), t2);

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if( !SAGE::isnan(cor_pair[w].first) )
        {
          cor(t1, t2)[w] = cor_pair[w].first;

          my_least_pair_count(t1, t2)[w] = std::min(my_least_pair_count(t1, t2)[w], cor_pair[w].second);
        }
        else if( my_conservative )
          cor(t1, t2)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

#if 0
  cout << "after least_pair_count in sub_module() :" << endl;
  print_vec_matrix(my_least_pair_count, cout);
#endif
    
  const Matrix2D< vector<internal_real_type> >* ik = NULL;
  const Matrix2D< vector<internal_real_type> >* kj = NULL;
  const Matrix2D< vector<internal_real_type> >* kl = NULL;
  
  if     ( p_t == UX )
  {
    ik = &my_sum_ik_ux[p][i];
    kj = &my_sum_kj_ux[p][j];
    kl = &my_sum_kl_ux[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_ux_type(sub, cor, ik, kj, kl);

    return;
  }

  if( p_t == UY )
  {
    ik = &my_sum_ik_uy[p][i];
    kj = &my_sum_kj_uy[p][j];
    kl = &my_sum_kl_uy[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_uy_type(sub, cor, ik, kj, kl);

    return;
  }

  if( p_t == VX )
  {
    ik = &my_sum_ik_vx[p][i];
    kj = &my_sum_kj_vx[p][j];
    kl = &my_sum_kl_vx[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_vx_type(sub, cor, ik, kj, kl);

    return;
  } 

  if( p_t == VY )
  {
    ik = &my_sum_ik_vy[p][i];
    kj = &my_sum_kj_vy[p][j];
    kl = &my_sum_kl_vy[p];

    assert( ik != NULL && kj != NULL && kl != NULL );

    sub_vy_type(sub, cor, ik, kj, kl);

    return;
  }
}

void
SubACCalM::sub_ux_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1_sub = 0, t1_cor = 0; t1_sub < sub.rows(); ++t1_sub )
  {
    if( t1_sub != 0 && t1_sub % trait_count() == 0 )
      ++t1_cor;

    for( size_t t2_sub = 0, t2_cor = 0; t2_sub < sub.cols(); ++t2_sub )
    {
      if( t2_sub != 0 && t2_sub % trait_count() == 0 )
        ++t2_cor;

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(  cor(t1_cor, t2_cor)[w]) && !SAGE::isnan((*ik)(t1_cor, t2_cor)[w])
            && !SAGE::isnan((*kj)(t1_cor, t2_cor)[w]) && !SAGE::isnan((*kl)(t1_cor, t2_cor)[w]) )
        {
          sub(t1_sub, t2_sub)[w] =     cor(t1_cor, t2_cor)[w] - (*ik)(t1_cor, t2_cor)[w]
                                   - (*kj)(t1_cor, t2_cor)[w] + (*kl)(t1_cor, t2_cor)[w];

          my_extended_least_pair_count(t1_sub, t2_sub)[w] = std::min(my_extended_least_pair_count(t1_sub, t2_sub)[w],
                                                                     my_least_pair_count(t1_cor, t2_cor)[w]);
        }
        else if( my_conservative )
          sub(t1_sub, t2_sub)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void
SubACCalM::sub_uy_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1_sub = 0, t1_cor = 0; t1_sub < sub.rows(); ++t1_sub )
  {
    if( t1_sub != 0 && t1_sub % trait_count() == 0 )
      ++t1_cor;

    for( size_t t2_sub = 0, t2_cor = 0; t2_sub < sub.cols(); ++t2_sub, ++t2_cor )
    {
      if( t2_cor == trait_count() )
        t2_cor = 0;

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(  cor(t1_cor, t2_cor)[w]) && !SAGE::isnan((*ik)(t1_cor, t2_cor)[w])
            && !SAGE::isnan((*kj)(t1_cor, t2_cor)[w]) && !SAGE::isnan((*kl)(t1_cor, t2_cor)[w]) )
        {
          sub(t1_sub, t2_sub)[w] =     cor(t1_cor, t2_cor)[w] - (*ik)(t1_cor, t2_cor)[w]
                                   - (*kj)(t1_cor, t2_cor)[w] + (*kl)(t1_cor, t2_cor)[w];

          my_extended_least_pair_count(t1_sub, t2_sub)[w] = std::min(my_extended_least_pair_count(t1_sub, t2_sub)[w],
                                                                     my_least_pair_count(t1_cor, t2_cor)[w]);
        }
        else if( my_conservative )
          sub(t1_sub, t2_sub)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void
SubACCalM::sub_vx_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1_sub = 0, t1_cor = 0; t1_sub < sub.rows(); ++t1_sub, ++t1_cor )
  {
    if( t1_cor == trait_count() )
      t1_cor = 0;

    for( size_t t2_sub = 0, t2_cor = 0; t2_sub < sub.cols(); ++t2_sub)
    {
      if( t2_sub != 0 && t2_sub % trait_count() == 0 )
        ++t2_cor;

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(  cor(t1_cor, t2_cor)[w]) && !SAGE::isnan((*ik)(t1_cor, t2_cor)[w])
            && !SAGE::isnan((*kj)(t1_cor, t2_cor)[w]) && !SAGE::isnan((*kl)(t1_cor, t2_cor)[w]) )
        {
          sub(t1_sub, t2_sub)[w] =     cor(t1_cor, t2_cor)[w] - (*ik)(t1_cor, t2_cor)[w]
                                   - (*kj)(t1_cor, t2_cor)[w] + (*kl)(t1_cor, t2_cor)[w];

          my_extended_least_pair_count(t1_sub, t2_sub)[w] = std::min(my_extended_least_pair_count(t1_sub, t2_sub)[w],
                                                                     my_least_pair_count(t1_cor, t2_cor)[w]);
        }
        else if( my_conservative )
          sub(t1_sub, t2_sub)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

void
SubACCalM::sub_vy_type(      Matrix2D< vector<internal_real_type> >& sub,
                       const Matrix2D< vector<double> >&             cor,
                       const Matrix2D< vector<internal_real_type> >* ik,
                       const Matrix2D< vector<internal_real_type> >* kj,
                       const Matrix2D< vector<internal_real_type> >* kl  )
{
  for( size_t t1_sub = 0, t1_cor = 0; t1_sub < sub.rows(); ++t1_sub, ++t1_cor )
  {
    if( t1_cor == trait_count() )
      t1_cor = 0;

    for( size_t t2_sub = 0, t2_cor = 0; t2_sub < sub.cols(); ++t2_sub, ++t2_cor )
    {
      if( t2_cor == trait_count() )
        t2_cor = 0;

      for( size_t w = 0; w < my_weight_count; ++w )
      {
        if(    !SAGE::isnan(  cor(t1_cor, t2_cor)[w]) && !SAGE::isnan((*ik)(t1_cor, t2_cor)[w])
            && !SAGE::isnan((*kj)(t1_cor, t2_cor)[w]) && !SAGE::isnan((*kl)(t1_cor, t2_cor)[w]) )
        {
          sub(t1_sub, t2_sub)[w] =     cor(t1_cor, t2_cor)[w] - (*ik)(t1_cor, t2_cor)[w]
                                   - (*kj)(t1_cor, t2_cor)[w] + (*kl)(t1_cor, t2_cor)[w];

          my_extended_least_pair_count(t1_sub, t2_sub)[w] = std::min(my_extended_least_pair_count(t1_sub, t2_sub)[w],
                                                                     my_least_pair_count(t1_cor, t2_cor)[w]);
        }
        else if( my_conservative )
          sub(t1_sub, t2_sub)[w] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
}

// end of SubACCalM Implementation

} // end of namespace FCOR
} // end of namespace SAGE
