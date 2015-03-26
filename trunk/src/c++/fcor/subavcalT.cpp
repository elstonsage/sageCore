//****************************************************************************
//* File:      subavcalT.cpp                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This source file calculates the asymptotic variance of        *
//*            a multipedigree correlation(XY & XY).                         *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/subavcalT.h"

namespace SAGE {
namespace FCOR {
    
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubAVCalT                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic variance of a correlation(XY & XY).     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of SubAVCalT
// ---------------------------------------------------------------------------

SubAVCalT::SubAVCalT(const SubPairSetData* p, const SubCorrelationCal* c)
{ 
  my_pairsetdata = p;
  my_correlation = c;
}

void
SubAVCalT::compute_asymptotic_variance(const pairset_by_pedigree_type& pairset,
                                       const CorrelationInfo&          xy,
                                       TriangleMatrix<double>&         av,
                                       TriangleMatrix<size_t>&         replaced)
{
  if( my_pairsetdata == NULL || my_correlation == NULL )
    return;

  my_corinfo_xy   = &xy;
  my_pairset_xy   = &pairset;

  my_conservative = true;

  init();
  compute(av);

  bool parser_conservative = my_pairsetdata->fcor_parser()->conservative();

  if( parser_conservative )
    return;

  my_conservative = false;

  TriangleMatrix<double> av_replaced(av.size(), std::numeric_limits<double>::quiet_NaN());

  init();
  compute(av_replaced);
  

  for( size_t t = 0; t < av.linear_size(); ++t )
    if( av(t) != av_replaced(t) )
    {
      av(t) = av_replaced(t);
      replaced(t) = 1;
    }

  return;
}

void
SubAVCalT::init()
{
  compute_we();

  my_weight_w.resize(0);
  my_weight_w.resize(pedigree_count());

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    my_weight_w[p].resize(0);
    my_weight_w[p].resize(xy_pair_count(p));

    for( size_t i = 0; i < xy_pair_count(p); ++i )
      build_weight_triangle(my_weight_w[p][i], i, p);
  }

  my_sum_kl.resize(0);
  my_sum_ikj.resize(0);

  my_sum_kl.resize(pedigree_count());
  my_sum_ikj.resize(pedigree_count());

  TriangleMatrix<internal_real_type> sum_p(trait_count(), 0.0);

  for( size_t p = 0; p < pedigree_count(); ++p )
  {
    sum_p.fill(0.0);
    pedigree_weight(sum_p, p);

    my_sum_ikj[p].resize(xy_pair_count(p));

    compute_sub(sum_p, my_sum_ikj[p], my_sum_kl[p], p);
  }
}

void
SubAVCalT::compute_we()
{
  my_we.resize(trait_count(), 0.0);

  // Compute 1/(W - sum(w^2)/W)
  for( size_t i = 0; i < my_we.linear_size(); ++i )
  {
    double cor            = my_corinfo_xy->correlation(i);
    internal_real_type w  = my_corinfo_xy->sum_weight(i);
    internal_real_type ww = my_corinfo_xy->sum_weight_square(i);

    if( !SAGE::isnan(cor) && w > std::numeric_limits<double>::epsilon() )
      my_we[i] = 1.0 / (w - ww / w);  
  }
}

void
SubAVCalT::compute_sub(const TriangleMatrix<internal_real_type>& sum_p,
                     vector< TriangleMatrix<internal_real_type> >& ikj,
                             TriangleMatrix<internal_real_type>&   kl, 
                       size_t p)
{
  kl.fill(0);
  kl.resize(trait_count(), 0.0);

  TriangleMatrix<internal_real_type> local_sum_w (trait_count(), 0.0);
  TriangleMatrix<internal_real_type> local_sum_ww(trait_count(), 0.0);
    
  for( size_t i = 0; i < xy_pair_count(p); ++i )
  {
    ikj[i].fill(0);
    ikj[i].resize(trait_count(), 0.0);

    local_sum_w.fill(0);
    local_sum_w.resize(trait_count(), 0.0);

    for( size_t j = 0; j < xy_pair_count(p); ++j )
    {
//      const CorrelationInfo& c = find_class(i, j, p);

//      if( c == my_correlation->corinfo().end() )
//        continue;

      for( size_t t1 = 0; t1 < ikj[i].size(); ++t1 )
        for( size_t t2 = 0; t2 <= t1; ++t2 )
        {
          double cor = get_corinfo(i, j, p, t1 + trait_count(), t2);

          if( !SAGE::isnan(cor) && !SAGE::isnan(kl(t1, t2)) && !SAGE::isnan(ikj[i](t1, t2)) )
          {
            // Compute:
            //  wi  = my_weight_w[p][i]
            //  wj  = my_weight_w[p][j]
            //  kl += wi*wj*kl
            kl(t1, t2) += my_weight_w[p][i](t1, t2) * my_weight_w[p][j](t1, t2) * cor;

            // Compute:
            //  ikj   += my_weight_w[p][j] * ikj
            ikj[i](t1, t2) += my_weight_w[p][j](t1, t2) * cor;

            local_sum_w(t1, t2)  += my_weight_w[p][j](t1, t2);
            local_sum_ww(t1, t2) += my_weight_w[p][i](t1, t2) * my_weight_w[p][j](t1, t2);
          }
          else if( my_conservative )
          {
            kl(t1, t2)     = std::numeric_limits<double>::quiet_NaN();
            ikj[i](t1, t2) = std::numeric_limits<double>::quiet_NaN();
          }
        }
    }

    for( size_t m = 0; m < ikj[i].linear_size(); ++m )
    {
      local_sum_w[m] = my_corinfo_xy->sum_weight(m) - sum_p[m] + local_sum_w[m];

      if( !SAGE::isnan(ikj[i][m]) && local_sum_w[m] > std::numeric_limits<double>::epsilon() )
        ikj[i][m] /= local_sum_w[m];
    }
  }

  for( size_t m = 0; m < kl.linear_size(); ++m )
  {
    internal_real_type sum_ww = my_corinfo_xy->sum_weight(m) * my_corinfo_xy->sum_weight(m);
    internal_real_type sum_pp = sum_p[m] * sum_p[m];

    local_sum_ww[m] = sum_ww - sum_pp + local_sum_ww[m];

    if( !SAGE::isnan(kl[m]) && local_sum_ww[m] > std::numeric_limits<double>::epsilon() )
      kl[m] /= local_sum_ww[m];
  }
} 

void
SubAVCalT::compute(TriangleMatrix<double>& ac)
{
  TriangleMatrix<internal_real_type> m1 (trait_count(), 0.0);
  TriangleMatrix<internal_real_type> m23(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> m4 (trait_count(), 0.0);
  TriangleMatrix<internal_real_type> m57(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> m68(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> m9 (trait_count(), 0.0);

  TriangleMatrix<internal_real_type> sub_xx(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> sub_xy(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> sub_yx(trait_count(), 0.0);
  TriangleMatrix<internal_real_type> sub_yy(trait_count(), 0.0);

  for( size_t p = 0; p < pedigree_count(); ++p )
    for( size_t i = 0; i < xy_pair_count(p); ++i )
      for( size_t j = 0; j < xy_pair_count(p); ++j )
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

  TriangleMatrix<internal_real_type> cov(trait_count());
  TriangleMatrix<internal_real_type> cov_mod(trait_count());

  for( size_t t1 = 0; t1 < m1.size(); ++t1 )
    for( size_t t2 = 0; t2 <= t1; ++t2 )
    {
      cov(t1, t2) = my_corinfo_xy->covariance(t1, t2);

      double xx = my_corinfo_xy->covariance(t1, t1);
      double yy = my_corinfo_xy->covariance(t2, t2);

      if( !SAGE::isnan(xx) && !SAGE::isnan(yy) )
        cov_mod(t1, t2) = sqrt(xx * yy);
      else if( my_conservative )
        cov_mod(t1, t2) = std::numeric_limits<double>::quiet_NaN();
    }

  for( size_t k = 0; k < m1.linear_size(); ++k )
  {
    double m0k = my_corinfo_xy->correlation(k);

    if(    !SAGE::isnan(m1[k])      && !SAGE::isnan(m23[k]) && !SAGE::isnan(m4[k])
        && !SAGE::isnan(m57[k])     && !SAGE::isnan(m68[k]) && !SAGE::isnan(m9[k])
        && !SAGE::isnan(cov_mod[k]) && !SAGE::isnan(m0k)
        && fabs(cov[k]) > std::numeric_limits<double>::epsilon()  )
    {
      m1 [k] *= 2.0 * my_we[k] * my_we[k];
      m23[k] *= 2.0 * my_we[k] * my_we[k];
      m4 [k] *= 2.0 * my_we[k] * my_we[k];
      m57[k] *= 2.0 * my_we[k] * my_we[k] * cov_mod[k];
      m68[k] *= 2.0 * my_we[k] * my_we[k] * cov_mod[k];
      m9 [k] *=       my_we[k] * my_we[k] * cov_mod[k] * cov_mod[k];

      // Compute: m57, m68, m9 /= cov
      m57[k] /= cov[k];
      m68[k] /= cov[k];
      m9[k]  /= ( cov[k] * cov[k] );

      m0k  *= m0k / 4.0;

      // Compute ac = m0*(m1+m23+m23+m4 - 2*(m57+m68+m57+m68) + 4*m9)
      ac[k] =   m0k * ( m1[k]  + m23[k] + m23[k] + m4[k]
              - 2.0 * ( m57[k] + m68[k] + m57[k] + m68[k] )
              + 4.0 *   m9[k] );
    }
  }
}

void
SubAVCalT::module1_8(const TriangleMatrix<internal_real_type>& s1,
                     const TriangleMatrix<internal_real_type>& s2,
                           TriangleMatrix<internal_real_type>& result,
                     size_t i, size_t j, size_t p) const
{
  // Compute: result += s1 * s2 * my_weight_w[p][i] * my_weight_w[p][j]
  for( size_t m = 0; m < result.linear_size(); ++m )
  {
    if(    !SAGE::isnan(result[m])
        && !SAGE::isnan(s1[m]) && !SAGE::isnan(s2[m]) )
      result[m] += s1[m] * s2[m] * my_weight_w[p][i][m] * my_weight_w[p][j][m];
    else if( my_conservative )
      result[m] = std::numeric_limits<double>::quiet_NaN();
  }
}

void
SubAVCalT::module9(const TriangleMatrix<internal_real_type>& s_xx,
                   const TriangleMatrix<internal_real_type>& s_yy, 
                   const TriangleMatrix<internal_real_type>& s_xy,
                   const TriangleMatrix<internal_real_type>& s_yx, 
                         TriangleMatrix<internal_real_type>& result,
                   size_t i, size_t j, size_t p) const
{
  // Compute: result += (s_xx*s_yy + s_xy*s_yx) * wi * wj
  for( size_t m = 0; m < result.linear_size(); ++m )
  {
    if(    !SAGE::isnan(result[m])
        && !SAGE::isnan(s_xx[m])  && !SAGE::isnan(s_yy[m])
        && !SAGE::isnan(s_xy[m])  && !SAGE::isnan(s_yx[m]) )
      result[m] +=   ( s_xx[m] * s_yy[m] + s_xy[m] * s_yx[m] )
                   *   my_weight_w[p][i][m]
                   *   my_weight_w[p][j][m];
    else if( my_conservative )
      result[m] = std::numeric_limits<double>::quiet_NaN();
  }
}

void
SubAVCalT::sub_module(TriangleMatrix<internal_real_type>& sub,
                      size_t i, size_t j, size_t p, pair_type p_t) const
{
  sub.fill(0.0);
  sub.resize(trait_count(), 0.0);

  TriangleMatrix<double> cor(trait_count(), 0.0);

//  const CorrelationInfo& c = find_class(i, j, p);

//  if( c != my_correlation->corinfo().end() )
    for( size_t t1 = 0; t1 < trait_count(); ++t1 )
      for( size_t t2 = 0; t2 <= t1; ++t2 )
        cor(t1, t2) = get_corinfo(i, j, p, t1 + trait_count(), t2);

  const TriangleMatrix<internal_real_type>* ik = &my_sum_ikj[p][i];
  const TriangleMatrix<internal_real_type>* kj = &my_sum_ikj[p][j];
  const TriangleMatrix<internal_real_type>* kl = &my_sum_kl[p];
  
  if( p_t == XX )
  {
    for( size_t t1 = 0; t1 < cor.size(); ++t1 )
     for( size_t t2 = 0; t2 <= t1; ++t2 )
       if(    !SAGE::isnan(cor(t1, t1))   && !SAGE::isnan((*ik)(t1, t1))
           && !SAGE::isnan((*kj)(t1, t1)) && !SAGE::isnan((*kl)(t1, t1)) )
         sub(t1, t2) =    cor (t1, t1) - (*ik)(t1, t1)
                       - (*kj)(t1, t1) + (*kl)(t1, t1);
       else if( my_conservative )
         sub(t1, t2) = std::numeric_limits<double>::quiet_NaN();    

    return;
  }

  if( p_t == XY )
  {
    for( size_t t1 = 0; t1 < cor.size(); ++t1 )
     for( size_t t2 = 0; t2 <= t1; ++t2 )
       if(    !SAGE::isnan(cor(t1, t2))   && !SAGE::isnan((*ik)(t1, t2))
           && !SAGE::isnan((*kj)(t1, t2)) && !SAGE::isnan((*kl)(t1, t2)) )
         sub(t1, t2) =    cor (t1, t2) - (*ik)(t1, t2)
                       - (*kj)(t1, t2) + (*kl)(t1, t2);
       else if( my_conservative )
         sub(t1, t2) = std::numeric_limits<double>::quiet_NaN();    

    return;
  }

  if( p_t == YX )
  {
    for( size_t t1 = 0; t1 < cor.size(); ++t1 )
     for( size_t t2 = 0; t2 <= t1; ++t2 )
       if(    !SAGE::isnan(cor(t2, t1))   && !SAGE::isnan((*ik)(t2, t1))
           && !SAGE::isnan((*kj)(t2, t1)) && !SAGE::isnan((*kl)(t2, t1)) )
         sub(t1, t2) =    cor (t2, t1) - (*ik)(t2, t1)
                       - (*kj)(t2, t1) + (*kl)(t2, t1);
       else if( my_conservative )
         sub(t1, t2) = std::numeric_limits<double>::quiet_NaN();    

    return;
 } 

  if( p_t == YY )
  {
    for( size_t t1 = 0; t1 < cor.size(); ++t1 )
     for( size_t t2 = 0; t2 <= t1; ++t2 )
       if(    !SAGE::isnan(cor(t2, t2))   && !SAGE::isnan((*ik)(t2, t2))
           && !SAGE::isnan((*kj)(t2, t2)) && !SAGE::isnan((*kl)(t2, t2)) )
         sub(t1, t2) =    cor (t2, t2) - (*ik)(t2, t2)
                       - (*kj)(t2, t2) + (*kl)(t2, t2);
       else if( my_conservative )
         sub(t1, t2) = std::numeric_limits<double>::quiet_NaN();    

    return;
  }
}

void
SubAVCalT::build_weight_triangle(TriangleMatrix<double>& weight, size_t i, size_t p)
{
  const pedigree_member_type* member1;

  member1 = (*my_pairset_xy)[p][i].member_pair.first;

  weight.resize(trait_count(), (*my_pairset_xy)[p][i].pair_weight);

  for( size_t t1 = 0; t1 < trait_count(); ++t1 )
  {
    if( SAGE::isnan( member1->pedigree()->info().trait(member1->index(), t1) ) )
      for( size_t t2 = 0; t2 <= t1; ++t2 )
        weight(t1, t2) = 0.0;
  }
}

// end of SubAVCalT Implementation

} // end of namespace FCOR
} // end of namespace SAGE
