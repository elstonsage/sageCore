//=============================================================================
// File:    sib_matrix.h
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//          Added full-half combined weight build routine.         yjs  Jun.09
//
// Notes:   This file contains definition for following data structures.
//
//            class SibpalWeights
//            
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/sib_matrix.h"

namespace SAGE   {
namespace SIBPAL {

// This is only a stop-gap.  Testing rank requires looking at the actual
// pairs since we may drop to a lower, but valid number of sibs and yet
// still not be of full rank.
bool full_sib_rank(size_t s)
{
  return false;
#if 0
  size_t n = rint( sqrt(1.0+8.0*s)+1.0) / 2.0;

  size_t p = n*(n-1)/2;

  return (s == p);
#endif
}

void
sib_weights::weight_matrix(const sib_cluster& ship, matrix& W,
                           double p2, double p1, double p0, 
                           weight_status_type& status)
{
  weight_matrix(ship, W, sib_matrix_pattern(p2,p1,p0,0.), status);
}

void
sib_weights::weight_matrix(const sib_cluster& ship, matrix& W,
                           const sib_matrix_pattern& p,
                           weight_status_type& status)
{
  size_t n = ship.valid_pair_count();

  if( status == BESTW )
    status = NORMALW;

  if(!n)
  {
    W.setstate(matrix::badbit);
    return;
  }

  W.clear();
  W.resize_nofill(n,n);

  sib_matrix_pattern pp;

  bool full = full_sib_rank(n);

//  if(full && status == PALBASE::INVERSEW)
//  {
//    int s = ship.valid_sib_count();
//    pp = magic_inverse(s, p);
//  }
//  else
    pp = p;

  for( size_t i = 0; i < n; ++i )
    W(i, i) = pp[2];

  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = i+1; j < n; ++j )
    {
      int c = sibs_shared( ship[i].rels().pair, ship[j].rels().pair );
      W(i, j) = W(j, i) = pp[c];
    }
  }

  if( !full && status == INVERSEW )
  {
    matrix W2 = W;

    SVDinverse(W2,W);
  }
}

void
sib_weights::weight_matrix_combined(const sib_cluster& ship, matrix& W, double c,
                                    double p2_ff, double p1_ff, double p0_ff,
                                    double p2_hh, double p1_hh, double p0_hh,
                                    double p2_fh, double p1_fh, double p0_fh, 
                                    weight_status_type& status)
{
  sib_matrix_pattern p_ff(p2_ff, p1_ff, p0_ff, 0.);
  sib_matrix_pattern p_hh(p2_hh, p1_hh, p0_hh, 0.);
  sib_matrix_pattern p_fh(p2_fh, p1_fh, p0_fh, 0.);

  weight_matrix_combined(ship, W, c, p_ff, p_hh, p_fh, status);
}

void
sib_weights::weight_matrix_combined(const sib_cluster& ship, matrix& W, double c,
                                    const sib_matrix_pattern& p_ff,
                                    const sib_matrix_pattern& p_hh,
                                    const sib_matrix_pattern& p_fh,
                                    weight_status_type& status)
{
  size_t n = ship.valid_pair_count();

  if( status == BESTW )
    status = NORMALW;

  if(!n)
  {
    W.setstate(matrix::badbit);
    return;
  }

  W.clear();
  W.resize_fill(n, n, 0.0);

  for( size_t i = 0; i < n; ++i )
  {
    if( ship[i].is_fsib_pair() )
    {
      W(i, i) = p_ff[2];
    }
    else if( ship[i].is_hsib_pair() )
    {
      W(i, i) = c * p_hh[2];
    }
  }

  for( size_t i = 0; i < n; ++i )
  {
    for( size_t j = i+1; j < n; ++j )
    {
      switch( sibs_shared( ship[i].rels().pair, ship[j].rels().pair ) )
      {
        case 2:  break;  // Dealt already above

        case 1:  // One sib in common

          if( ship[i].is_fsib_pair() && ship[j].is_fsib_pair() )
          {
            W(i, j) = W(j, i) = p_ff[1];
          }
          else if( ship[i].is_hsib_pair() && ship[j].is_hsib_pair() )
          {
            W(i, j) = W(j, i) = p_hh[1];
          }
          else
          {
            W(i, j) = W(j, i) = p_fh[1];
          }

          break;

        case 0:  // No sibs in common

          if( ship[i].is_fsib_pair() && ship[j].is_fsib_pair() )
          {
            W(i, j) = W(j, i) = p_ff[0];
          }
          else if( ship[i].is_hsib_pair() && ship[j].is_hsib_pair() )
          {
            W(i, j) = W(j, i) = p_hh[0];
          }
          else
          {
            W(i, j) = W(j, i) = p_fh[0];
          }

          break;

        default: // Unrelated pairs = 0.0
          break;
      }
    }
  }

  if( status == INVERSEW )
  {
    matrix W2 = W;

    SVDinverse(W2, W);
  }

  return;
}


void
sib_weights::normal_diff_weight_matrix(const sib_cluster &ship, matrix &W, double p,
                                       weight_status_type& status)
{
  weight_matrix(ship, W, 1, p, 0, status);
}

void
sib_weights::normal_sum_weight_matrix(const sib_cluster &ship, matrix &W, double r,
                                      weight_status_type& status)
{
  double p1 = 0.5*(1.0+3*r)/(1.0+r);
  double p0 = 2.0*r/(1.0+r);
  p1 = p1*p1;
  p0 = p0*p0;

  weight_matrix(ship, W, 1, p1, p0, status);
}

void
sib_weights::normal_prod_weight_matrix(const sib_cluster &ship, matrix &W, double r,
                                       weight_status_type& status)
{
  double p1 = r*(1.0+r)/(1.0+r*r);
  double p0 = 2.0*r*r/(1.0+r*r);

  weight_matrix(ship, W, 1, p1, p0, status);
}
/*
void
sib_weights::weighted_sum_matrix(const sib_cluster& ship, matrix& W,
                                 const sib_matrix_pattern& p1,
                                 double ss1,
                                 const sib_matrix_pattern& p2,
                                 double ss2,
                                 weight_status_type& status)
{
  // Given pattern weight matricies W_p1 and W_p2 compute:
  // M^-1 = (1/ss1 * W_p1^-1 + 1/ss2 * W_p2^-1)
  // Return M      if not inverse
  //        M^-1   otherwise

  int s = ship.valid_sib_count();

  if( full_sib_rank(s) )
    weighted_sum_matrix_full_rank(ship, W, p1, ss1, p2, ss2, status);

  if( status == PALBASE::BESTW )
    status = PALBASE::INVERSEW;

  weight_status_type istatus = PALBASE::INVERSEW;
  weight_matrix(ship, W,  p1, istatus);
  weight_matrix(ship, W2, p2, istatus);

  W  /= ss1;
  W2 /= ss2;

  if(status == PALBASE::INVERSEW)
  {
    W += W2;
  }
  else
  {
    W2 += W;

#if TRY_SYMINVERSE
    SYMinverse(W2,W);
    if(!W)
#endif
      SVDinverse(W2,W);
  }
}

void
sib_weights::weighted_sum_matrix_full_rank(const sib_cluster& ship, matrix& W,
                                           const sib_matrix_pattern& p1,
                                           double ss1,
                                           const sib_matrix_pattern& p2,
                                           double ss2,
                                           weight_status_type& status)
{
  // Given pattern weight matricies W_p1 and W_p2 compute:
  // M^-1 = (1/ss1 * W_p1^-1 + 1/ss2 * W_p2^-1)
  // Return M      if not inverse
  //        M^-1   otherwise

  int s = ship.valid_sib_count();


  // Construct seperate inverses
  sib_matrix_pattern pi1 = magic_inverse(s, p1);
  sib_matrix_pattern pi2 = magic_inverse(s, p2);

  for(size_t i = 0; i < 3; ++i)
  {
    pi1[i] /= ss1;
    pi2[i] /= ss2;
  }

  sib_matrix_pattern pi( pi1[2]+pi2[2], pi1[1]+pi2[1], pi1[0]+pi2[0] );

  weight_status_type reverse_status;

  if( status == PALBASE::BESTW )
    status = PALBASE::INVERSEW;

  if( status == PALBASE::NORMALW )
    reverse_status = PALBASE::INVERSEW;
  else
    reverse_status = PALBASE::NORMALW;

  // Construct final weight matrix
  weight_matrix(ship, W, pi, reverse_status);
}
*/

void marker_matrix(const sib_cluster &ship, matrix &y, size_t m, double w0, double w1, double w2)
{
  size_t pair_count = ship.valid_pair_count();

  if( !pair_count )
  {
    y.setstate(matrix::failbit);
    return;
  }

  y.clear();
  y.resize_fill(pair_count,2,0);

  for(size_t j = 0; j < pair_count; ++j)
  {
    y(j,0)=ship[j].weighted_share(m, w0, w1, w2);
    y(j,1)=ship[j].prob_share(m,1);

#ifdef DEBUG_MARKER_VECTOR
    cout << ship[j].sibs().first->pedigree()->name() << ":"
         << ship[j].sibs().first->name() << " = " << y(j,0) << endl;
    cout << ship[j].sibs().first->pedigree()->name() << ":"
         << ship[j].sibs().second->name() << " = " << y(j,1) << endl;
#endif
  }
}

/*
sib_matrix_pattern
magic_inverse(size_t s, const sib_matrix_pattern& p)
{
  // Given a symmetric pattern matrix X of rank equal s*(s-1)^2 with rows
  // and columns representing pairs (1..s,1..s) where no pair is repeated
  // and no self pairs.  The entries of X are determined by how many
  // elements of the row and column pair are the same c={0, 1, 2}.  A 
  // vector of reals p[0..2] provides the entries in each of the cells:
  //
  //  X_ii = p[2]
  //  X_ij = p[1] iff i and j represent a pairs with one index in common
  //  X_ij = p[0] iff i and j represent a pairs with no index in common
  //
  // This function returns the pattern elements of the inverse of X in
  // constant time by solving a set of linear equations by Cramer's rule.
  // If X is ill-conditioned, then a vector of qNaN is returned.

  double d11 = p[2];
  double d21 =                       p[1];
  double d31 =                                                         p[0];
  double d12 =           (2.0*s-4.0)*p[1];
  double d22 = p[2] +        (s-2.0)*p[1] +                    (s-3.0)*p[0];
  double d32 =                   4.0*p[1] +                (2.0*s-8.0)*p[0];
  double d13 =                               (0.5*(s-4.0)*(s-1.0)+1.0)*p[0];
  double d23 =               (s-3.0)*p[1] +        0.5*(s-3.0)*(s-4.0)*p[0];
  double d33 = p[2] +    (2.0*s-8.0)*p[1] +        0.5*(s-5.0)*(s-4.0)*p[0];

  double d   = d11*d22*d33+d12*d23*d31+d13*d21*d32
             -(d13*d22*d31+d12*d21*d33+d11*d23*d32);

  if( fabs(d) <= 0 )
    return sib_matrix_pattern();

  double a =  (d22*d33-d23*d32)/d;
  double b =  (d23*d31-d21*d33)/d;
  double c =  (d21*d32-d22*d31)/d;

  return sib_matrix_pattern(a,b,c);
}
*/

} // end of namespace SIBPAL
} // end of namespace SAGE
