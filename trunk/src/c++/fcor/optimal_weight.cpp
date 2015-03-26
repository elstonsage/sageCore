//****************************************************************************
//* File:      optimal_weight.cpp                                            *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Mar 04 *
//*                                                                          *
//* Notes:     This file implements class for finding the optimal weight.    *
//*                                                                          *
//* Copyright (c) 2004 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/optimal_weight.h"

namespace SAGE {
namespace FCOR {

// ---------------------------------------------------------------------------
// Out-of-Line Implementation of optimal_weight_finder
// ---------------------------------------------------------------------------

optimal_weight_finder::optimal_weight_finder()
{}

void
optimal_weight_finder::estimate_variance_by_quadradic(pairset_result& result) const
{
  // Estimate Cor & Var by quadradic
  //
  for( size_t i = 0; i < result.corr.rows(); ++i )
  {
    for( size_t j = 0; j < result.corr.cols(); ++j )
    {
      double V0 = result.std_err(i, j).standard_error[UNIFORM];
      double V2 = result.std_err(i, j).standard_error[PAIR_WISE];
      double V1 = result.std_err(i, j).standard_error[MEAN];

      double cor0 = result.corr(i, j).correlation[UNIFORM];
      double cor2 = result.corr(i, j).correlation[PAIR_WISE];
      double cor1 = result.corr(i, j).correlation[MEAN];

      result.corr(i, j).correlation[WEIGHT_COUNT]       = cor0; //1.5;
      result.std_err(i, j).standard_error[WEIGHT_COUNT] = V0;   //0.0;

      //result.corr(i, j).correlation[WEIGHT_COUNT + QUAD_SV]       = cor2; //1.5;
      //result.std_err(i, j).standard_error[WEIGHT_COUNT + QUAD_SV] = V2;   //0.0;

      result.std_err(i, j).w1 = 0.0;

      if( V1 < V0 && V1 < V2 )
      {
        double b0 =      V0;
        double b1 = -3.0*V0 + 4.0*V1 -     V2;
        double b2 =  2.0*V0 - 4.0*V1 + 2.0*V2;

        double Vquad = b0 - ((b1*b1) / (4.0*b2)) ;

        double w1 = -b1 / (2.0*b2);

        double a0 =      cor0;
        double a1 = -3.0*cor0 + 4.0*cor1 -     cor2;
        double a2 =  2.0*cor0 - 4.0*cor1 + 2.0*cor2;

        double Cquad = a0 + a1*w1 + a2*(w1*w1) ;

        result.corr(i, j).correlation[WEIGHT_COUNT]       = Cquad;
        result.std_err(i, j).standard_error[WEIGHT_COUNT] = Vquad;

        //result.corr(i, j).correlation[WEIGHT_COUNT + QUAD_SV]       = Cquad;
        //result.std_err(i, j).standard_error[WEIGHT_COUNT + QUAD_SV] = Vquad;

        result.std_err(i, j).w1 = w1;
      }
      else if( V0 < V2 )
      {
        result.corr(i, j).correlation[WEIGHT_COUNT]       = cor2; //-1.5;
        result.std_err(i, j).standard_error[WEIGHT_COUNT] = V2;   // 0.0;

        //result.corr(i, j).correlation[WEIGHT_COUNT + QUAD_SV]       = cor0; //-1.5;
        //result.std_err(i, j).standard_error[WEIGHT_COUNT + QUAD_SV] = V0;   // 0.0;

        result.std_err(i, j).w1 = 1.0;
      }
    }
  }
}

void
optimal_weight_finder::estimate_weighted_covariance(size_t                            trait_count,
                                                    const pairset_result&             result1,
                                                    const pairset_result&             result2,
                                                    const Matrix2D< vector<double> >& ac_vec_matrix,
                                                          Matrix2D< double >&         ac_matrix) const
{
#if 0
  cout << result1.std_err.rows() << "," << result1.std_err.cols() << " "
       << result2.std_err.rows() << "," << result2.std_err.cols() << endl;
#endif
  for( size_t r1t1 = 0, r = 0; r1t1 < trait_count; ++r1t1 )
  {
    for( size_t r1t2 = 0; r1t2 < trait_count; ++r1t2, ++r )
    {
      double w_x  = result1.std_err(r1t1, r1t2).w1;
      double v_xp = result1.std_err(r1t1, r1t2).standard_error[PAIR_WISE];
      double v_xu = result1.std_err(r1t1, r1t2).standard_error[UNIFORM];

      for( size_t r2t1 = 0, c = 0; r2t1 < trait_count; ++r2t1 )
      {
        for( size_t r2t2 = 0; r2t2 < trait_count; ++r2t2, ++c )
        {
          double w_y  = result2.std_err(r2t1, r2t2).w1;
          double v_yp = result2.std_err(r2t1, r2t2).standard_error[PAIR_WISE];
          double v_yu = result2.std_err(r2t1, r2t2).standard_error[UNIFORM];

          double cov_xpyp = ac_vec_matrix(r, c)[PAIR_WISE];
          double cov_xuyu = ac_vec_matrix(r, c)[UNIFORM];

          double cor_xpyp = cov_xpyp / sqrt(v_xp * v_yp);
          double cor_xuyu = cov_xuyu / sqrt(v_xu * v_yu);

          double r_xy = 0.5 * (cor_xpyp + cor_xuyu);

          double cov_xpyu = r_xy * sqrt(v_xp * v_yu);
          double cov_xuyp = r_xy * sqrt(v_xu * v_yp);

          double cov_xy = 0.0;

          if( w_x == 0.0 ) // w_x is unform
          {
            if( w_y == 0.0 ) // w_y is uniform
            {
              cov_xy = cov_xuyu;
            }
            else if( w_y == 1.0 ) // w_y is pair-wise
            {
              cov_xy = cov_xuyp;
            }
            else
            {
              cov_xy = w_y*cov_xuyp + (1.0-w_y)*cov_xuyu;
            }
          }
          else if( w_x == 1.0 ) // w_x is pair-wise
          {
            if( w_y == 0.0 ) // w_y is uniform
            {
              cov_xy = cov_xpyu;
            }
            else if( w_y == 1.0 ) // w_y is pair-wise
            {
              cov_xy = cov_xpyp;
            }
            else
            {
              cov_xy = w_y*cov_xpyp + (1.0-w_y)*cov_xpyu;
            }
          }
          else
          {
            if( w_y == 0.0 ) // w_y is uniform
            {
              cov_xy = w_x*cov_xpyu + (1.0-w_x)*cov_xuyu;
            }
            else if( w_y == 1.0 ) // w_y is pair-wise
            {
              cov_xy = w_x*cov_xpyp + (1.0-w_x)*cov_xuyp;
            }
            else
            {
              cov_xy =   w_x       * w_y       * cov_xpyp
                       + w_x       * (1.0-w_y) * cov_xpyu
                       + (1.0-w_x) * w_y       * cov_xuyp
                       + (1.0-w_x) * (1.0-w_y) * cov_xuyu;
            }
          }
#if 0
  cout << "w_x = " << w_x << ", v_xp = " << v_xp << ", v_xu = " << v_xu
       << ", w_y = " << w_y << ", v_yp = " << v_yp << ", v_yu = " << v_yu
       << ", cov_xpyp = " << cov_xpyp << ", cov_xuyu = " << cov_xuyu
       << endl;
#endif
          ac_matrix(r, c) = cov_xy;
        }
      }
    }
  }
#if 0
  cout << result1.std_err(0, 0).standard_error[WEIGHT_COUNT] << "  "
       << result2.std_err(0, 0).standard_error[WEIGHT_COUNT] << endl;
  print_vec_matrix(ac_vec_matrix, cout, "ac_vec_matrix");
  print_matrix(ac_matrix, cout, "ac_matrix");
#endif
}

void
optimal_weight_finder::find_optimal_weights(const pairset_vector&        pairsets,
                                            const pairset_result_vector& results)
{
  my_weight_matrix.resize(0);

  for( size_t r = 0; r < pairsets.size(); ++r )
  {
    weight_matrix_by_pedigree optimal_weights_by_pedigree;

    find_optimal_weight(pairsets[r], results[r], optimal_weights_by_pedigree);

    my_weight_matrix.push_back(optimal_weights_by_pedigree);
  }
}

void
optimal_weight_finder::find_optimal_weight(const pairset_by_pedigree_type&  pairset,
                                           const pairset_result&            result,
                                                 weight_matrix_by_pedigree& weight)
{
  // Construct optimal weight_matrix for each pedigree.
  //
  for( size_t ped = 0; ped < pairset.size(); ++ped )
  {
    const pairset_type& pset = pairset[ped];
#if 0
  cout << "  " << ped << "(" << pset.size() << ")";
#endif
    weight_matrix optimal_weight;

    if( !pset.size() )
    {
      optimal_weight.first.resize (0,0);
      optimal_weight.second.resize(0,0);
      weight.push_back(optimal_weight);

      continue;
    }

    double p = 1.0;
    double u = 1.0 / double( pset.size() );
#if 0
  cout << " : pair = " << p << ", uniform = " << u << endl;
#endif
    if( p == u )
    {
      optimal_weight.first.resize (result.corr.rows(), result.corr.cols(), p);
      optimal_weight.second.resize(result.corr.rows(), result.corr.cols(), 1.);
      weight.push_back(optimal_weight);
      continue;
    }

    optimal_weight.first.resize (result.corr.rows(), result.corr.cols(), std::numeric_limits<double>::quiet_NaN());
    optimal_weight.second.resize(result.corr.rows(), result.corr.cols(), std::numeric_limits<double>::quiet_NaN());

    for( size_t i = 0; i < result.corr.rows(); ++i )
    {
      for( size_t j = 0; j < result.corr.cols(); ++j )
      {
        double V0 = result.std_err(i, j).standard_error[UNIFORM];
        double V2 = result.std_err(i, j).standard_error[PAIR_WISE];
        double V1 = result.std_err(i, j).standard_error[MEAN];

        double w1ped = p;
        double min_v = 1.;

        if( V1 < V0 && V1 < V2 )
        {
          double b1 = -3.0*V0 + 4.0*V1 -     V2;
          double b2 =  2.0*V0 - 4.0*V1 + 2.0*V2;

          double w1 = -b1 / (2.0*b2);

          w1ped = u + w1*(1.0 - u);
          min_v = 0.5;
        }
        else if( V0 < V2 )
        {
          w1ped = u;
          min_v = 0.;
        }

        optimal_weight.first(i, j) = w1ped;
        optimal_weight.second(i,j) = min_v;
      }
    }

    weight.push_back(optimal_weight);
#if 0
    //cout << long_relationship_name(...) << " : " << ped << endl;
    debug_view(optimal_weight);
#endif
  }
}

void
optimal_weight_finder::debug_view(const weight_matrix& m) const
{
  for( size_t t1 = 0; t1 < m.first.rows(); ++t1 )
  {
    for( size_t t2 = 0; t2 < m.first.cols(); ++t2 )
    {
      cout << setprecision(10);
      cout.setf(ios_base::fixed,ios_base::floatfield);
      cout << m.first(t1, t2) << ", " << m.second(t1, t2) << "	";
    }
    cout << endl;
  }
  cout << endl;
}

} // end of namespace FCOR
} // end of namespace SAGE
