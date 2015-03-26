#ifndef GLS3_H
#define GLS3_H

//=============================================================================
//
// Author:  Kevin Jacobs
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structure.
//
//            class UnivariateGeneralizedLeastSquares3
//
// Copyright (c) 2001   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/sib_matrix.h"

namespace SAGE   {
namespace SIBPAL {

class UnivariateGeneralizedLeastSquares3
{
  public:

    enum inverse_type       { SVD = 0, LUP };
      
    UnivariateGeneralizedLeastSquares3(size_t m = 0);

    ~UnivariateGeneralizedLeastSquares3();

    inverse_type      inverse_method()                     const;
    void          set_inverse_method(inverse_type d);
    void          set_inverse_method(const string& l);

    bool              leverage_adjustment() const;
    void          set_leverage_adjustment(bool adjust);

    void reset_results();
    void reset();
    void reset(size_t m, size_t n = 1);

    void add_block(matrix& y, const matrix& A);
    void add_block(matrix& y, const matrix& W, const matrix& A, const weight_status_type& status = INVERSEW);
    void add_block_kron(matrix& y, const matrix& W, const matrix& A, const weight_status_type& status = INVERSEW);

    void compute_covariance(const matrix& X, matrix& C);
    bool compute();

    void build_residuals(const matrix& y, const matrix& A, matrix& r) const;

    // For independence working model
    void update_robust_variance_block(matrix& y, const matrix& A);
    
    // For specified covariance matrix
    void update_robust_variance_block(matrix& y, const matrix& W, const matrix& A, const weight_status_type& status = INVERSEW);
    bool compute_robust_variance();

    size_t get_svd_return_code() const;

    size_t                  parameters;
    size_t                  variates;
    size_t                  observation_count;
    size_t                  cluster_count;
    size_t                  svd_return_code;
    matrix                  beta;
    matrix                  hat;
    matrix                  information;
    matrix                  Variance;
    matrix                  SandwichVariance;
    matrix                  residual_variance;
    matrix                  total_variance;

    // Decompositions
    SAGE::SVD svd;
    SAGE::LUP lup;
    
  private:

    inverse_type            my_inverse_method;
    bool                    my_leverage_adjustment;
      
    // Temps
    matrix          B;             // partial W^-.5
    matrix          BA;            // partial B A
    matrix          By;            // partial B y
    stable_matrix   yWy;           // weighted sum of squares
    stable_matrix   AWA;           // sum(A' W A)
    stable_matrix   AWy;           // sum(A' W y)
    stable_matrix   OPG;           // Outer product gradient: sum(A' W E(y) E(y') W A)
    matrix temp1,temp2,temp3;
};

typedef UnivariateGeneralizedLeastSquares3 GLS3;

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif
