#include <vector>
#include "numerics/svd.h"
#include "numerics/kahan.h"

namespace SAGE {

double cnorm(const Matrix2D<double>& m, size_t col, size_t row)
{
  KahanAdder<double> ssum = 0.0;
  
  for(size_t i = row; i < m.rows(); ++i)
    ssum += m(i,col)*m(i,col);
    
  return sqrt(ssum);
}

void cholesky(const Matrix2D<double>& A, Matrix2D<double>& G)
{
  assert(A.rows() == A.cols());

  G = A;

  size_t m = A.cols();

  for(size_t k = 0; k < m; ++k)
  {
    for(size_t j = k + 1; j < m; ++j)
      for(size_t i = j; i < m; ++i)
        G(j,i) -= G(k,i)*G(k,j)/G(k,k);

    for(size_t i = 0; i < k; ++i)
      G(k,i) = 0.0;

    double rkk = sqrt(G(k,k));
    for(size_t i = k; i < m; ++i)
      G(k,i) /= rkk;
  }
}

Matrix2D<double>& msqrt(const Matrix2D<double>& XX, Matrix2D<double>& X)
{
  X.clear();

  assert(XX.rows() == XX.cols());
  Matrix2D<double> U,S,V;
  cholesky(XX, X);
  Transpose(X);
  svd(X, U, S, V);

  if(!U || !S || !V)
  {
    X.setstate(Matrix2D<double>::failbit);
    return X;
  }

  Multiply(U, S, V);
  Transpose(U);
  Multiply(V, U, X);
  return X;
}

Matrix2D<double> msqrt(const Matrix2D<double>& XX)
{
  assert(XX.rows() == XX.cols());
  Matrix2D<double> X;
  return msqrt(XX, X);
}

void qr(const Matrix2D<double>& A, Matrix2D<double>& Q, Matrix2D<double>& R)
{
    Q.clear();
    R.clear();

    size_t n = A.rows();
    size_t m = A.cols();
    size_t minmn = std::min(m,n);

    Q.fill(0);
    Q.resize(n,m, 0.0);
    R.fill(0);
    R.resize(minmn,minmn, 0.0);
    Q=A;

    Matrix1D<double> sigma(minmn);

    /* REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS. */
    for (size_t j = 0; j < minmn; ++j) 
    {
        /* COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE    */
        /* J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR. */

	double ajnorm = cnorm(Q, j, j);

        sigma[j] = 0.0;

	if (ajnorm == 0.0)
            continue;

	if (Q(j, j) < 0.0) 
	    ajnorm = -ajnorm;

        sigma[j] = -ajnorm;

	for (size_t i = j; i < m; ++i) 
	    Q(i, j) /= ajnorm;

	Q(j, j) += 1;

        /* APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS */
        /* AND UPDATE THE NORMS. */

	if (n < j)
          continue;

	for (size_t k = j+1; k < n; ++k) 
	{
	    KahanAdder<double> sum = 0.0;
	    for (size_t i = j; i < m; ++i)
		sum += Q(i, j) * Q(i, k);

	    double temp = sum / Q(j, j);

	    for (size_t i = j; i < m; ++i)
		Q(i, k) -= temp * Q(i, j);
	}
    }

    R = Q.extract(minmn, minmn);
    for(size_t j = 0; j < minmn; ++j)
      R(j,j) = sigma[j];

    Matrix1D<double> work(m);
    for (size_t j = 1; j < minmn; ++j)
      for (size_t i = 0; i < j; ++i)
      {
        Q(i, j) = 0.0;                 
        R(j, i) = 0.0;
      }

    if( m > n+1)
      for (size_t j = n; j < m; ++j)    
      {
        for (size_t i = 0; i < m; ++i) 
            Q(i, j) = 0.0;
        Q(j, j) = 1.0;    
      }

    for (long k = minmn-1; k >= 0; --k) 
    {
	for (size_t i = k; i < m; ++i) 
	{
	    work[i] = Q(i, k);
	    Q(i, k) = 0.0;
	}

	Q(k, k) = 1.0;

	if (work[k] == 0.0)
	    continue;

	for (size_t j = k; j < m; ++j) 
	{
	    KahanAdder<double> sum = 0.0;
	    for (size_t i = k; i < m; ++i)
		sum += Q(i, j) * work[i];

	    double temp = sum / work[k];
	    for (size_t i = k; i < m; ++i)
		Q(i, j) = Q(i, j) - temp * work[i];
	}
    }
}             

}
