#include <iostream>
#include <iomanip>
#include <vector>
#include <functional>
#include "numerics/clapack.h"
#include "numerics/fmatrix.h"

using namespace std;

namespace SAGE {

typedef FortranMatrix<double> matrix;

//
// SVD
//
size_t SVD::compute(const matrix& A)
{
  if( !A || !A.rows() || !A.cols() )
  {
    U.setstate(matrix::failbit);
    S.setstate(matrix::failbit);
    Vt.setstate(matrix::failbit);
    return 0;
  }

  U.clear();
  S.clear();
  Vt.clear();

  matrix AA  = A;
  char jobu  = 'A';
  char jobvt = 'A';
  int m      = AA.rows();
  int n      = AA.cols();
  int lda    = AA.lda();

  U.resize_nofill(m,m);
  Vt.resize_nofill(n,n);

  std::vector<double>  sigma( std::min(m,n) );

  int ldu = U.lda();
  int ldv = Vt.lda();
  int info;

  int minmn = std::min(A.rows(), A.cols());
  int maxmn = std::max(A.rows(), A.cols());
  int lwork = 8*std::max(3*minmn+maxmn,5*minmn);

  if( (int) work.size() < lwork )
    work.resize(lwork);

  F77_CALL(dgesvd)(&jobu, &jobvt, 
          &m, &n, AA.raw_storage(), &lda,
          &sigma[0], 
          U.raw_storage(), &ldu,
          Vt.raw_storage(), &ldv,
          &work[0], &lwork,
          &info);

  int return_code = 1;
  S.fill(0);
  S.resize_fill(m,n,0);
  for(int i = 0; i < std::min(m,n); ++i)
  {
    if( fabs(sigma[i]) < 1.0e-10 )
      S(i,i) = 0.;
    else
    {
      S(i,i) = sigma[i];

      if( fabs(sigma[i]) < 1.0e-7 )
        return_code = 2;
    }
  }

  positive_diagonal.resize(0);

  if( info != 0 )
  {
    U.setstate(matrix::failbit);
    S.setstate(matrix::failbit);
    Vt.setstate(matrix::failbit);

    return 0;
  }

  double tol = std::max(std::sqrt(std::numeric_limits<double>::epsilon()) * S(0,0), 0.0);

  for( size_t i = 0; i < std::min(S.cols(),S.rows()); ++i )
  {
    if( S(i,i) > tol )
      positive_diagonal.push_back(i);
  }

  return return_code;
}

matrix& SVD::inverse_of(const matrix& A, matrix& Ai)
{
  compute(A);
  return inverse(Ai);
}

matrix SVD::inverse_of(const matrix& A)
{
  compute(A);
  matrix temp;
  return inverse(temp);
}

matrix SVD::inverse()
{
  matrix temp;
  return inverse(temp);
}

matrix& SVD::inverse(matrix& Ai)
{
  if(fail())
  {
    Ai.setstate(matrix::failbit);
    return Ai;
  }

  Ai.clear();

  temp = U;

  for(matrix::size_type j = 0; j < std::min(S.cols(),S.rows()); ++j) 
  {
    // multiply with diagonal 1/W
    double weight = S(j, j);
    // invert only non singular weight, and multiply row
    if(weight != 0)
      for(matrix::size_type i = 0; i < temp.rows(); ++i)
        temp(i, j) /= weight;
    else
      // for singular weight, replace 1/0 with 0, so zero out row.
      for(matrix::size_type i = 0; i < temp.rows(); ++i)
        temp(i, j) = 0;
  }

  multiply(temp, Vt, Ai);
  Ai.transpose();  

  return Ai;
}

matrix& SVD::sqrt(matrix& R)
{
  if(fail())
  {
    R.setstate(matrix::failbit);
    return R;
  }

  R.clear();

  temp = U;

  for(matrix::size_type j = 0; j < std::min(S.cols(),S.rows()); ++j) 
  {
    // multiply with diagonal W
    double weight = ::sqrt(S(j, j));
    // invert only non singular weight, and multiply row
    for(matrix::size_type i = 0; i < temp.rows(); ++i)
      temp(i, j) *= weight;
  }

  multiply(temp, Vt, R);
  R.transpose();  

  return R;
}

matrix& SVD::inverse_sqrt(matrix& Ri)
{
  if(fail())
  {
    Ri.setstate(matrix::failbit);
    return Ri;
  }

  Ri.clear();

  temp = U;

  for(matrix::size_type j = 0; j < std::min(S.cols(),S.rows()); ++j) 
  {
    // multiply with diagonal 1/W^.5
    double weight = ::sqrt(S(j, j));
    // invert only non singular weight, and multiply row
    if(weight != 0)
      for(matrix::size_type i = 0; i < temp.rows(); ++i)
        temp(i, j) /= weight;
    else
      // for singular weight, replace 1/0 with 0, so zero out row.
      for(matrix::size_type i = 0; i < temp.rows(); ++i)
        temp(i, j) = 0;
  }

  multiply(temp, Vt, Ri);
  Ri.transpose();  

  return Ri;
}

// general inverse
matrix& SVD::general_inverse_of(const matrix& A, matrix& Ai)
{
  compute(A);
  return general_inverse(Ai);
}

matrix SVD::general_inverse_of(const matrix& A)
{
  compute(A);
  matrix temp;
  return general_inverse(temp);
}

matrix SVD::general_inverse()
{
  matrix temp;
  return general_inverse(temp);
}

matrix& SVD::general_inverse(matrix& Ai)
{
  Ai.clear();

  matrix V = Vt;
  V.transpose();

  matrix new_S(positive_diagonal.size(), positive_diagonal.size());
  matrix new_U(U.rows(), positive_diagonal.size());
  matrix new_V(V.rows(), positive_diagonal.size());

  for( size_t i = 0; i < S.rows(); ++i )
  {
    for( size_t j = 0; j < positive_diagonal.size(); ++j )
    {
      new_S(j,j) = S(positive_diagonal[j], positive_diagonal[j]);
      new_V(i,j) = V(i,positive_diagonal[j]);
      new_U(i,j) = U(i,positive_diagonal[j]);
    }
  }

  matrix new_Ut = new_U;
  new_Ut.transpose();

  for( size_t i = 0; i < new_S.rows(); ++i )
  {
    // multiply with diagonal 1/W
    double weight = new_S(i,i);
    for( size_t j = 0; j < new_Ut.cols(); ++j )
      new_Ut(i,j) /= weight;
  }

  multiply(new_V, new_Ut, Ai);

  return Ai;
}

matrix& SVD::general_sqrt(matrix& R)
{
  R.clear();

  matrix V = Vt;
  V.transpose();

  matrix new_S(positive_diagonal.size(), positive_diagonal.size());
  matrix new_U(U.rows(), positive_diagonal.size());
  matrix new_V(V.rows(), positive_diagonal.size());

  for( size_t i = 0; i < S.rows(); ++i )
  {
    for( size_t j = 0; j < positive_diagonal.size(); ++j )
    {
      new_S(j,j) = S(positive_diagonal[j], positive_diagonal[j]);
      new_V(i,j) = V(i,positive_diagonal[j]);
      new_U(i,j) = U(i,positive_diagonal[j]);
    }
  }

  matrix new_Ut = new_U;
  new_Ut.transpose();

  for( size_t i = 0; i < new_S.rows(); ++i )
  {
    // multiply with diagonal sqrt(W)
    double weight = ::sqrt(new_S(i,i));
    for( size_t j = 0; j < new_Ut.cols(); ++j )
      new_Ut(i,j) *= weight;
  }

  multiply(new_V, new_Ut, R); 

  return R;
}

matrix& SVD::general_inverse_sqrt(matrix& Ri)
{
  Ri.clear();

  matrix V = Vt;
  V.transpose();

  matrix new_S(positive_diagonal.size(), positive_diagonal.size());
  matrix new_U(U.rows(), positive_diagonal.size());
  matrix new_V(V.rows(), positive_diagonal.size());

  for( size_t i = 0; i < S.rows(); ++i )
  {
    for( size_t j = 0; j < positive_diagonal.size(); ++j )
    {
      new_S(j,j) = S(positive_diagonal[j], positive_diagonal[j]);
      new_V(i,j) = V(i,positive_diagonal[j]);
      new_U(i,j) = U(i,positive_diagonal[j]);
    }
  }

  matrix new_Ut = new_U;
  new_Ut.transpose();

  for( size_t i = 0; i < new_S.rows(); ++i )
  {
    // multiply with diagonal 1/W^.5
    double weight = ::sqrt(new_S(i,i));
    for( size_t j = 0; j < new_Ut.cols(); ++j )
      new_Ut(i,j) /= weight;
  }

  multiply(new_V, new_Ut, Ri); 

  return Ri;
}

double SVD::determinant() const
{
  if(fail() || U.rows() != Vt.cols() || !S.cols() )
    return std::numeric_limits<double>::quiet_NaN();

  double product = 1;
  for(size_t k = 0; k < S.cols(); ++k)     
    product *= S(k, k);                    
  return product;                          
}

double SVD::well_condition() const
{
  if(fail())
    return std::numeric_limits<double>::quiet_NaN();

  double min_w = fabs(S(0,0));
  double max_w = min_w;
  for(matrix::size_type k = 1; k < S.cols(); ++k) 
  {
    // find max and min absolute diagonals in W.
    double weight = fabs(S(k, k));
    min_w = std::min(min_w, weight);
    max_w = std::max(max_w, weight);
  }
  return min_w / max_w;
}

matrix& SVDinverse(const matrix& A, matrix& Ai)
{
  SVD svd(A);
  return svd.inverse(Ai);
}

matrix SVDinverse(const matrix& A)
{
  SVD svd(A);
  return svd.inverse();
}

//
// GQR
//
void GQR::compute(const matrix& A, const matrix& B, bool reduced)
{
  if(!A || !A.rows() || !A.cols() || !B || !B.cols() || B.rows() != A.rows())
  {
    Q.setstate(matrix::failbit);
    R.setstate(matrix::failbit);
    T.setstate(matrix::failbit);
    Z.setstate(matrix::failbit);
    return;
  }

  Q.clear();
  R.clear();
  Z.clear();

  matrix AA  = A;
  matrix BB  = B;

  int n = AA.rows();
  int m = AA.cols();
  int p = BB.cols();
  int lda = AA.lda();

  std::vector<double> tau_A( std::min(m,n) );
  std::vector<double> tau_B( std::min(n,p) );

  int ldb = BB.lda();

  int minmn   = std::min(m,n);
  int maxmn   = std::max(m,n);
  int maxmnp = std::max(maxmn,p);
  int nb1    = n*8;
  int lwork = maxmnp*nb1;
  int info;

  if((int) work.size() < lwork)
    work.resize(lwork);

  // Compute GQR
  F77_CALL(dggqrf)(&n, &m, &p,
                   AA.raw_storage(), &lda,
                   &tau_A[0], 
                   BB.raw_storage(), &ldb,
                   &tau_B[0], 
                   &work[0], &lwork, 
                   &info);

  if(info != 0)
  {
    Q.setstate(matrix::failbit);
    R.setstate(matrix::failbit);
    T.setstate(matrix::failbit);
    Z.setstate(matrix::failbit);
    return;
  }

  // Form Q
  Q = AA;

  int dim = m;
  if(!reduced)
  {
    Q.resize_nofill(n,n,0);
    dim = n;
  }
  int ldq = Q.lda();

  F77_CALL(dorgqr)(&n, &dim, &minmn,
                   Q.raw_storage(), &ldq,
                   &tau_A[0], 
                   &work[0], &lwork, 
                   &info);

  if(info != 0)
  {
    Q.setstate(matrix::failbit);
    R.setstate(matrix::failbit);
    T.setstate(matrix::failbit);
    Z.setstate(matrix::failbit);
    return;
  }

  // Form Z
  Z = BB;
  int ldz = Z.lda();

  F77_CALL(dorgrq)(&p, &n, &p,
                    Z.raw_storage(), &ldz,
                    &tau_B[0], 
                    &work[0], &lwork, 
                    &info);

  if(reduced)
    Z.resize_nofill(m,n);

  if(info != 0)
  {
    Q.setstate(matrix::failbit);
    R.setstate(matrix::failbit);
    T.setstate(matrix::failbit);
    Z.setstate(matrix::failbit);
    return;
  }

  // Form R
  dim = m;
  if(!reduced)
    dim = n;

  R.resize_fill(dim,m,0);
  for(int i = 0; i < m; ++i)
    for(int j = i; j < m; ++j)
      R(i,j) = AA(i,j);

  // Form T
  if(reduced)
  {
    T.resize_fill(m,m,0);
    for(int i = 0; i < m; ++i)
      for(int j = i+p-n; j < m; ++j)
        T(i,j) = BB(i,j);
  }
  else
  {
    T.resize_fill(n,p,0);
    for(int i = 0; i < n; ++i)
      for(int j = i+p-n; j < p; ++j)
        T(i,j) = BB(i,j);
  }
}

//
// LUP
//
void LUP::compute(const matrix& A)
{
  if(!A || !A.rows() || !A.cols())
  {
    LU.setstate(matrix::failbit);
    return;
  }

  LU.clear();
  LU = A;

  int rows = LU.rows();
  int cols = LU.cols();
  int lda  = LU.lda();
  int info;

  pivots.resize( std::min(rows,cols) );
  F77_CALL(dgetrf)(&rows, &cols, LU.raw_storage(), &lda, &pivots[0], &info);

  if(info != 0)
    LU.setstate(matrix::failbit);
}

matrix& LUP::inverse(matrix& Ai)
{
  if(fail())
  {
    Ai.setstate(matrix::failbit);
    return Ai;
  }

  Ai.clear();
  Ai = LU;
  int cols = Ai.cols();
  int lda  = Ai.lda();
  int info;

#if 0
  int lwork = -1;

  if(work.size() < 1)
    work.resize(1);

  F77_CALL(dgetri)(&cols, Ai.raw_storage(), &lda, &pivots[0], &work[0], &lwork, &info);

  lwork = cols*work[0];
#else
  int lwork = 128*Ai.cols();
#endif

  if( lwork > (int) work.size())
    work.resize(lwork);

  F77_CALL(dgetri)(&cols, Ai.raw_storage(), &lda, &pivots[0], &work[0], &lwork, &info);

  if(info != 0)
    Ai.setstate(matrix::failbit);

  return Ai;
}

matrix& LUP::inverse_of(const matrix& A, matrix& Ai)
{
  compute(A);
  return inverse(Ai);
}

matrix LUP::inverse_of(const matrix& A)
{
  compute(A);
  matrix temp;
  return inverse(temp);
}

matrix LUP::inverse()
{
  matrix temp;
  return inverse(temp);
}

double LUP::determinant() const
{
  if(fail() || LU.rows() != LU.cols() )
    return std::numeric_limits<double>::quiet_NaN();

  double product = 1;
  for(size_t k = 0; k < LU.cols(); ++k)     
    product *= LU(k, k);                    
  return product;                          
}

matrix& LUinverse(const matrix& A, matrix& Ai)
{
  LUP lup(A);
  return lup.inverse(Ai);
}

matrix LUinverse(const matrix& A)
{
  LUP lup(A);
  return lup.inverse();
}

matrix& Cholesky(const matrix& A, matrix& UL, bool upper)
{
  if(!A || !A.rows() || A.rows() != A.cols() )
  {
    UL.setstate(matrix::failbit);
    return UL;
  }

  UL.clear();

  UL = A;
  char upl0 = upper? 'U' : 'L';
  int n     = UL.rows();
  int lda   = UL.lda();
  int info;

  F77_CALL(dpotrf)(&upl0, &n, UL.raw_storage(), &lda, &info);

  if(info != 0)
  {
    UL.setstate(matrix::failbit);
    return UL;
  }

  if(upper)
  {
    for(size_t j = 0; j < UL.cols(); ++j)
      for(size_t i = j + 1; i < UL.rows(); ++i)
        UL(i,j) = 0.0;
  }
  else
  {
    for(size_t j = 0; j < UL.cols(); ++j)
      for(size_t i = 0; i < j; ++i)
        UL(i,j) = 0.0;
  }

  return UL;
}

matrix& TRIinverse(const matrix& UL, matrix& ULi, bool upper)
{
  ULi = UL;
  char upl0 = upper? 'U' : 'L';
  char diag = 'N';
  int n     = UL.rows();
  int lda   = UL.lda();
  int info;

  ULi.clear();
  
  F77_CALL(dtrtri)(&upl0, &diag, &n, ULi.raw_storage(), &lda, &info);

  if(info != 0)
    ULi.setstate(matrix::failbit);
  return ULi;
}

matrix TRIinverse(const matrix& UL, bool upper)
{
  matrix temp;
  return TRIinverse(UL, temp, upper);
}

matrix& SYMinverse(const matrix& S, matrix& Si)
{
  if(!S || !S.rows() || S.rows() != S.cols() )
  {
    Si.setstate(matrix::failbit);
    return Si;
  }

  Si.clear();

  Si = S;
  char upl0 = 'U';
  int n     = Si.rows();
  int lda   = Si.lda();
  int info;

  F77_CALL(dpotrf)(&upl0, &n, Si.raw_storage(), &lda, &info);

  if(info != 0)
  {
    Si.setstate(matrix::failbit);
    return Si;
  }

  F77_CALL(dpotri)(&upl0, &n, Si.raw_storage(), &lda, &info);

  if(info != 0)
  {
    Si.setstate(matrix::failbit);
    return Si;
  }

  for(size_t j = 0; j < Si.cols(); ++j)
    for(size_t i = 0; i < j; ++i)
      Si(j,i) = Si(i,j);


  return Si;
}

matrix SYMinverse(const matrix& S)
{
  matrix temp;
  return SYMinverse(S, temp);
}

matrix& Msqrt(const matrix& XX, matrix& X, SVD& svd)
{
  if(!XX || !XX.rows() || XX.rows() != XX.cols())
  {
    X.setstate(matrix::failbit);
    return X;
  }

  X.clear();

  Cholesky(XX, X);

  if( X.good() )
  {
    svd.compute(X);

    if(!svd.fail())
    {
      multiply(svd.S, svd.Vt, svd.temp);
      XTZ(svd.Vt, svd.temp, X);
      return X;
    }
  }    
  X.setstate(matrix::failbit);
  return X;
}

matrix Msqrt(const matrix& XX)
{
  matrix X;
  SVD svd;
  return Msqrt(XX, X, svd);
}

matrix& Minverse_sqrt(const matrix& XX, matrix& X, SVD& svd)
{
  if(!XX || !XX.rows() || XX.rows() != XX.cols())
  {
    X.setstate(matrix::failbit);
    return X;
  }

  X.clear();

  Cholesky(XX, X, false);
  TRIinverse(X, X, false);
  
  if( X.good() )
  {
    svd.compute(X);

    if(!svd.fail())
    {
      multiply(svd.S, svd.Vt, svd.temp);
      XTZ(svd.Vt, svd.temp, X);
      return X;
    }
  }    
  X.setstate(matrix::failbit);
  return X;
}

matrix Minverse_sqrt(const matrix& XX)
{
  matrix X;
  SVD svd;
  return Minverse_sqrt(XX, X, svd);
}

matrix& forward_substitute(matrix& L, matrix &B)
{
  if( L.rows() != B.rows() || L.rows() != L.cols() || B.cols() != 1)
  {
    B.setstate(matrix::failbit);
    return B;
  }

  B.clear();
  size_t n = L.cols();

  B(0,0) /= L(0,0);
  for(size_t i = 1; i < n; ++i)
  {
    double sum = 0.0;
    for(size_t j = 0; j < i; ++j)
      sum += L(i,j)*B(j,0);
    B(i,0) = (B(i,0)-sum)/L(i,i);
  }
  return B;
}

matrix& back_substitute(matrix& U, matrix &B)
{
  if( U.rows() != B.rows() || U.rows() != U.cols() || B.cols() != 1)
  {
    B.setstate(matrix::failbit);
    return B;
  }

  B.clear();
  size_t n = U.cols();

  B(n-1,n-1) /= U(n-1,n-1);
  for(long i = n-1; i >= 0; --i)
  {
    double sum = 0.0;
    for(size_t j = i+1; j < n; ++j)
      sum += U(i,j)*B(j,0);
    B(i,0) = (B(i,0)-sum)/U(i,i);
  }
  return B;  
}

int Eigen(matrix& b, vector<double>& w)
{
  vector<double> work;

  char jobz = 'V';
  char uplo = 'U';
  int n     = b.rows();
  int lda   = b.lda();
  int lwork = 8*(3*n-1);
  int info;

  w.resize(n);
  work.resize(lwork);

  F77_CALL(dsyev)(&jobz, &uplo, &n, b.raw_storage(), &lda, &w[0], &work[0], &lwork, &info);
        
  return info;
}

} // end of namespace SAGE

#ifdef NEEDS_FORTRAN_MAIN

// Dummy function to keep F2C happy
extern "C" int MAIN_()
{
  return 1;
}

// Dummy function to keep F2C happy
extern "C" int MAIN__()
{
  return 1;
}

#endif
