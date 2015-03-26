#ifndef SVD_H
#define SVD_H

#include "numerics/matrix.h"
#include "numerics/matfunctions.h"

// svd(u, w, v) -- Computes the singular value decomposition of THIS matrix,
// and returns the 3 matrices U, W, V, such that:
//      A = U * W * V.transpose().      O(n^3).
//
// Small diagonals in W should be zeroed out, to clearly mark singularities
// in the matrix. The range space are columns of U such that the
// corresponding diagonals in W are non zero. The null space are columns of
// V such that the corresponding diagonals in W are zero. The vectors in U
// and V are orthonormal. If the matrix is symmetric, then U and V are equal
// and they contain the eigenvectors, W contains the eigenvalues.

// The eigenvalues are positive if the matrix is positive definite. The
// matrices U, W, V are returned for the user to mark out singularities
// based on some lower bound on the well-condition number, such as 1.0e-6.
// These matrices are cached so that subsequent computations such as rank,
// determinant, condition number, inverse, back substitution, will not need
// to compute the expensive SVD again.
// See Numerical Recipees for reference.

template<class T>
void singular_value_decomposition(const Matrix2D<T>& A,
                                        Matrix2D<T>& u, 
                                        Matrix2D<T>& w, 
                                        Matrix2D<T>& v)
{
    // dimension of vectors in U,
    // dimension of vectors in V.
    // fill *this into u.
    size_t m = std::max(A.rows(), A.cols()); 
    size_t n = A.cols();                    
    u.resize(m, n);                             
    u.fill(0);
    u.update(A);

    // clear the extra rows
    for(size_t i = A.rows(); i < m; ++i)   
      for(size_t j = 0; j < n; ++j)
        u(i,j) = 0;

    w.resize(n, n);
    w.fill(0);
    v.resize(n, n);

    // do all the work and store results in place
    svdcmp(u, m, n, w, v);
}

template<class T>
void svd(const Matrix2D<T>& A,
               Matrix2D<T>& u, 
               Matrix2D<T>& w, 
               Matrix2D<T>& v)
{
  singular_value_decomposition(A,u,w,v);
}

// singularities -- Returns the number of singular weights which make
// the well-condition number less than given lower bound threshold. O(n).
// Default lower bound on the well-condition number is 1.0e-6.
// These singular weights are marked by setting them to zero.
// The inverse of these zero weights will be set to zero instead of
// infinite, which corresponds to finding least-norm solution.

template<class T>
int singularities(Matrix2D<T>& w, 
                  float lbound)
{
  T max_w = 0;

  // find max absolute weight in diagonal
  for(size_t k = 0; k < w.cols(); ++k) 
  {
    T weight = fabs(w(k, k));
    if(max_w < weight) max_w = weight;
  }

  T lbound_w = lbound * max_w;

  // count the number of singular diagonals
  size_t count = 0;

  // find near-singular weights and zero them out
  for(size_t k = 0; k < w.cols(); ++k) 
  {
    T weight = fabs(w(k, k));
    if(weight < lbound_w) 
    {
      w(k, k) = 0;
      ++count;                 
    }
  }
  return count;
}


// rank(u, w, v) -- Returns the rank of a rectangular matrix,
// which is the number of nonzero diagonal weights in W. O(n).
// Elimination is done on the matrix and singularities are marked away
// if W is empty, otherwise the decomposition cached in (U, W, V) is used.

template<class T>
int rank(Matrix2D<T>& w)
{
  int count = 0;
  for(size_t k = 0; k < w.cols(); ++k)
    if(w(k, k) != 0) 
      ++count;
  return count;
}

// 1/condition_number(u, w, v) -- Returns the inverse of the condition number
// of a rectangular matrix, or degree of well-condition in range [0,1].
// The well-condition number is the ratio of the smallest weight to the
// largest weight in W. O(n). Ill-conditioned matrices with very small
// well-condition number should not be inverted without marking off the
// singularities.
// Elimination is done on the matrix if W is empty, otherwise the
// decomposition cached in (U, W, V) is used.

template<class T>
T well_condition(Matrix2D<T>& w)
{
  T min_w = fabs(w(0,0));
  T max_w = min_w;
  for(size_t k = 1; k < w.cols(); ++k) 
  {
    // find max and min absolute diagonals in W.
    T weight = fabs(w(k, k));
    if(min_w > weight)
      min_w = weight;
    if(max_w < weight)
      max_w = weight;
  }
  return min_w / max_w;
}

// determinant(u, w, v) -- Returns the determinant of a square matrix,
// which is the product of the diagonal weights in W, O(n).
// Determinant of a rectangular matrix is not well defined, 0 is returned.
// Elimination is done on the matrix and singularies are marked away
// if W is empty, otherwise the decomposition cached in (U, W, V) is used.

template<class T>
T determinant(Matrix2D<T>& u, Matrix2D<T>& w, Matrix2D<T>& v)
{
  if(u.rows() != u.cols()) {
//    cerr << "Return 0 for determinant of a rectangular matrix." << endl;
    return 0;
  }
  T product = 1;
  for(size_t k = 0; k < w.cols(); ++k)                // product of diagonals in W.
    product *= w(k, k);                         // small weights should be zeroed
  return product;                               // out, so that determinant = 0.
}

// -- Returns the inverse of a square matrix,
// which is: V * (1/W) * U.transpose(), O(n^3).
// Small weights in W should be zeroed out to clearly mark singularities
// in the matrix. If the matrix is rectangular or  singular, the
// pseudo-inverse A^(-1) is returned, x = A^(-1) * y, which
// corresponds to finding least-norm and least-square error x,
// such that: A * x = y.
// Elimination is done on the matrix and singularies are marked away
// if W is empty, otherwise the decomposition cached in (U, W, V) is used.

template<class T>
Matrix2D<T> inverse(Matrix2D<T>& u, Matrix2D<T>& w, Matrix2D<T>& v)
{
  Matrix2D<T> result = v;

  size_t i, j;
  for(i = 0; i < result.rows(); ++i) 
  {
    // multiply with diagonal 1/W
    T weight = w(i, i);
    // invert only non singular weight, and multiply row
    if(weight != 0)
      for(j = 0; j < result.cols(); ++j)
        result(i, j) /= weight;
    else
      // for singular weight, replace 1/0 with 0, so zero out row.
      for(j = 0; j < result.cols(); ++j)
        result(i, j) = 0;
  }

  Matrix2D<double> temp;
  Transpose(result);
  Transpose(u);
  Multiply(result, u, temp);
  result = temp;

  // throw away extra zero block on the right
#if NOT_YET
  if(this->rows() < this->cols())          
    result = result.extract(this->cols(),this->rows());
#endif

  return result;
}

#if 0
// -- Solves for the least-norm and least-square-error solution x, such that
// A x = y, O(n^2). The vectors x are represented as column vectors.
// Elimination is done on the matrix and singularies are marked away
// if W is empty, otherwise the decomposition cached in (U, W, V) is used.

template<class T>
Matrix2D<T> back_substitution(const Matrix2D<T>& u, 
                              const Matrix2D<T>& w, 
                              const Matrix2D<T>& v, 
                              const Matrix2D<T>& y)
{
  // solution matrix
  Matrix2D<T> x;

  if(this->rows() < this->cols()) 
  { 
    // augment y with extra rows of zeros, so that it match cols of u.transpose.
    Matrix2D<T> yy(u.rows(), y.cols(), 0);
    yy.update(y);
    x = u.transpose() * yy;
  } 
  else
    x = u.transpose() * y;

  size_t i, j;
  for(i = 0; i < x.rows(); ++i) 
  {    
    // multiply with diagonal 1/W
    T weight = w(i, i);
    if(weight != 0)                  // invert only non singular
      for(j = 0; j < x.cols(); ++j)  // weight, and multiply row
        x(i, j) /= weight;
    else
      for(j = 0; j < x.cols(); ++j)  // for singular weight, replace
        x(i, j) = 0;                  // 1/0 with 0, so zero out row.
  }
  x = Multiply(v,x);                  // premultiply with v.
  return x;
}
#endif

template<class A, class B>
A sign(const A& a, const B& b)
{
  return (b>=0)? a : -a;
}

template<class A>
A pythag(const A& a, const A& b)
{
  if( fabs(a) > fabs(b) )
  {
    A c = fabs(b/a);
    return fabs(a)*sqrt(1+c*c);
  }
  else if(b != 0)
  {
    A c = fabs(a/b);
    return fabs(b)*sqrt(1+c*c);
  }
  return 0;
}

// svdcmp 
// adapted from Numerical Recipees.
// Indices ranges in [0,n-1] instead of [1,n].

template<class T>
void svdcmp(Matrix2D<T>& a, int m, int n, Matrix2D<T>& w, Matrix2D<T>& v)
{
  int flag,i,its,j,jj,k,l=0,nm=0;
  double c,f,h,s,x,y,z;
  double anorm=0, g=0, scale=0;
  std::vector<T> rv1(n,0);

  // 1. Householder reduction to bidiagonal form
  for (i=0;i<n;i++) {
    l=i+1;
    rv1[i]=(scale*g);
    g=s=scale=0;
    if (i < m) {
      for (k=i;k<m;k++) scale += fabs(a(k,i));
      if (scale!=0.0) {
	for (k=i;k<m;k++) {
	  a(k,i) /= scale;
	  s += a(k,i)*a(k,i);
	}
	f=a(i,i);
	g = -sign(sqrt(s),f);
	h=f*g-s;
	a(i,i)=(f-g);
	if (i != n) {
	  for (j=l;j<n;j++) {
	    for (s=0,k=i;k<m;k++) s += a(k,i)*a(k,j);
	    f=s/h;
	    for (k=i;k<m;k++) a(k,j) += (f*a(k,i));
	  }
	}
	for (k=i;k<m;k++) a(k,i) = (a(k,i) * scale);
      }
    }
    w(i,i)=(scale*g);
    g=s=scale=0;
    if (i < m && i != n) {
      for (k=l;k<n;k++) scale += fabs(a(i,k));
      if (scale!=0.0) {
	for (k=l;k<n;k++) {
	  a(i,k) =  (a(i,k) / scale);
	  s += a(i,k)*a(i,k);
	}
	f=a(i,l);
	g =  -sign(sqrt(s),f);
	h=f*g-s;
	a(i,l)=(f-g);
	for (k=l;k<n;k++) rv1[k]=(a(i,k)/h);
	if (i != m) {
	  for (j=l;j<m;j++) {
	    for (s=0,k=l;k<n;k++) s += a(j,k)*a(i,k);
	    for (k=l;k<n;k++) a(j,k) += (s*rv1[k]);
	  }
	}
	for (k=l;k<n;k++) a(i,k) = (a(i,k) * scale);
      }
    }
    anorm= std::max(anorm,(fabs(w(i,i))+fabs(rv1[i])));
  }

  // 2. Accumulation of right-hand transform V.
  for (i=n-1;i>=0;i--) {
    if (i < n) {
      if (g!=0.0) {
	for (j=l;j<n;j++)
	  v(i,j)= ((a(i,j)/a(i,l))/g);		// double div to avoid underflow
	for (j=l;j<n;j++) {
	  for (s=0,k=l;k<n;k++) s += a(i,k)*v(j,k);
	  for (k=l;k<n;k++) v(j,k) += (s*v(i,k));
	}
      }
      for (j=l;j<n;j++) v(j,i)=v(i,j)=0;
    }
    v(i,i)=1;
    g=rv1[i];
    l=i;
  }

  // 3. Accumulation of left-hand transform U.
  for (i=n-1;i>=0;i--) {
    l=i+1;
    g=w(i,i);
    if (i < n)
      for (j=l;j<n;j++) a(i,j)=0;
    if (g!=0.0) {
      g=1/g;
      if (i != n) {
	for (j=l;j<n;j++) {
	  for (s=0,k=l;k<m;k++) s += a(k,i)*a(k,j);
	  f=(s/a(i,i))*g;
	  for (k=i;k<m;k++) a(k,j) += (f*a(k,i));
	}
      }
      for (j=i;j<m;j++) a(j,i) *= g;
    } else {
      for (j=i;j<m;j++) a(j,i)=0;
    }
    a(i,i) += 1.0;
  }

  // 4. Diagonalization of bidiagonal form
  for (k=n-1;k>=0;k--) { // loop over singular values
    for (its=1;its<=30;its++) {	// loop over allowed iterations
      flag=1;
      for (l=k;l>=0;l--) {			// test for splitting
	nm=l-1;
	if ((fabs(rv1[l])+anorm) == anorm) {
	  flag=0;
	  break;
	}
	if ((fabs(w(nm,nm))+anorm) == anorm) break;
      }
      if (flag) {
	c=0;
	s=1;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  rv1[i]=(c*rv1[i]);
	  if ((fabs(f)+anorm) == anorm) break;
	  g=w(i,i);
	  h= pythag(f,g);
	  w(i,i)=h;
	  h=1/h;
	  c=g*h;
	  s=(-f*h);
	  for (j=0;j<m;j++) {
	    y=a(j,nm);
	    z=a(j,i);
	    a(j,nm)=(y*c+z*s);
	    a(j,i)=(z*c-y*s);
	  }
	}
      }
      z=w(k,k);
      if (l == k) {				// convergence
	if (z < 0) {				// singular val is made non negative
	  w(k,k) = -z;
	  for (j=0;j<n;j++) v(k,j)=(-v(k,j));
	}
	break;
      }

      if (its == 30) 
      {
        w.setstate(Matrix2D<T>::failbit);
        v.setstate(Matrix2D<T>::failbit);
        return;
      }
       
      x=w(l,l);					// shift from bottom 2x2 minor
      nm=k-1;
      y=w(nm,nm);
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y);
      g=pythag(f,1.);
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
      c=s=1;					// next QR transformation
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w(i,i);
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y=y*c;
	for (jj=0;jj<n;jj++) {
	  x=v(j,jj);
	  z=v(i,jj);
	  v(j,jj)=(x*c+z*s);
	  v(i,jj)=(z*c-x*s);
	}
	z=pythag(f,h);
	w(j,j)=z;
	if (z!=0.0) {
	  z=1/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	for (jj=0;jj<m;jj++) {
	  y=a(jj,j);
	  z=a(jj,i);
	  a(jj,j)=(y*c+z*s);
	  a(jj,i)=(z*c-y*s);
	}
      }
      rv1[l]=0;
      rv1[k]=f;
      w(k,k)=x;
    }
  }
}

#endif
