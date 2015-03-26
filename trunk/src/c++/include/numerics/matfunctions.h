#ifndef MAT_FUNCTIONS_H
#define MAT_FUNCTIONS_H

//
//  Cut Down Matrix Functions.
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   1.1 gcw Initial modified version                Apr  9 96
//             1.2 gcw Cut Down for Relpal                     Jan  6 97
//
//  Copyright (c) 1996  R.C. Elston
//

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
#include "numerics/matrix.h"

//
//  Function - Transpose.  Transposes either a matrix or a vector.
//

template<class A, class B>
inline Matrix2D<B>& Transpose(const Matrix2D<A>& m, Matrix2D<B>& result)
{
  result.resize(m.cols(), m.rows());

      for(size_t i = 0; i < m.rows(); ++i)
        for(size_t j = 0; j < m.cols(); ++j)
          result(j,i) = m(i,j);

  //std::copy(m.begin(), m.end(), result.begin(Col_Order));

  result.setstate(m.rdstate());

  return result;
}

template<class A>
inline Matrix2D<A>& Transpose(Matrix2D<A>& m)
{
  Matrix2D<A> temp;
  Transpose(m, temp);
  m = temp;
  return m;
}

//
//  Function XXT - Returns X * Transpose(X).  Useful for regressions, etc,
//	and can be done much more quickly due to symmetry.
//

template <class A, class B>
inline Matrix2D<B>& XXT(const Matrix1D<A>& v, Matrix2D<B>& s)
{
  if(v.fail()) { s.setstate(MatrixBase::failbit); return s; }

  typename Matrix1D<A>::size_type n = v.size(), i, j;

  if(v.dir() == Row_Order) 
  {
    s.resize(1,1);
    s(0,0) = 0;
    for(i = 0; i < n; ++i)
      s(0,0) += v[i] * v[i];
  }
  else
  {
    s.resize(v.size(),v.size());
    for(i = 0; i < n; ++i)
      for(j = i; j < n; ++j)
        s(j, i) = s(i, j) = v[i] * v[j];
  }
  return s;
}

//inline Matrix2D<double> XXT(const Matrix1D<double>& v)
//{ Matrix2D<double> s; return XXT(v, s); }

//
//  Function XTX - Returns Transpose(X) * X.  Useful for regressions, etc,
//	and can be done much more quickly due to symmetry.
//

template <class A, class B>
inline Matrix2D<B>& XTX(const Matrix2D<A>& m, Matrix2D<B>& s)
{
  if(m.fail()) { s.setstate(MatrixBase::failbit); return s; }

  s.set(0);
  s.resize(m.cols(), m.cols(), 0);

  typename Matrix2D<A>::size_type n = m.cols(), i, j, k;

  for(i = 0; i < n; ++i)
    for(j = i; j < n; ++j)
      for(k = 0; k < m.rows(); k++)
        s(j,i) = s(i,j) += m(k,i) * m(k, j);

  return s;
}

//inline Matrix2D<double> XTX(const Matrix2D<double>& m)
//{ Matrix2D<double> s; return XTX(m, s); }

template <class A, class B, class C>
inline Matrix1D<C>& XTZ(const Matrix2D<A>& m, const Matrix1D<B>& v, Matrix1D<C>& s)
{
  if(m.fail() || v.fail() || v.size() != m.rows())
  { s.setstate(MatrixBase::failbit); return s; }

  s.set(0);
  s.resize(m.cols(), 0);

  typename Matrix2D<A>::size_type i, j;

  for(i = 0; i < m.rows(); ++i)
    for(j = 0; j < m.cols(); ++j)
        s[j] += m(i,j) * v[i];

  return s;
}

//  Function LUPDecompose. Given a n x n matrix A, solves the equation
//      PA = LU, where L is a lower triangular matrix, and U an upper
//      triangular matrix, and P is a permutation matrix of A.  (P is chosen
//      based on a greedy algorithm during the decomposition)

//  Algorithms taken from Introduction to Algorithms, by Thomas H. Cormen,
//    Charles E. Leiserson and Ronald L. Rivest.  McGraw-Hill, 1990. Pgs 754-761

//    struct LUP stores the l, u, and p matrices.  In the case of
//	LUDecompostion, P is returned as an identity matrix.

template <class T>
struct LUP
{
  typedef int state;

  Matrix1D<typename Matrix1D<T>::size_type> p;
  Matrix2D<T>                      l, u;

  LUP(const typename Matrix2D<T>::size_type n=0) : p(n), l(n,n), u(n,n) {}
  LUP(const Matrix2D<T>& m, const double err = 1.e-15)
    { if(!&m || m.fail())
      { setstate(MatrixBase::badbit); return; }
      LUPDecompose(m, *this, err); }
  LUP(const LUP& lup) { p = lup.p; l = lup.l; u = lup.u; }

  LUP& operator () (typename Matrix2D<T>::size_type n = 0)
    { p.resize(n); l.resize(n,n); u.resize(n,n); return *this; }

  LUP& setstate(state s)
    { l.setstate(s); u.setstate(s);
      p.setstate(s); return *this; }
  LUP& clear()
    { l.clear(); u.clear();
      p.clear(); return *this; }
};

template <class A, class B>
inline LUP<B>& LUPDecompose(const Matrix2D<A>& a, LUP<B>& lup, const double err = 1.e-15)
{
  if(!&a || a.fail() || a.rows() != a.cols())
  {
    return lup.setstate(MatrixBase::badbit);
  }

  typename Matrix2D<A>::size_type n = a.rows(), i, j, k, kprime = 0;
  B p;
  B temp;

  lup.clear();
  lup.u = a;
  lup.l.resize(n, n);
  lup.p.resize(n);
  
  for(i = 0; i < n; ++i) lup.p[i] = i;

  for(k = 0; k < n-1; k++)
  {
    p = 0;
    for(i = k; i < n; ++i)
    {
      temp = (lup.u(i,k) > 0) ? (B)lup.u(i,k) : (B)-lup.u(i,k);
      if(temp > p)
      { p = temp; kprime = i; }
    }
    
    if(p == 0)
      return lup.setstate(MatrixBase::failbit);

    j = lup.p[k]; lup.p[k] = lup.p[kprime]; lup.p[kprime] = j;

    for(i = 0; i < n; ++i)
    {
      temp             = lup.u(k, i);
      lup.u(k, i)      = lup.u(kprime, i);
      lup.u(kprime, i) = temp;
    }

    for(i = k + 1; i < n; ++i)
    {
      lup.u(i, k) = lup.u(i, k) / lup.u(k, k);

      for(j = k + 1; j < n; ++j)
        lup.u(i,j) = lup.u(i,j) - lup.u(i, k) * lup.u(k,j);
    }
  }

  for(i = 0; i < n; ++i)
  {
    if(lup.u(i, i) <= err && lup.u(i, i) >= -err)
      return lup.setstate(MatrixBase::failbit);

    for(j = 0; j < n; ++j)
      if(i > j)      { lup.l(i, j) = lup.u(i, j); lup.u(i,j) = 0; }
      else if(i < j) { lup.l(i, j) = 0; }
      else           { lup.l(i, j) = 1; }
  }
  return lup;
}

//  Function LUPSolve - Solves set of linear equations Ax=b where A has been
//    put in LUP form, thus LUx = Pb.
//

//  Algorithms taken from Introduction to Algorithms, by Thomas H. Cormen,
//    Charles E. Leiserson and Ronald L. Rivest.  McGraw-Hill, 1990. Pg. 749-754

template <class A, class B, class C>
inline Matrix1D<C>& LUPSolve(const LUP<A>& lup, const Matrix1D<B>& b, Matrix1D<C>& x)
{
  if(lup.p.fail())
  { x.setstate(MatrixBase::failbit); return x; }

  typename Matrix1D<C>::size_type n = lup.l.rows();

  if(n != b.size()) { x.setstate(MatrixBase::failbit); return x; }

  typename Matrix1D<C>::size_type i, j;
  Matrix1D<C>		 y(n);

  x.resize(n);

  for(i = 0; i < n; ++i)
  {
    y[i] = b[lup.p[i]];
    for(j=0; j < i; ++j)
      y[i] -= lup.l(i, j) * y[j];
  }

  for(i = n-1; i > 0; i--)
  {
    for(j = i+1; j < n; ++j)
      y[i] -= lup.u(i, j) * x[j];
    x[i] = y[i] / lup.u(i, i);
  } 
  for(j = 1; j < n; ++j)
    y[0] -= lup.u(0, j) * x[j];
  x[0] = y[0] / lup.u(0, 0);

  return x;
}

template <class A, class B>
inline Matrix1D<B> 
LUPSolve(const LUP<A>& lup, const Matrix1D<B>& b)
{ 
  Matrix1D<B> x;
  return LUPSolve(lup, b, x);
}

//
// Function Inverse - Given a n x n matrix A, returns the inverse.  
//

template <class A, class B>
inline Matrix2D<B>& Inverse(const Matrix2D<A>& a, Matrix2D<B>& inv,
                          const double err = 0)
{
  if(a.rows() != a.cols()) { inv.setstate(MatrixBase::failbit); return inv; }

  typename Matrix1D<A>::size_type n = a.rows();
  typename Matrix1D<A>::size_type i;

  inv.resize(n, n);

  LUP<A> lup(a, err);

  if(lup.p.fail() || !(lup.u(n-1,n-1)))
  { inv.setstate(MatrixBase::failbit); return inv; }

  Matrix1D<A> en(n, 0), x(n);

  for(i = 0; i < n; ++i)
  {
    en[i] = 1;
    LUPSolve(lup, en, x);
    en[i] = 0;

    std::copy(x.begin(), x.end(), inv.begin(i, Col_Order));
  }

  return inv;
}

template <class A>
inline
 Matrix2D<A> Inverse(const Matrix2D<A>& a, const double err = 0)
{ 
  Matrix2D<A> inv; 
  return Inverse(a, inv, err); 
}

template <class A>
inline A Det(const Matrix2D<A>& a, double err = 0.0)
{
  if(a.rows() != a.cols()) 
    return std::numeric_limits<A>::quiet_NaN();

  typename Matrix1D<A>::size_type n = a.rows();
  typename Matrix1D<A>::size_type i;

  LUP<A> lup(a, err);

  if(lup.p.fail() || !(lup.u(n-1,n-1)))
    return std::numeric_limits<A>::quiet_NaN();

  A det = 1.0;
  for(i=0; i < n; ++i)
    det *= lup.u(i,i);
  return det;
}

template <class A, class B, class C>
inline Matrix2D<C>& Multiply(const Matrix2D<A>& n, const Matrix2D<B>& m, Matrix2D<C>& r)
{
  if(m.rows() != n.cols() || n.fail() || m.fail())
  {  r.setstate(MatrixBase::failbit); return r; }

  r.resize(n.rows(), m.cols());
  r.set(0);
  
  typename Matrix2D<A>::size_type ip = r.rows(), jp = r.cols(), kp = n.cols(); 

  typename Matrix2D<A>::size_type i, j, k;
  
  for(i = 0; i < ip; i++)
    for(j = 0; j < jp; j++)
    {
      typename Matrix2D<A>::value_type s = 0.0;
      for(k = 0; k < kp; k++)
        s += n(i, k) * m(k,j);
      r(i, j) = s;
    }

  return r;
}

template <class A, class B, class C>
inline Matrix2D<C>& XTZ(const Matrix2D<A>& n, const Matrix2D<B>& m, Matrix2D<C>& r)
{
  if(m.rows() != n.rows() || n.fail() || m.fail())
  {  r.setstate(MatrixBase::failbit); return r; }

  r.resize(n.cols(), m.cols());
  r.set(0);
  
  typename Matrix2D<A>::size_type ip = r.rows(), jp = r.cols(), kp = n.rows(); 

  typename Matrix2D<A>::size_type i, j, k;
  
  for(i = 0; i < ip; i++)
    for(j = 0; j < jp; j++)
    {
      typename Matrix2D<A>::value_type s = 0.0;
      for(k = 0; k < kp; k++)
        s += n(k, i) * m(k,j);
      r(i, j) = s;
    }

  return r;
}

template <class A, class B, class C>
inline 
Matrix2D<C>& Multiply(const Matrix2D<A>& m, const Matrix1D<B>& v, Matrix2D<C>& r)
{
  if(m.cols() != ((v.dir() == Col_Order) ? v.size() : 1) ||
     v.fail() || m.fail())
  { r.setstate(MatrixBase::failbit); return r; }

  r.resize(m.rows(), (v.dir()==Row_Order) ? v.size() : 1);
  r.set(0);

  if(v.dir() == Row_Order)
    for(typename Matrix2D<A>::iterator i = r.begin(); i < r.end(); ++i)
      *i = m(i.rows(),0) * v[i.cols()];
  else
    for(typename Matrix2D<A>::const_iterator i = m.begin(); i < m.end(); ++i)
      r(i.rows(),0) += *i * v[i.cols()];

  return r;
}

// FIXME:  This can be done better!
template <class A, class B, class C>
inline 
Matrix2D<C>& Multiply(const Matrix1D<A>& m, const Matrix2D<B>& v, Matrix2D<C>& r)
{
  Matrix2D<A> M(m);

  return Multiply( M, v, r );
}

//  operator + - Matrix + scalar or Matrix + Matrix.

template<class A, class B>
inline Matrix2D<A>& operator += (Matrix2D<A>& m, const Matrix2D<B>& s)
{
  if(m.rows() != s.rows() || m.cols() != s.cols())
  { m.setstate(MatrixBase::failbit); return m; }

  typename Matrix2D<B>::const_iterator i = s.begin();
  for(typename Matrix2D<A>::iterator j = m.begin(); j < m.end(); ++j)
    *j += *i++;
  return m;
}

template <class A, class B>
inline Matrix1D<A>& operator += (Matrix1D<A>& m, const Matrix1D<B>& s)
{
  if(m.size() != s.size())
  { m.setstate(MatrixBase::failbit); return m; }

  typename Matrix1D<A>::const_iterator i = s.begin();
  for(typename Matrix1D<A>::iterator j = m.begin(); j < m.end(); ++j, ++i)
    *j += *i;
  return m;
}

template <class A, class B>
inline Matrix2D<A>& operator *= (Matrix2D<A>& m, B s)
{
  if(!m.rows() || !m.cols())
  { m.setstate(MatrixBase::failbit); return m; }

  for(typename Matrix2D<A>::iterator j = m.begin(); j < m.end(); ++j)
    *j *= s;
  return m;
}

template <class A, class B>
inline Matrix2D<A>& operator /= (Matrix2D<A>& m, B s)
{
  if(!m.rows() || !m.cols())
  { m.setstate(MatrixBase::failbit); return m; }

  for(typename Matrix2D<A>::iterator j = m.begin(); j < m.end(); ++j)
    *j /= s;
  return m;
}

template <class A, class B>
inline Matrix1D<A>& operator -= (Matrix1D<A>& m, const Matrix1D<B>& s)
{
  if(m.size() != s.size())
  { m.setstate(MatrixBase::failbit); return m; }

  typename Matrix1D<B>::const_iterator i = s.begin();
  for(typename Matrix1D<A>::iterator j = m.begin(); j < m.end(); ++j, ++i)
    *j -= *i;
  return m;
}

template <class A, class B>
inline Matrix2D<A>& operator -= (Matrix2D<A>& m, const Matrix2D<B>& s)
{
  if(m.rows() != s.rows() || m.cols() != s.cols())
  { m.setstate(MatrixBase::failbit); return m; }

  typename Matrix2D<B>::const_iterator i = s.begin();
  for(typename Matrix2D<A>::iterator j = m.begin(); j < m.end(); ++j, ++i)
    *j -= *i;
  return m;
}

template <class A, class B>
inline Matrix2D<A> operator + (Matrix2D<A>& m, const Matrix2D<A>& s)
{
  Matrix2D<A> M = m;
  return (M += s);
}

template <class A, class B>
inline Matrix2D<A> operator - (Matrix2D<A>& m, const Matrix2D<B>& s)
{
  Matrix2D<A> M = m;
  return (M -= s);
}

template <class T>
inline void print_matrix(const Matrix2D<T> &w, std::ostream &o)
{
  o << '[';
  for(int i=0; i < w.rows(); ++i)
  {
    if(i)
      o << ';' << std::endl << ' ';

    for(int j=0; j < w.cols(); ++j)
      if( fabs(w(i,j)) < 10e-15)
        o << std::setw(5) << 0.0 << " ";
      else
        o << std::setw(5) << w(i,j) << " ";
  }
  o << " ]" << std::endl;
}

template <class T>
inline void print_vector(const Matrix1D<T> &w, std::ostream &o)
{
  std::string space;
  if(w.dir() == Row_Order)
    space = "  ";
  else
    space = "; ";

  o << "[ ";
  for(size_t i = 0; i < w.size(); ++i)
  {
    if( fabs(w[i]) < 10e-10)
      o << std::setw(5) << 0.0;
    else
      o << std::setw(5) << w[i];
    if(i+1 < w.size()) o << space;
  }
  o << " ]" << std::endl;
}

template <class T>
inline void print_matrix_splus(const Matrix2D<T> &w, std::ostream &o)
{
  o << "matrix(c(";
  for(int i=0; i < w.rows(); ++i)
    for(int j=0; j < w.cols(); ++j)
    {
      if( fabs(w(i,j)) < 10e-10)
        o << std::setw(5) << 0.0;
      else
        o << std::setw(5) << w(i,j);
      if( j+1 < w.cols() || i+1 < w.rows() )
        o << ",";
    }

  o << "),ncol=" << w.cols() << ",byrow=T)" << std::endl;
}

template <class T>
inline void print_vector_splus(const Matrix1D<T> &w, std::ostream &o)
{
  std::string end;

  if(w.dir() == Row_Order)
  {
    o << "t(";
    end = ")";
  }

  o << "c(";
  for(size_t i = 0; i < w.size(); ++i)
  {
    if( fabs(w[i]) < 10e-10)
      o << std::setw(5) << 0.0;
    else
      o << std::setw(5) << w[i];
    if(i+1 < w.size()) 
      o << ",";
  }
  o << ")" << end << std::endl;
}

template <class T>
inline void print_matrix_python(const Matrix2D<T> &w, std::ostream &o)
{
  o << "array([";
  for(int i=0; i < w.rows(); ++i)
  {
    o << "[";
    for(int j=0; j < w.cols(); ++j)
    {
      if( fabs(w(i,j)) < 10e-10)
        o << std::setw(5) << 0.0;
      else
        o << std::setw(5) << w(i,j);
      if( j+1 < w.cols() )
        o << ",";
    }
    o << "]";
    if( i+1 < w.rows() )
      o << ",";
  }
  o << "])" << std::endl;
}

template <class T>
inline void print_vector_python(const Matrix1D<T> &w, std::ostream &o)
{
  std::string end;

  if(w.dir() == Row_Order)
  {
    o << "t(";
    end = ")";
  }

  o << "array([";
  for(size_t i = 0; i < w.size(); ++i)
  {
    if( fabs(w[i]) < 10e-10)
      o << std::setw(5) << 0.0;
    else
      o << std::setw(5) << w[i];
    if(i+1 < w.size()) 
      o << ",";
  }
  o << "])" << end << std::endl;
}



#endif
