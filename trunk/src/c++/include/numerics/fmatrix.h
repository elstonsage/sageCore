#ifndef FORTRAN_MATRIX
#define FORTRAN_MATRIX

#include <iostream>
#include <iomanip>
#include <string>
#include <memory>
#include <valarray>
#include "globals/config.h"
#include "numerics/matrix.h"

namespace SAGE {

// Currently this is triggered spuriously during certain
// transpose and resize operations.
#undef CHECK_BOUNDS

template<class T>
class FortranMatrix : public MatrixBase
{
  // Implements Fortran-style 2-dimensional matrices.
  //
  // -----------------------------
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |------+------+------+------|
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // |      |      |      |      | 
  // -----------------------------

public:

  typedef MatrixBase base_type;

  typedef std::vector<T>                  storage_type;
  
  typedef typename storage_type::value_type        value_type;
  typedef typename storage_type::size_type         size_type;
  typedef typename storage_type::reference         reference;
  typedef typename storage_type::const_reference   const_reference;
  typedef typename storage_type::pointer           pointer;
  typedef typename storage_type::const_pointer     const_pointer;
  typedef typename storage_type::difference_type   difference_type;

private:

  storage_type      my_storage;      // Vector of storage
  size_type         my_rows;         // Rows currently stored
  size_type         my_cols;         // Columns currently stored
  size_type         my_lda;          // Leading dimension (max column length)

public:
  FortranMatrix()
           : my_storage(), my_rows(0), my_cols(0), my_lda(0) { }

  explicit
  FortranMatrix( size_type n ) 
           : my_storage(n*n), my_rows(n), my_cols(n), my_lda(n) { }

  FortranMatrix( size_type m, size_type n) 
           : my_storage(m*n, T()), my_rows(m), my_cols(n), my_lda(m) { }

  FortranMatrix( size_type m, size_type n, size_type lda, const T& value = T() ) 
           : my_storage(n*std::max(m,lda), value), my_rows(m), my_cols(n), 
                   my_lda(std::max(m,lda)) { }

  FortranMatrix(const FortranMatrix& x) : my_storage(x.my_storage),
         my_rows(x.my_rows), my_cols(x.my_cols), my_lda(x.my_lda) { }
  
 ~FortranMatrix() { }

#ifndef BROKEN_TEMPLATE_FRIENDS
  template<class X> friend class FortranMatrix;
#endif

  template<class X>
  FortranMatrix<T>& operator=(const FortranMatrix<X>& x)
  {
    size_type size_needed = x.my_rows * x.my_cols;
    
    if( size_needed > my_storage.size() )
    {
      // We do need to realloc
      my_storage = x.my_storage;
      my_rows    = x.my_rows;
      my_cols    = x.my_cols;
      my_lda     = x.my_lda;
    }
    else
    {
      // We don't need to realloc
      my_rows = x.my_rows;
      my_cols = x.my_cols;
      my_lda  = my_storage.size()/my_cols;

      assert( my_lda >= my_rows && my_cols*my_lda <= my_storage.size() );
      
      for(size_type j = 0; j < my_cols; ++j)
        for(size_type i = 0; i < my_rows; ++i)
          (*this)(i,j) = x(i,j);
    }
      
    return *this;
  }

  // Be real careful with these!
        T*   raw_storage()       { return &my_storage[0];    }
  const T*   raw_storage() const { return &my_storage[0];    }

//size_type         size() const { return my_storage.size(); }
  size_type         rows() const { return my_rows;           }
  size_type         cols() const { return my_cols;           }
  size_type          lda() const { return my_lda;            }
  size_type row_capacity() const { return my_lda;            }
  size_type col_capacity() const { return my_cols;           }
  bool             empty() const { return (my_rows == 0 || my_cols == 0); }


  void          resize_fill(size_type m, size_type n, T c = T() )
  {  
      // FIXME: LDN should be a parameter
      size_type ldn = m;
              
      size_type new_size = n*ldn;
      if(new_size > my_storage.size() )
        my_storage.resize(new_size, c);
      std::fill(my_storage.begin(), my_storage.end(), c);
      my_rows = m;
      my_cols = n;
      if(my_cols > 0)
        my_lda  = my_storage.size()/my_cols;
      else
        my_lda = 0;
  }

  void          resize_nofill(size_type m, size_type n, size_type ldn = 0)
  {  
      ldn = std::max(m, ldn);
              
      size_type new_size = n*ldn;
      size_type old_size = my_storage.size();
      if(new_size > old_size )
        my_storage.resize(new_size);
      my_rows = m;
      my_cols = n;
      if(my_cols > 0)
        my_lda  = my_storage.size()/my_cols;
      else
        my_lda  = 0;
  }

  void          resize(size_type m, size_type n, T c = T() )
  {
      // If the rows requested is greater than our leading dimension
      // we must resize.
      // NOTE: This section is exception safe since everything is done
      //       on a temporary until the last step.
      if( m > lda() )
      {
        FortranMatrix<T> temp(m, n, m, c);
        extract_to(temp, std::min(m,rows()), std::min(n,cols()));
        swap(temp);
        return;
      }

      // Add requested columns 
      // NOTE: This step should be done first for exception safety.
      // ALSO: This step should be smarter!
      if( n > cols() && lda()*n > my_storage.size() )
        my_storage.resize( lda() * n, c );

      // Grow already allocated rows and initialize the entries
      // NOTE: This section should not throw, but if it does,
      //       we are still safe since the affected elements must be 
      //       re-initialized before they can be used.
      if( m > rows() )
      {
        size_type old_rows = rows();
        my_rows = m;
        for(size_type j = 0; j < cols(); ++j)
          for(size_type i = old_rows; i < rows(); ++i)
            (*this)(i,j) = c;
      }

      my_rows = m;
      my_cols = n;
  }

  // Reserve row capacity: i.e. resize lda
/*
// This function has been commented out because the instantiation
// of FortranMatrix<T> uses the parameter 'c', but 'c' is not available.
// This is probably a legitimate bug, but fixing it currently is a bit
// complicated. --sgross 31-oct-05

  void reserve(size_type m)
  {
      if( m > lda() )
      {
        FortranMatrix<T> temp( rows(), cols(), m, T());
        temp.update( *this, 0, 0, rows(), cols() );
        swap(temp);
        return;
      }
  }
*/
  
        reference operator() (size_type m, size_type n)       
  { 
#ifdef CHECK_BOUNDS
    if(m >= my_lda || m >= rows() || n >= cols() || my_lda*n+m >= my_storage.size())
    {
      cout << "Bounds error: " << m << "x" << n << "." << endl 
           << "        this: " << rows() << "x" << cols() << ", lda = " << lda() 
           << ", offset = " << my_lda*n+m << ", size = " << my_storage.size() << endl;
    }
#endif
               
    //assert(my_lda > m && m < rows() && n < cols() && (my_lda*n+m) < my_storage.size());
    return my_storage[my_lda*n+m]; 
  }
  const_reference operator() (size_type m, size_type n) const 
  { 
#ifdef CHECK_BOUNDS
    if(m >= my_lda || m >= rows() || n >= cols() || my_lda*n+m >= my_storage.size())
    {
      cout << "Bounds error: " << m << "x" << n << "." << endl 
           << "        this: " << rows() << "x" << cols() << ", lda = " << lda() 
           << ", offset = " << my_lda*n+m << ", size = " << my_storage.size() << endl;
    }
#endif
    //assert(my_lda > m && m < rows() && n < cols() && (my_lda*n+m) < my_storage.size());
    return my_storage[my_lda*n+m]; 
  }

        reference at(size_type m, size_type n)       { return my_storage.at(my_lda*n+m); }
  const_reference at(size_type m, size_type n) const { return my_storage.at(my_lda*n+m); }
  
  void swap(FortranMatrix<T>& x)
  {
    my_storage.swap(x.my_storage);
    std::swap(my_rows, x.my_rows);
    std::swap(my_cols, x.my_cols);
    std::swap(my_lda,  x.my_lda);
  }
  
  FortranMatrix<T> extract(size_type rows, size_type cols,
                                      size_type top = 0, size_type left = 0) const
  {
    FortranMatrix<T> temp(rows-top, cols-left);
    extract_to(temp, rows, cols, top, left);
    return temp;
  }
  
  void extract_to(FortranMatrix<T>& x, size_type erows, size_type ecols,
                  size_type top = 0, size_type left = 0) const
  {
    const size_type r = std::min(x.rows(), std::min( rows(), top  + erows ));
    const size_type c = std::min(x.cols(), std::min( cols(), left + ecols ));
    
    for(size_type j = left; j < c; ++j)
      for(size_type i = top; i < r; ++i)
        x( i - top, j - left ) = (*this)(i,j);
  }

  void copy_to(FortranMatrix<T>& x, size_type erows, size_type ecols,
                  size_type top = 0, size_type left = 0) const
  {
    const size_type r = std::min(rows(), std::min( x.rows(), top  + erows ));
    const size_type c = std::min(cols(), std::min( x.cols(), left + ecols ));
    
    for(size_type j = 0; j < c; ++j)
      for(size_type i = 0; i < r; ++i)
        x( i + top, j + left ) = (*this)(i,j);
  }
  
  void fill(T c = T())
  {
    my_storage.assign(my_storage.size(), c);
  }

  // Unary operators
  FortranMatrix operator-()
  {
    FortranMatrix temp = (*this);
    for(size_type j = 0; j < cols(); ++j)
      for(size_type i = 0; i < rows(); ++i)
        temp(i,j) = -(*this)(i,j);
    return temp;
  }

  // Computed assignment
  FortranMatrix& operator+= (const T& t) { return binary_assign_op< std::plus<T>        >(t); }
  FortranMatrix& operator-= (const T& t) { return binary_assign_op< std::minus<T>       >(t); }
  FortranMatrix& operator*= (const T& t) { return binary_assign_op< std::multiplies<T>  >(t); }
  FortranMatrix& operator/= (const T& t) { return binary_assign_op< std::divides<T>     >(t); }
  FortranMatrix& operator%= (const T& t) { return binary_assign_op< std::modulus<T>     >(t); }

  FortranMatrix& operator+= (const FortranMatrix& v) { return binary_assign_op< std::plus<T>  >(v); }
  FortranMatrix& operator-= (const FortranMatrix& v) { return binary_assign_op< std::minus<T> >(v); }
  FortranMatrix& operator*= (const FortranMatrix& v) 
  {
    FortranMatrix temp;
    multiply(*this, v, temp);
    swap(temp);
    return *this;
  }

  FortranMatrix& transpose()
  {
      if( !(*this) )
        return *this;

      // Case 1:  Fast transposes where memory layout does not change
      // Case 1a:   Fast transpose for a row vector with lda == 1:
      // Case 1b:   Fast transpose for a column vector:
      //
      //                            |------|
      //                            | x ...|
      // -----------------          | x ...|
      // |x x x x x  ....|  <===>   | x ...|
      // -----------------          | x ...|
      //                            | x ...|
      //                            |------|
      //
      // Case 2: Simple/direct in-place transpose using excess space between
      //         cols and lda
      // Case 3: Expensive re-allocate and copy

      if( cols() < 2 )
      {
        std::swap(my_rows, my_cols);
        my_lda = 1;
      }
      else if( lda() < 2 )
      {
        std::swap(my_rows, my_cols);
        my_lda = my_rows;
      }
#if 0
      else if( cols() <= lda() )
      {
        for(size_type j = 0; j < cols(); ++j)
          for(size_type i = j+1; i < rows(); ++i)
            std::swap( (*this)(i,j), (*this)(j,i) );
        std::swap(my_rows, my_cols);
      }
#endif
      else
      {
        FortranMatrix temp(cols(), rows());

        for(size_type j = 0; j < cols(); ++j)
          for(size_type i = 0; i < rows(); ++i)
            temp(j,i) = (*this)(i,j);

        swap(temp);
      }
      return *this;
  }

private:
  template<class X,class F> 
  FortranMatrix( const FortranMatrix<X>& v, const F )
  {
    if(!v)
    {
      setstate(failbit);
      return;
    }

    resize_nofill(v.rows(), v.cols());

    for(size_type j = 0; j < cols(); ++j)
      for(size_type i = 0; i < rows(); ++i)
        (*this)(i,j) = f(v(i,j));
  }

  template<class X,class F> 
  FortranMatrix( const FortranMatrix<X>& u, const FortranMatrix<X>& v, const F f)
  {
    if(!v || !u || v.rows() != u.rows() || v.cols() != u.cols() )
    {
      setstate(failbit);
      return;
    }

    resize_nofill(v.rows(), v.cols());

    for(size_type j = 0; j < cols(); ++j)
      for(size_type i = 0; i < rows(); ++i)
        (*this)(i,j) = f(u(i,j), v(i,j));
  }

  template<class X,class F> 
  FortranMatrix( const FortranMatrix<X>& v, const X& x, const F f)
  {
    if(!v)
    {
      setstate(failbit);
      return;
    }

    resize_nofill(v.rows(), v.cols());

    for(size_type j = 0; j < cols(); ++j)
      for(size_type i = 0; i < rows(); ++i)
        (*this)(i,j) = f(v(i,j), x);
  }

  template<class X,class F> 
  FortranMatrix( const X& x, const FortranMatrix<X>& v, const F f)
  {
    if(!v)
    {
      setstate(failbit);
      return;
    }

    resize_nofill(v.rows(), v.cols());

    for(size_type j = 0; j < cols(); ++j)
      for(size_type i = 0; i < rows(); ++i)
        (*this)(i,j) = f(x, v(i,j));
  }

  template<class F>
  inline 
  FortranMatrix& binary_assign_op(const T& t) 
  {
      if( !(*this) )
        return *this;
        
      for(size_type j = 0; j < cols(); ++j)
        for(size_type i = 0; i < rows(); ++i)
          (*this)(i,j) = F()( (*this)(i,j), t);

      return *this;
  }  

  template<class F>
  inline 
  FortranMatrix& binary_assign_op(const FortranMatrix& v)
  {
      if( !(*this) || !v || rows() != v.rows() || cols() != v.cols() )
      {
        setstate(failbit);
        return *this;
      }        

      for(size_type j = 0; j < cols(); ++j)
        for(size_type i = 0; i < rows(); ++i)
          (*this)(i,j) = F()( (*this)(i,j), v(i,j) );

      return *this;
  }

  template <class TT>
  friend FortranMatrix<TT> operator*(const FortranMatrix<TT>&, const FortranMatrix<TT>&);
  template <class TT>
  friend FortranMatrix<TT> operator*(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<TT> operator*(const TT&, const FortranMatrix<TT>&);

  template <class TT>
  friend FortranMatrix<TT> operator/(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<TT> operator/(const TT&, const FortranMatrix<TT>&);

  template <class TT>
  friend FortranMatrix<TT>& 
  multiply(const FortranMatrix<TT>& v1, const FortranMatrix<TT>& v2, FortranMatrix<TT>& x);

  template <class TT>
  friend FortranMatrix<TT> operator+(const FortranMatrix<TT>&, const FortranMatrix<TT>&);
  template <class TT>
  friend FortranMatrix<TT> operator+(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<TT> operator+(const TT&, const FortranMatrix<TT>&);
  template <class TT>

  friend FortranMatrix<TT> operator-(const FortranMatrix<TT>&, const FortranMatrix<TT>&);
  template <class TT>
  friend FortranMatrix<TT> operator-(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<TT> operator-(const TT&, const FortranMatrix<TT>&);

#if NOT_YET
  template <class TT>
  friend FortranMatrix<bool> operator==(const FortranMatrix<TT>&, const FortranMatrix<TT>&);
  template <class TT>
  friend FortranMatrix<bool> operator==(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<bool> operator==(const TT&, const FortranMatrix<TT>&);

  template <class TT>
  friend FortranMatrix<bool> operator!=(const FortranMatrix<TT>&, const FortranMatrix<TT>&);
  template <class TT>
  friend FortranMatrix<bool> operator!=(const FortranMatrix<TT>&, const TT&);
  template <class TT>
  friend FortranMatrix<bool> operator!=(const TT&, const FortranMatrix<TT>&);
#endif

  template <class TT>
  friend FortranMatrix<TT>& transpose(const FortranMatrix<TT>&, FortranMatrix<TT>&);
};

template <class T>
inline FortranMatrix<T> 
operator+(const FortranMatrix<T>& v, const T& t) 
{ 
  return FortranMatrix<T>(v, t, std::plus<T>()); 
}

template <class T>
inline FortranMatrix<T> 
operator+(const T& t, const FortranMatrix<T>& v) 
{ 
  return FortranMatrix<T>(t, v, std::plus<T>());
}

template <class T>
inline FortranMatrix<T> 
operator+(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2) 
{ 
  return FortranMatrix<T>(v1, v2, std::plus<T>()); 
}

template <class T>
inline FortranMatrix<T> 
operator-(const FortranMatrix<T>& v, const T& t) 
{ 
  return FortranMatrix<T>(v, t, std::minus<T>()); 
}

template <class T>
inline FortranMatrix<T> 
operator-(const T& t, const FortranMatrix<T>& v) 
{ 
  return FortranMatrix<T>(t, v, std::minus<T>());
}

template <class T>
inline FortranMatrix<T> 
operator-(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2) 
{ 
  return FortranMatrix<T>(v1, v2, std::minus<T>()); 
}

template <class T>
inline FortranMatrix<T> 
operator*(const T& t, const FortranMatrix<T>& v) 
{ 
  return FortranMatrix<T>(t, v, std::multiplies<T>());
}

template <class T>
inline FortranMatrix<T> 
operator*(const FortranMatrix<T>& v, const T& t) 
{ 
  return FortranMatrix<T>(v, t, std::multiplies<T>()); 
}

template <class T>
inline FortranMatrix<T> 
operator/(const T& t, const FortranMatrix<T>& v) 
{ 
  return FortranMatrix<T>(t, v, std::divides<T>());
}

template <class T>
inline FortranMatrix<T> 
operator/(const FortranMatrix<T>& v, const T& t) 
{ 
  return FortranMatrix<T>(v, t, std::divides<T>()); 
}

template <class T>
inline FortranMatrix<T>& 
multiply(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2, FortranMatrix<T>& o) 
{
  if( !v1 || !v2 || v1.cols() != v2.rows() )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(v1.rows(), v2.cols());

  typename FortranMatrix<T>::size_type i, j, k;

  for(i = 0; i < o.rows(); ++i)
    for(j = 0; j < o.cols(); ++j)
    {
      T dot(0.0);
      for(k = 0; k < v1.cols(); ++k)
        dot += v1(i,k)*v2(k,j);
      o(i,j) = dot;
    } 
  return o;
}

template <class T>
inline FortranMatrix<T>& 
XZT(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2, FortranMatrix<T>& o) 
{
  if( !v1 || !v2 || v1.cols() != v2.cols() )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(v1.rows(), v2.rows());

  typename FortranMatrix<T>::size_type i, j, k;

  for(i = 0; i < o.rows(); ++i)
    for(j = 0; j < o.cols(); ++j)
    {
      T dot = 0.0;
      for(k = 0; k < v1.cols(); ++k)
        dot += v1(i,k)*v2(j,k);
      o(i,j) = dot;
    } 
  return o;
}

template <class T>
inline FortranMatrix<T>& 
XTZ(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2, FortranMatrix<T>& o) 
{
  if( !v1 || !v2 || v1.rows() != v2.rows() )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(v1.cols(), v2.cols());

  typename FortranMatrix<T>::size_type i, j, k;

  for(i = 0; i < o.rows(); ++i)
    for(j = 0; j < o.cols(); ++j)
    {
      T dot = 0.0;
      for(k = 0; k < v1.rows(); ++k)
        dot += v1(k,i)*v2(k,j);
      o(i,j) = dot;
    } 
  return o;
}

template <class T>
inline FortranMatrix<T>& 
XTX(const FortranMatrix<T>& v, FortranMatrix<T>& o) 
{
  if( !v )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(v.cols(), v.cols());

  typename FortranMatrix<T>::size_type i, j, k;

  for(i = 0; i < o.rows(); ++i)
    for(j = 0; j < o.cols(); ++j)
    {
      T dot = 0.0;
      for(k = 0; k < v.rows(); ++k)
        dot += v(k,i)*v(k,j);
      o(i,j) = dot;
    } 
  return o;
}

template <class T>
inline FortranMatrix<T>& 
XXT(const FortranMatrix<T>& v, FortranMatrix<T>& o) 
{
  if( !v )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(v.rows(), v.rows());

  typename FortranMatrix<T>::size_type i, j, k;

  for(i = 0; i < o.rows(); ++i)
    for(j = 0; j < o.cols(); ++j)
    {
      T dot(0.0);
      for(k = 0; k < v.cols(); ++k)
        dot += v(i,k)*v(j,k);
      o(j,i) = dot;
    } 
  return o;
}

template <class T>
inline FortranMatrix<T> 
operator*(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2) 
{
  FortranMatrix<T> temp;
  multiply(v1,v2,temp);
  return temp;
}

template<class T>
inline
FortranMatrix<T>& transpose(const FortranMatrix<T>& m, FortranMatrix<T>& o)
{
  if( !m )
  {
    o.setstate(FortranMatrix<T>::failbit);
    return o;
  }

  o.resize_nofill(m.cols(), m.rows());
  typename FortranMatrix<T>::size_type i,j;

  for(j = 0; j < m.cols(); ++j)
    for(i = 0; i < m.rows(); ++i)
      o(j,i) = m(i,j);

  return o;
}

template<class T>
inline
FortranMatrix<T> transpose(const FortranMatrix<T>& m)
{
  FortranMatrix<T> temp;
  transpose(m, temp);
  return temp;
}

template <class T>
inline FortranMatrix<T> 
zero(size_t n)
{
  FortranMatrix<T> temp(n,n,0);
  return temp;
}

template <class T>
inline FortranMatrix<T> 
zero(size_t n, size_t m)
{
  FortranMatrix<T> temp(n,m,0);
  return temp;
}

template <class T>
inline FortranMatrix<T> 
eye(size_t n)
{
  FortranMatrix<T> temp(n,n,0);
  for(size_t i = 0; i < n; ++i)
    temp(i,i) = 1;
  return temp;
}

template <class T>
inline FortranMatrix<T> 
kron(const FortranMatrix<T>& v1, const FortranMatrix<T>& v2) 
{
  FortranMatrix<T> k(v1.rows()*v2.rows(),v1.cols()*v2.cols());
  FortranMatrix<T> temp(v2.rows(),v2.cols());
  
  for(size_t j = 0; j < v1.cols(); ++j)
    for(size_t i = 0; i < v1.rows(); ++i)
    {
      temp  = v2;
      temp *= v1(i,j);
      temp.copy_to(k, v2.rows(), v2.cols(), i*v2.rows(), j*v2.cols());
    }
  return k;
}

class SVD
{
public:
  SVD() : threshold(1e-6) { }
  SVD(const FortranMatrix<double>& A) : threshold(1e-6) { compute(A); }
  bool fail() const                   { return !(U && S && Vt); }

  size_t compute(const FortranMatrix<double>& A);

  FortranMatrix<double>  inverse_of(const FortranMatrix<double>& A);
  FortranMatrix<double> &inverse_of(const FortranMatrix<double>& A, FortranMatrix<double>& Ai);
  FortranMatrix<double>  inverse();
  FortranMatrix<double> &inverse(FortranMatrix<double>& Ai);

  FortranMatrix<double> &sqrt(FortranMatrix<double>& S);
  FortranMatrix<double> &inverse_sqrt(FortranMatrix<double>& Si);

  FortranMatrix<double>  general_inverse_of(const FortranMatrix<double>& A);
  FortranMatrix<double> &general_inverse_of(const FortranMatrix<double>& A, FortranMatrix<double>& Ai);
  FortranMatrix<double>  general_inverse();
  FortranMatrix<double> &general_inverse(FortranMatrix<double>& Ai);

  FortranMatrix<double> &general_sqrt(FortranMatrix<double>& S);
  FortranMatrix<double> &general_inverse_sqrt(FortranMatrix<double>& Si);

  double threshold;
  
  double determinant() const;
  double well_condition() const;

  FortranMatrix<double> U;
  FortranMatrix<double> S;
  FortranMatrix<double> Vt;
  FortranMatrix<double> temp;

  std::vector<size_t> positive_diagonal;
    
private:
  std::vector<double> work;
};  

FortranMatrix<double>
SVDinverse(const FortranMatrix<double>& A);

FortranMatrix<double>& 
SVDinverse(const FortranMatrix<double>& A, FortranMatrix<double>& Ai);

inline
double well_condition(const FortranMatrix<double>& A)
{
  SVD svd(A);
  return svd.well_condition();
}

class GQR
{
public:
  GQR() { }
  GQR(const FortranMatrix<double>& A, const FortranMatrix<double>& B, bool reduced = true) 
  { 
    compute(A,B, reduced); 
  }

  bool fail() const { return !(Q && R && Z && T); }

  void compute(const FortranMatrix<double>& A, const FortranMatrix<double>& B, bool reduced = true);

  FortranMatrix<double> Q;
  FortranMatrix<double> R;
  FortranMatrix<double> T;
  FortranMatrix<double> Z;
    
private:
  FortranMatrix<double> AA,BB;
  std::vector<double> work;
};  

class LUP
{
public:
  LUP() { }
  LUP(const FortranMatrix<double>& A) { compute(A); }
  bool fail() const                   { return !LU.good(); }

  void compute(const FortranMatrix<double>& A);
  FortranMatrix<double>  inverse_of(const FortranMatrix<double>& A);
  FortranMatrix<double> &inverse_of(const FortranMatrix<double>& A, FortranMatrix<double>& Ai);
  FortranMatrix<double>  inverse();
  FortranMatrix<double> &inverse(FortranMatrix<double>& Ai);

  double determinant() const;
  
  std::vector<int> pivots;
  FortranMatrix<double> LU;
    
private:
  std::vector<double> work;
};  

FortranMatrix<double>
LUinverse(const FortranMatrix<double>& A);

FortranMatrix<double>& 
LUinverse(const FortranMatrix<double>& A, FortranMatrix<double>& Ai);

FortranMatrix<double>& 
Cholesky(const FortranMatrix<double>& A, FortranMatrix<double>& UL, bool upper = true);

FortranMatrix<double>& 
TRIinverse(const FortranMatrix<double>& UL, FortranMatrix<double>& ULi, bool upper = true);

FortranMatrix<double>
TRIinverse(const FortranMatrix<double>& UL, bool upper = true);

FortranMatrix<double>& 
SYMinverse(const FortranMatrix<double>& S, FortranMatrix<double>& Si);

FortranMatrix<double>
SYMinverse(const FortranMatrix<double>& S);

FortranMatrix<double>& 
Msqrt(const FortranMatrix<double>& XX, FortranMatrix<double>& X, SVD& svd);

FortranMatrix<double>
Msqrt(const FortranMatrix<double>& XX);

FortranMatrix<double>& 
Minverse_sqrt(const FortranMatrix<double>& XX, FortranMatrix<double>& X, SVD& svd);

FortranMatrix<double>
Minverse_sqrt(const FortranMatrix<double>& XX);

FortranMatrix<double>& 
forward_substitute(FortranMatrix<double>& L, FortranMatrix<double> &B);

FortranMatrix<double>& 
back_substitute(FortranMatrix<double>& U, FortranMatrix<double> &B);

int
Eigen(FortranMatrix<double>& b, std::vector<double>& w);

template <class T>
inline bool matrix_equal(const FortranMatrix<T> &m1, const FortranMatrix<T> &m2, 
                         T epsilon = std::numeric_limits<T>::epsilon())
{
  if(!m1 || !m2 || m1.rows() != m2.rows() || m1.cols() != m2.cols())
    return false;

  for(size_t j = 0; j < m1.cols(); ++j)
    for(size_t i = 0; i < m1.rows(); ++i)
      if( abs( m1(i,j) - m2(i,j) ) > epsilon )
        return false;
  return true;
}

template <class T>
inline void print_matrix(const FortranMatrix<T> &w, std::ostream &o, std::string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';


  for(int i=0; i < (int)w.rows(); ++i)
  {
    if(i)
      o << ';' << std::endl << std::setw(offset) << "";

    for(int j=0; j < (int)w.cols(); ++j)
      if( fabs(w(i,j)) < 10e-15)
        o << std::setw(5) << 0.0 << " ";
      else
        o << std::setw(5) << w(i,j) << " ";
  }
  o << " ]" << std::endl;
}

template <class T>
inline void print_matrix_splus(const FortranMatrix<T> &w, std::ostream &o, std::string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';


  for(int i=0; i < (int)w.rows(); ++i)
  {
    if(i)
      o << ';' << std::endl << std::setw(offset) << "";

    for(int j=0; j < (int)w.cols(); ++j)
      if( fabs(w(i,j)) < 10e-15)
        o << std::setw(5) << 0.0 << " ";
      else
        o << std::setw(5) << w(i,j) << " ";
  }
  o << " ]" << std::endl;
}
                                                                        

template <class T>
inline void print_matrix_python(const FortranMatrix<T> &w, std::ostream &o, std::string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';


  for(int i=0; i < (int)w.rows(); ++i)
  {
    if(i)
      o << ';' << std::endl << std::setw(offset) << "";

    for(int j=0; j < (int)w.cols(); ++j)
      if( fabs(w(i,j)) < 10e-15)
        o << std::setw(5) << 0.0 << " ";
      else
        o << std::setw(5) << w(i,j) << " ";
  }
  o << " ]" << std::endl;
}

template <class T>
inline void print_matrix_first10(const FortranMatrix<T> &w, std::ostream &o, std::string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';

  int row = 10;
  if( (int)w.rows() < 10 )
    row = (int)w.rows();

  int col = 10;
  if( (int)w.cols() < 10 )
    col = (int)w.cols();

  for(int i=0; i < row; ++i)
  {
    if(i)
      o << ';' << std::endl << std::setw(offset) << "";

    for(int j=0; j < col; ++j)
      if( fabs(w(i,j)) < 10e-15)
        o << std::setw(5) << 0.0 << " ";
      else
        o << std::setw(5) << w(i,j) << " ";
  }
  o << " ]" << std::endl;
}

}

#endif
