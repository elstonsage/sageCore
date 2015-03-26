#ifndef SAGE_MATRIX_H
#define SAGE_MATRIX_H

//
//  matrix.h - Basic Matrices. 
//
//  Author: Geoff Wedig  (wedig@darwin.cwru.edu)
//          Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.1  gcw Final modifications before LSF.       Mar 12 1996
//             1.2  gcw Major modifications.  New functions
//                      for matrices, renamed many classes,
//                      fixed problems, simplified
//                      iterators and made LSF compatible     Sep 25 1996
//             1.3  gcw final Matrix implementation for
//                      Relpal 1.0                            Dec 12 1996
//             1.99 kbj Major update before converting to
//                      FORTRAN storage format.               Jun  4 2000
//
//  Copyright (c) 1996  R.C. Elston
//

#include <iterator>
#include <algorithm>
#include <utility>
#include <vector>
#include <functional>
#include <cmath>
#include <cassert>

//  Notes:
//    1.  This file includes the basic 1 and 2 dimensional matrix classes
//        (Matrix1D and Matrix2D respectively) and the iterators for them.
//        Short of the boolean comparisons all operators (multiplication,
//        addition) and functions are included in in the file matfunc.h. 
//        Special matrix forms (identity, permutation matrices, etc) are
//        included in specmat.h.
//    2.  While the iterators provided for the 2 dimensional matrices are
//        reasonably efficient, there is some additional overhead.  The
//        extra overhead is neglible in most cases, but in cases where
//        indices are used for several matrices, or repeated coursing over
//        the matrix is required (for example, matrix multiplication
//        requires several repeats over the matrices involved, and indices
//        are concurrent for several matrices) it may be better to use
//        standard indices.

//
//  Base defines for matrices
//

enum __matrix_dir { Row_Order = 0, Col_Order = 1 };
enum __iter_type  { Forward_Iter = 1, Reverse_Iter = 0 };

//
// MatrixBase is the basic error handling of a matrix.  If we are using the
//   LSF for our matrices, we would like the matrices to use the common LSF
//   error handling (our matrices will gain many other benefits as well, but
//   these are the ones we are interested in right now)
//

class MatrixBase
{
public:

// Constructor/Destructor

  MatrixBase();
 ~MatrixBase();

  enum m_state { goodbit, failbit, badbit };

  void     clear();
  void     clear(const int state);
  void  setstate(const int flag);

  int       good() const;
  int       fail() const;
  int        bad() const;
  int    rdstate() const;
  operator void*() const;
  int  operator!() const;

private:

  int             _state;              // Flags
};

//
//  class MatrixDefines - Basic storage dependent typedefs and definitions
//      of terms necessary for all matrices.
//

template <class T>
class MatrixDefines
{
public:
  typedef T                    value_type;
  typedef typename std::vector<T>::size_type size_type;
  typedef ptrdiff_t            distance;
  typedef T&                   reference;
  typedef const T&             const_reference;
};

//
//  class LimitedVector - A limited vector class.  Provides raw storage with
//      limited functionality. We only allow those functions of vectors we
//      wish. This is mainly to keep tighter control on the allocation of
//      memory.  A number of functions have also been added to make the code
//      for 2D matrices easier to write (minimize, insert, etc.)
//

template <class T>
class LimitedVector
{
public:

  typedef typename MatrixDefines<T>::value_type            value_type;
  typedef typename MatrixDefines<T>::size_type             size_type;
  typedef typename MatrixDefines<T>::distance              distance;
  typedef typename MatrixDefines<T>::reference             reference;
  typedef typename MatrixDefines<T>::const_reference       const_reference;
  typedef typename std::vector<T>::iterator                iterator;
  typedef typename std::vector<T>::const_iterator          const_iterator;
  typedef typename std::vector<T>::reverse_iterator        reverse_iterator;
  typedef typename std::vector<T>::const_reverse_iterator  const_reverse_iterator;

// Constructors/Destructor

  LimitedVector()                                  : array(0,T())   {};
  LimitedVector(size_type n)                       : array(n,T())   {};
  LimitedVector(size_type n, const T& t)           : array(n,t)     {};
  LimitedVector(const LimitedVector<T>& a)         : array(a.array) {};
  LimitedVector(const std::vector<T>& v)           : array(v)       {};

  ~LimitedVector() {};

// Input

  LimitedVector& operator = (const LimitedVector<T>& a)
    { array = a.array; return *this; }

  void resize(const size_type i, const T& t = T())
    { array.reserve(i);
      if (i > array.size()) array.insert(array.end(),(size_t)i - array.size(), t);
      if (i < array.size()) array.erase(array.begin() + i,array.end());    }

  void erase(const size_type i)
    { if (i <  array.size()) array.erase(begin()+i);      }

  void insert(const size_type i, const T& t = T())
    { array.reserve(array.size() + 1);
      if (i <= array.size()) array.insert(begin()+i, t);  }

  void minimize()
    { if (array.capacity() > array.size())
      { std::vector<T> v(array.size()); 
        iterator i = array.begin();
        iterator j = v.begin();
        for( ; i != array.end(); ++i, ++j)
           *i = *j;
        v.swap(array);                                  } } 

  void set(const T& t = T())
    { for (iterator i=begin(); i<end(); ++i) *i = t;    }

  void swap(LimitedVector<T>& m2) { array.swap(m2.array); }
    
// Output

  size_type     size () const { return array.size();      }
  size_type capacity () const { return array.capacity();  }

// Iterators

  iterator               begin()              { return array.begin();  }
  const_iterator         begin()      const   { return array.begin();  }
  iterator                 end()              { return array.end();    }
  const_iterator           end()      const   { return array.end();    }

  reverse_iterator       rbegin()             { return array.rbegin(); }
  const_reverse_iterator rbegin()     const   { return array.rbegin(); }
  reverse_iterator         rend()             { return array.rend();   }
  const_reverse_iterator   rend()     const   { return array.rend();   }

// References

  reference              operator [] (size_type n)       { return array[n]; }
  const_reference        operator [] (size_type n) const { return array[n]; }

  reference              front()       { return array.front(); }
  const_reference        front() const { return array.front(); }
  reference              back()        { return array.back();  }
  const_reference        back()  const { return array.back();  }

  LimitedVector slice(size_type begin_index, size_type end_index) const
  {
    LimitedVector s(end_index-begin_index);
    iterator i = begin() + begin_index;
    iterator j = s.begin();
    for( ; i != begin()+end_index; ++i, ++j)
      *i = *j;
    return s;
  }

  LimitedVector slice(const_iterator begin_iter, const_iterator end_iter) const
  {
    LimitedVector s( (size_t) (end_iter-begin_iter) );
    iterator i = begin_iter;
    iterator j = s.begin();
    for( ; i != end_iter; ++i, ++j)
      *i = *j;
    return s;
  }

protected:

  std::vector<T> array;    
};

template <class T, class S>
bool operator==(const LimitedVector<T>& x, const LimitedVector<S>& y);

template <class T, class S>
bool operator<(const LimitedVector<T>& x, const LimitedVector<S>& y);

//
//  class Matrix1D - 1 dimensional vector class.  Includes a row/column
//      to show direction, otherwise very similar to LimitedVector.
//

template <class T>
class Matrix1D : public MatrixBase, public LimitedVector<T>
{
public:

  typedef LimitedVector<T> base;
  typedef __matrix_dir     direction;
  typedef typename MatrixDefines<T>::size_type        size_type;

// Constructors/Destructor

  Matrix1D(const direction d = Row_Order)
    : base()        { _dir = d; }
  Matrix1D(const size_type n, const T& t=T(), direction d = Row_Order)
    : base(n,t)     { _dir = d; }
  Matrix1D(const size_type n, const direction d)
    : base(n,T())   { _dir = d; }
  Matrix1D(const Matrix1D<T>& a)
    : base(a.array) { _dir = a._dir; clear(); setstate(a.rdstate()); }
  Matrix1D(const std::vector<T>& v, const direction d = Row_Order)
    : base(v)       { _dir = d; }

  ~Matrix1D() {};

// Input

  Matrix1D& operator= (const Matrix1D<T>& a)
    { _dir = a._dir; LimitedVector<T>::array = a.array; clear(); setstate(a.rdstate());
      return *this; }

  direction convert(direction d) { return _dir = d; }
  direction convert()            { if(_dir == Row_Order) 
                                     return convert(Col_Order);
                                   else return convert(Row_Order); }

  direction dir() const          { return _dir; }

protected:

  direction _dir; 

};

template <class T, class S>
inline bool operator==(const Matrix1D<T>& x, const Matrix1D<S>& y);

template <class T, class S>
inline bool operator<(const Matrix1D<T>& x, const Matrix1D<S>& y);

//
// class Matrix2D - 2 dimentional matrix class. 
//

template <class T, class MatType> class Matrix2D_Iterator;
template <class T>                class Const_Matrix2D_Iterator;

template <class T>
class Matrix2D : public MatrixBase
{
protected:
  typedef LimitedVector<T>      MatRow;
  typedef LimitedVector<MatRow> Matrix;

public:

  typedef typename MatrixDefines<T>::value_type        value_type;
  typedef typename MatrixDefines<T>::size_type         size_type;
  typedef typename MatrixDefines<T>::distance          distance;
  typedef typename MatrixDefines<T>::reference         reference;
  typedef typename MatrixDefines<T>::const_reference   const_reference;
  typedef          Matrix2D_Iterator<T, Matrix2D<T> >  iterator;
  typedef          Const_Matrix2D_Iterator<T>          const_iterator;
  typedef __matrix_dir                        it_dir;
  typedef __iter_type                         it_type;
  typedef std::pair<size_type, size_type>          mat_size;

// Constructors / Destructor

  Matrix2D() { clear(); row = col = 0; };
  Matrix2D(const size_type n, const size_type m, const T& t = T())
      { clear(); row = n; col = m; 
        LimitedVector<T> mat(m,t); matrix.resize(n,mat); }
  Matrix2D(const Matrix1D<T>&);
  Matrix2D(const Matrix2D<T>& m) { row = col = 0; (*this) = m; }

  Matrix2D& operator= (const Matrix2D<T>& m)
  { 
    if(&m != this)
    {
      resize(m.rows(), m.cols()); 
      const_iterator i = m.begin();
      iterator       j = begin();
      
      for(;i != m.end(); ++i, ++j)
        *j = *i;
         
      clear(); setstate(m.rdstate()); 
    }
    return *this;
  }
  
  bool operator<(const Matrix2D &) const { return false; }
  
  ~Matrix2D() {}

// Input

  void resize(const size_type n, const size_type m, const T& t=T())
      { if (n > rows()) { resizerow(m, t); resizecol(n, t); }
        else            { resizecol(n, t); resizerow(m, t); } }

  void resizerow(const size_type m, const T& t=T())
      { if (m != cols())
        for (typename Matrix::iterator i=matrix.begin(); i<matrix.end(); ++i)
          (*i).resize(m,t);
        col = m;
      }
  void resizecol(const size_type n, const T& t=T())
      { if(n != rows())
        { LimitedVector<T> mat(cols(),t);
          matrix.resize(n,mat); }
        row = n;
      }

  void set(const T& t= T() )
      { if (matrix.size() && matrix[0].size())
        { LimitedVector<T> b(col, t); matrix.set(b); } }

  void fill(const T& t = T() )
  {
     set(t);
  }

  void minimize()
      { if (matrix.capacity() > rows())
        {
          Matrix m2(rows(),MatRow(cols()));
          for (iterator i = begin(); i != end(); ++i)
            m2[i.rows()][i.cols()] = *i;
          matrix.swap(m2);
        }
      }
  void swap(Matrix2D<T>& m2)
      { size_type t = row; row = m2.row; m2.row = t;
                  t = col; col = m2.col; m2.col = t;
        matrix.swap(m2.matrix);
      }
  void insertrow(const size_type i, const T& t = T())
      { if (i <= rows())
        { resizecol(rows() + 1, t);
          for (size_type j = rows()-1; j > i; --j)
            matrix[j].swap(matrix[j-1]);
        }
      }
  void insertcol(const size_type i, const T& t = T())
      { if (i <= cols())
        { for (size_type j = 0; j < rows(); ++j)
            matrix[j].insert(i, t);
          col++;
        }
      }
  void eraserow(const size_type i)
      { if (i < rows())
        { for (size_type j = i+1; j < rows(); ++j)
            matrix[j].swap(matrix[j-1]);
          resizecol(rows() - 1);
        }
      }
  void erasecol(const size_type i)
      { if (i < cols())
        { for (size_type j = 0; j < rows(); ++j)
            matrix[j].erase(i);
          col--;
        }
      }

// Output

  mat_size  size()     const { return make_pair(rows(),cols());}
  size_type rows()     const { return row; }
  size_type cols()     const { return col; }
  mat_size  capacity() const 
    { if (rows())
        return make_pair(matrix.capacity(), matrix[0].capacity() );
      else
        return make_pair(matrix.capacity(), (size_type) 0 );         }

// Iterators

  iterator       begin(const size_type n = 0, const it_dir  d = Row_Order,
                       const it_type t = Forward_Iter)
      { 
        if( !rows() || !cols() )
          return iterator(*this, 0, 0);
      
        if ((rows() && cols()) || (t != Reverse_Iter))
        { bool b = (t == Reverse_Iter);
          return iterator(*this, d ? b*(rows()-1) : n, 
                                 d ? n : b*(cols()-1), d, t); }
        else return   iterator(*this, rows() - 1, cols() - 1, d, t);
      }

  const_iterator begin(const size_type n = 0, const it_dir  d = Row_Order, 
                       const it_type t = Forward_Iter) const
      { 
        if( !rows() || !cols() )
          return const_iterator(*this, 0, 0);
          
        if ((rows() && cols()) || (t != Reverse_Iter))
        { bool b = (t == Reverse_Iter);
          return const_iterator(*this, d ? b*(rows()-1) : n,
                                       d ? n : b*(cols()-1), d, t); }
        else return const_iterator(*this, rows() - 1, cols() - 1, d, t);
      }

  iterator       begin(const it_dir  d, const it_type t = Forward_Iter)
      { return (t==Forward_Iter) ? begin(0, d, t)
                                 : begin(!(d)*rows() + d*cols() - 1, d, t); }

  const_iterator begin(const it_dir  d, const it_type t = Forward_Iter) const
      { return (t==Forward_Iter) ? begin(0, d, t)
                                 : begin(!(d)*rows() + d*cols() - 1, d, t); }

  iterator       begin(const it_type t)       { return begin(Row_Order, t); }
  const_iterator begin(const it_type t) const { return begin(Row_Order, t); }
  iterator       begin(const size_type n,
                       const it_type t)       { return begin(n,Row_Order,t);}
  const_iterator begin(const size_type n,
                       const it_type t) const { return begin(n,Row_Order,t);}

  iterator       end()
  { 
    if( !rows() || !cols() )
      return iterator(*this, 0, 0);
    return iterator(*this, rows(), 0); 
  }
  const_iterator end() const 
  { 
    if( !rows() || !cols() )
      return const_iterator(*this, 0, 0);
    return const_iterator(*this, rows(), 0); 
  }

  iterator       end(const size_type n, const it_dir d=Row_Order, 
                     const it_type t=Forward_Iter)
      { return begin(n + (t==Forward_Iter) * 2 - 1, d, t); }
  const_iterator end(const size_type n, const it_dir d=Row_Order, 
                     const it_type t=Forward_Iter) const 
      { return begin(n + (t==Forward_Iter) * 2 - 1, d, t); }
  iterator       end(const it_dir d, const it_type t = Forward_Iter)
       { return begin(t * (!(d) * rows() + d * cols()) + !(t) * -1, d, t); }
  const_iterator end(const it_dir d, const it_type t = Forward_Iter) const
       { return begin(t * (!(d) * rows() + d * cols()) + !(t) * -1, d, t); }

  iterator       end(const it_type t)          { return end(Row_Order, t); }
  const_iterator end(const it_type t) const    { return end(Row_Order, t); }
  iterator       end(const size_type n, const it_type t)
                                               { return end(n,Row_Order,t);}
  const_iterator end(const size_type n, const it_type t) const 
                                               { return end(n,Row_Order,t);}

// Reverse Iterators (rbegin and rend to be compatible with STL)

  iterator       rbegin(const size_type n,
                        const it_dir    d = Row_Order)
      { return begin(n, d, Reverse_Iter); }
  const_iterator rbegin(const size_type n,
                        const it_dir    d = Row_Order) const
      { return begin(n, d, Reverse_Iter); }
  iterator       rbegin(const it_dir d = Row_Order)
      { return begin(d, Reverse_Iter); }
  const_iterator rbegin(const it_dir d = Row_Order) const
      { return begin(d, Reverse_Iter); }

  iterator         rend(const size_type n = 0,
                        const it_dir    d = Row_Order)
      { return end(n, d, Reverse_Iter); }
  const_iterator   rend(const size_type n = 0,
                        const it_dir    d = Row_Order) const
      { return end(n, d, Reverse_Iter); }
  iterator         rend(const it_dir d)
      { return end(d, Reverse_Iter); }
  const_iterator   rend(const it_dir d) const
      { return end(d, Reverse_Iter); }

// References

  reference       operator () (const size_type n, const size_type m) 
      { assert(n < rows() && m < cols()); return matrix[n][m]; }
  const_reference operator () (const size_type n, const size_type m) const
      { assert(n < rows() && m < cols()); return matrix[n][m]; }

  reference           front()       { return matrix.front().front(); }
  const_reference     front() const { return matrix.front().front(); }
  reference            back()       { return matrix.back().back(); }
  const_reference      back() const { return matrix.back().back(); }


  Matrix2D extract(size_type r, size_type c,
                   size_type top = 0, size_type  left = 0) const
  {
    Matrix2D m(r, c);
    for(size_t i = 0; i < r; ++i)
      for(size_t j = 0; j < c; ++j)
        m(i,j) = (*this)(i+top,j+left);
    return m;
  }

  void update(const Matrix2D& m, size_type top = 0, size_type left = 0)
  {
    size_type r = m.rows();
    size_type c = m.cols();

    for(size_t i = 0; i < r; ++i)
      for(size_t j = 0; j < c; ++j)
        (*this)(i+top,j+left) = m(i,j);
  }
    
protected:

  size_type row, col;
  Matrix   matrix;
};

template <class T, class S>
bool operator == (const Matrix2D<T>& m, const Matrix2D<S>& n);

//
// Matrix2D_Iterator - The Iterator classes for Matrix2D (iterators for
//      Matrix1D are the same as those for STL vectors).  Allows for
//      iteration by rows, columns, both forward and reverse, comparisons,
//      random access, etc.
//

template <class T, class MatType>
class Matrix2D_Iterator
{
public:

  typedef T&                            Reference;
  typedef T&                            reference;
  typedef T*                            pointer;
  typedef Matrix2D_Iterator<T, MatType> identity;
  typedef identity                      return_type;
  typedef __matrix_dir                  it_dir;
  typedef __iter_type                   it_type;
  typedef size_t                        size_type;
  typedef ptrdiff_t                     distance;
  typedef ptrdiff_t                     difference_type;
  typedef std::random_access_iterator_tag    iterator_category;
  typedef T                             value_type;
  
// Constructors/Destructor
  Matrix2D_Iterator() : matrix(NULL) {};
  Matrix2D_Iterator(const identity& b) : matrix(b.matrix), row(b.row),
    col(b.col), dir(b.dir), iter_type(b.iter_type) {};
  Matrix2D_Iterator(MatType& mat, size_type n, size_type m,
                    it_dir d = Row_Order, it_type i = Forward_Iter)
    : matrix(&mat), row(n), col(m), dir(d), iter_type(i)
    { if (iter_type == Reverse_Iter) reverse(); };
   ~Matrix2D_Iterator() {};

// Incrementation Functions

  return_type& cplus(distance s = 1)
    { s += col * matrix->rows() + row; col = s / matrix->rows(); 
      row = s % matrix->rows(); return *this; }
  return_type& rplus(distance s = 1)
    { s += row * matrix->cols() + col; row = s / matrix->cols();
      col = s % matrix->cols(); return *this; }
  return_type& cminus(distance s = 1)
    { s = col * matrix->rows() + row-s; col = s / matrix->rows();
      row = s % matrix->rows(); return *this; }
  return_type& rminus(distance s = 1)
    { s = row * matrix->cols() + col-s; row = s / matrix->cols();
      col = s % matrix->cols(); return *this; }

// Operators

  return_type& operator ++ ()
    { if (dir) { if (++row == matrix->rows()) { row = 0; col++; }}
      else     { if (++col == matrix->cols()) { col = 0; row++; }}
      return *this; }
  return_type operator ++ (int) { return_type i(*this); ++(*this); return i; }

  return_type& operator -- ()
    { if (dir) { if (row == 0) { row = matrix->rows() - 1; col--; } }
      else     { if (col == 0) { col = matrix->rows() - 1; row--; } }
      return *this; }
  return_type operator -- (int) { return_type i(*this); --(*this); return i; }

  return_type operator + (const size_type s) const
    { return_type i(*this); return i += s; }

  return_type& operator += (size_type s)
    { if (dir) { if ((row += s) >= matrix->rows()) return cplus(0); }
      else     { if ((col += s) >= matrix->cols()) return rplus(0); }
      return *this; }

  return_type  operator - (const size_type s) const
    { identity i(*this); return i -= s; }

  return_type& operator -= (size_type s)
    { if (dir) { if (row < s) return cminus(s); row -= s; }
      else     { if (col < s) return rminus(s); col -= s; }
      return *this; }

// Converters - change the direction of an iterator.

  return_type& convert(it_dir d, it_type t)
    { convert(d); return convert(t); }
  return_type& convert(it_dir d)
    { dir = d; return *this; }
  return_type& convert()
    { return convert(Row_Order, Forward_Iter); }
  return_type& convert(it_type t)
    { if (iter_type==t) return *this; iter_type=t; reverse(); return *this; }

// Output

  size_type   rows() const
    { return (iter_type) ? row : matrix->rows() - row - 1; }
  size_type   cols() const
    { return (iter_type) ? col : matrix->cols() - col - 1; }
  MatType&     mat() const { return *matrix; }
  it_dir  curr_dir() const { return dir; }
  it_type     type() const { return iter_type; }

  Reference operator *  () const { return (*matrix)(rows(), cols()); }

  Reference operator [] (const distance s) const{ return *(*this + s); }

protected:

  friend class Const_Matrix2D_Iterator<T>;

  MatType   *matrix;
  size_type row, col;
  it_dir    dir;
  it_type   iter_type;

// method function

  void reverse()
    { row = matrix->rows() - row - 1; col = matrix->cols() - col - 1; }
};

template <class T>
class Const_Matrix2D_Iterator : 
    public Matrix2D_Iterator<const T, const Matrix2D<T> >
{
public:
  typedef Matrix2D<T>                              MatType;
  typedef Matrix2D_Iterator<const T, const Matrix2D<T> > base;
  typedef Matrix2D_Iterator<T, Matrix2D<T> >       iterator;
  typedef typename MatrixDefines<T>::size_type     size_type;
  typedef __matrix_dir                             it_dir;  
  typedef __iter_type                              it_type;     
  
  Const_Matrix2D_Iterator() :                            base()  {}
//Const_Matrix2D_Iterator(const base& b)               : base(b) {}
  Const_Matrix2D_Iterator(const iterator& b) 
    : base(*b.matrix, b.row, b.col, b.dir, b.iter_type)          {}
  Const_Matrix2D_Iterator(const MatType& mat, size_type n, size_type m,
                        it_dir d = Row_Order, it_type i = Forward_Iter)
    : base( (MatType &)mat, n, m, d, i)                                      {}
    
  Const_Matrix2D_Iterator& operator = (const iterator& b)
  { LimitedVector<T>::matrix = b.matrix; 
    LimitedVector<T>::row = b.row;
    LimitedVector<T>::col = b.col;
    LimitedVector<T>::dir = b.dir;
    LimitedVector<T>::iter_type = b.iter_type; return *this; }
};

template <class T, class MatType>
inline std::bidirectional_iterator_tag
iterator_category(const Matrix2D_Iterator<T,MatType>&) {
  return std::bidirectional_iterator_tag();
}

template <class T, class MatType>
inline T*
value_type(const Matrix2D_Iterator<T,MatType>&) {
  return (T*) 0;
}

template <class T, class MatType>
inline ptrdiff_t*
distance_type(const Matrix2D_Iterator<T,MatType>&) {
  return (ptrdiff_t*) 0;
}

template <class T>
inline std::bidirectional_iterator_tag
iterator_category(const Const_Matrix2D_Iterator<T> &) {
  return std::bidirectional_iterator_tag();
}

template <class T>
inline T*
value_type(Const_Matrix2D_Iterator<T>&) {
  return (T*) 0;
}

template <class T>
inline ptrdiff_t*
distance_type(const Const_Matrix2D_Iterator<T>&) {
  return (ptrdiff_t*) 0;
}

//
//  iterator comparison functions.  < == and - functions for Base Iterators.
//

template <class T1, class M1, class T2, class M2>
bool operator < (const Matrix2D_Iterator<T1, M1>& x, 
                 const Matrix2D_Iterator<T2, M2>& y);

template <class T, class R, class M>
bool operator <  (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y);

template <class T, class R, class M>
bool operator <  (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T>
bool operator <  (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T1, class M1, class T2, class M2>
bool operator == (const Matrix2D_Iterator<T1, M1>& x, const Matrix2D_Iterator<T2, M2>& y);

template <class T, class R, class M>
bool operator == (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y);

template <class T, class R, class M>
bool operator == (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T>
bool operator == (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T1, class M1, class T2, class M2>
bool operator != (const Matrix2D_Iterator<T1, M1>& x, const Matrix2D_Iterator<T2, M2>& y);

template <class T, class R, class M>
bool operator != (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y);

template <class T, class R, class M>
bool operator != (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T>
bool operator != (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T1, class M1, class T2, class M2>
bool operator - (Matrix2D_Iterator<T1, M1> x, Matrix2D_Iterator<T2, M2> y);

template <class T, class R, class M>
bool operator - (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y);

template <class T, class R, class M>
bool operator - (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y);

template <class T>
bool operator - (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y);

#include "matrix.ipp"

#endif
