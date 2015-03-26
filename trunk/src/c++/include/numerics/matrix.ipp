
inline
MatrixBase::MatrixBase()
{
    _state = 0;
}

inline
MatrixBase::~MatrixBase()
{
}

inline
void MatrixBase::clear()
{
  _state = 0;
}

inline
void MatrixBase::clear(const int state)
{
  _state &= ~state;
}

inline
void MatrixBase::setstate(const int flag)
{
  _state |= flag;
}

inline
int MatrixBase::good() const
{
  return _state == 0;
}

inline
int MatrixBase::fail() const
{
  return _state & (badbit | failbit);
}

inline
int MatrixBase::bad() const
{
  return _state & badbit;
}

inline
int MatrixBase::rdstate() const
{
  return _state;
}

inline
MatrixBase::operator void*() const
{
  return fail() ? (void *) 0 : (void *) (-1);
}

inline
int MatrixBase::operator!() const
{
  return fail();
}

template <class T, class S>
inline bool operator==(const LimitedVector<T>& x, const LimitedVector<S>& y)
{ return x.size() == y.size() && equal(x.begin(), x.end(), y.begin()); }

template <class T, class S>
inline bool operator<(const LimitedVector<T>& x, const LimitedVector<S>& y)
{ return lexicographical_compare(x.begin(), x.end(), y.begin(), y.end()); }

template <class T, class S>
inline bool operator==(const Matrix1D<T>& x, const Matrix1D<S>& y)
{ return x.size() == y.size()                 &&
         x.dir()  == y.dir()                  &&
         equal(x.begin(), x.end(), y.begin());    }

template <class T, class S>
inline bool operator<(const Matrix1D<T>& x, const Matrix1D<S>& y)
{ return lexicographical_compare(x.begin(), x.end(), y.begin(), y.end()); }

template <class T>
Matrix2D<T>::Matrix2D(const Matrix1D<T>& vect) : row(0), col(0)
{
  if (vect.dir() == Row_Order) matrix.resize(1,vect);
  else
  {
    resize(vect.size(),1);
    std::copy(vect.begin(), vect.end(), begin(Col_Order));
  }
  clear();
  setstate(vect.rdstate());
  row    = matrix.size();
  col    = matrix.front().size();
}

template <class T, class S>
inline bool operator == (const Matrix2D<T>& m, const Matrix2D<S>& n)
{
  return m.rdstate() == n.rdstate()            &&
         m.rows()    == n.rows()               &&
         m.cols()    == n.cols()               && 
         equal(m.begin(), m.end(), n.begin());
}

template <class T1, class M1, class T2, class M2>
inline bool operator < (const Matrix2D_Iterator<T1, M1>& x, 
                        const Matrix2D_Iterator<T2, M2>& y)
{
  MatrixDefines<int>::distance xr(x.rows()), xc(x.cols()),
                               yr(y.rows()), yc(y.cols());

  if (x.type())
  {
    if (x.curr_dir())
      return ((xc < yc) || (xc == yc && xr < yr));
    else
      return ((xr < yr) || (xr == yr && xc < yc));
  }
  if (x.curr_dir())
      return ((xc > yc) || (xc == yc && xr > yr));
    else
      return ((xr > yr) || (xr == yr && xc > yc));
}

template <class T, class R, class M>
inline bool operator <  (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x <  y; }

template <class T, class R, class M>
inline bool operator <  (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y)
{ return x <  (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T>
inline bool operator <  (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x <  
         (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T1, class M1, class T2, class M2>
inline bool operator == (const Matrix2D_Iterator<T1, M1>& x, const Matrix2D_Iterator<T2, M2>& y)
{ return &x.mat() == &y.mat() && 
         x.rows() == y.rows() &&
         x.cols() == y.cols(); }

template <class T, class R, class M>
inline bool operator == (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x == y; }

template <class T, class R, class M>
inline bool operator == (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y)
{ return x == (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T>
inline bool operator == (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x == 
         (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T1, class M1, class T2, class M2>
inline bool operator != (const Matrix2D_Iterator<T1, M1>& x, const Matrix2D_Iterator<T2, M2>& y)
{ return &x.mat() != &y.mat() || 
         x.rows() != y.rows() ||
         x.cols() != y.cols(); }

template <class T, class R, class M>
inline bool operator != (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x != y; }

template <class T, class R, class M>
inline bool operator != (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y)
{ return x != (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T>
inline bool operator != (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x != 
         (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T1, class M1, class T2, class M2>
inline bool operator - (Matrix2D_Iterator<T1, M1> x, Matrix2D_Iterator<T2, M2> y)
{ 
  typename Matrix2D<T1>::distance i;
  if (x.curr_dir())
    i = (x.cols() - y.cols()) * x.mat().rows() + x.rows() - y.rows();
  else
    i = (x.rows() - y.rows()) * x.mat().rows() + x.cols() - y.cols();
  return (x.type()) ? i : !i;
}

template <class T, class R, class M>
inline bool operator - (const Const_Matrix2D_Iterator<T>& x, const Matrix2D_Iterator<R,M>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x - y; }

template <class T, class R, class M>
inline bool operator - (const Matrix2D_Iterator<R,M>& x, const Const_Matrix2D_Iterator<T>& y)
{ return x - (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }

template <class T>
inline bool operator - (const Const_Matrix2D_Iterator<T>& x, const Const_Matrix2D_Iterator<T>& y)
{ return (Matrix2D_Iterator<const T, const Matrix2D<T> >) x - 
         (Matrix2D_Iterator<const T, const Matrix2D<T> >) y; }
