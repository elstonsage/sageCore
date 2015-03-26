#include "relpal/two_level_calculator.h"

namespace SAGE   {
namespace RELPAL {

calculator_base::calculator_base(size_t m)
               : parameters(m)
{}

calculator_base::~calculator_base()
{}

void
calculator_base::touch_up_matrix(matrix& m) const
{
  for( size_t i = 0; i < m.rows(); ++i )
    for( size_t j = 0; j < m.cols(); ++j )
      if( fabs(m(i, j)) < 1.0e-10 )
        m(i, j) = 0.0;
                    
  return;
}

int
calculator_base::get_tri_inverse(const trimatrix& w, trimatrix& wi) const
{
  assert( w.size() );

  wi = w;

  matrix ww;
  ww.resize(w.size(), w.size());
  for( size_t r = 0; r < w.size(); ++r )
    for( size_t c = 0; c < w.size(); ++c )
      ww(r,c) = w(r,c);

  matrix wwi;

  assert( ww && ww.rows() &&ww.cols() );

  SAGE::SVD svd;
  svd.compute(ww);
  svd.inverse(wwi);

  int new_info = 0;

  if( !wwi )
  {
    new_info = 1;
  }

  if( new_info != 1 )
  {
    for( size_t r = 0; r < ww.rows(); ++r )
      for( size_t c = 0; c < ww.cols(); ++c )
        wi(r,c) = wwi(r,c);
  }

  return new_info;
}

void
calculator_base::get_tri_kron(const trimatrix& b, trimatrix& wi) const
{
  assert( b.size() );

  wi.resize(b.size()*b.size());

  for( size_t r1 = 0, r = 0; r1 < b.size(); ++r1 )
  {
    for( size_t r2 = 0; r2 < b.size(); ++r2, ++r )
    {
      for( size_t c1 = 0, c = 0; c1 < b.size(); ++c1 )
      {
        for( size_t c2 = 0; c2 < b.size(); ++c2, ++c )
          wi(r,c) = 0.5 * b(r1,c1) * b(r2,c2);
      }
    }
  }

  return;
}

matrix&
calculator_base::XtriZ(const matrix& X, const trimatrix& Z, matrix& o) const
{
  if( !X || X.cols() != Z.size() )
  {
    o.setstate(FortranMatrix<double>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(X.rows(), Z.size());

  for( size_t i = 0; i < o.rows(); ++i )
    for( size_t j = 0; j < o.cols(); ++j )
    {
      double dot = 0.0;
      for( size_t k = 0; k < X.cols(); ++k )
        dot += X(i,k) * Z(k,j);
      o(i,j) = dot;
    }

  return o;
}

matrix&
calculator_base::triXZ(const trimatrix& X, const matrix& Z, matrix& o) const
{
  if( !Z || X.size() != Z.rows() )
  {
    o.setstate(FortranMatrix<double>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(X.size(), Z.cols());

  for( size_t i = 0; i < o.rows(); ++i )
    for( size_t j = 0; j < o.cols(); ++j )
    {
      double dot = 0.0;
      for( size_t k = 0; k < X.size(); ++k )
        dot += X(i,k) * Z(k,j);
      o(i,j) = dot;
    }

  return o;
}

matrix&
calculator_base::XTtriZ(const matrix& X, const trimatrix& Z, matrix& o) const
{
  if( !X || X.rows() != Z.size() )
  {
    o.setstate(FortranMatrix<double>::failbit);
    return o;
  }

  o.clear();
  o.resize_nofill(X.cols(), Z.size());

  for( size_t i = 0; i < o.rows(); ++i )
    for( size_t j = 0; j < o.cols(); ++j )
    {
      double dot = 0.0;
      for( size_t k = 0; k < X.rows(); ++k )
        dot += X(k,i) * Z(k,j);
      o(i,j) = dot;
    }

  return o;
}

} // end of namespace RELPAL
} // end of namespace SAGE
