#include "relpal/two_level_calculator.h"

namespace SAGE   {
namespace RELPAL {

relpal_least_square::relpal_least_square(size_t m)
                   : calculator_base(m)
{
  reset(m);

  beta.resize_fill(parameters, 1, QNAN);
  information.resize_fill(parameters, parameters, QNAN);
  Variance.resize_fill(parameters, parameters, QNAN);
  total_variance.resize_fill(1, 1, QNAN);
  residual_variance.resize_fill(1, 1, QNAN);
}

relpal_least_square::~relpal_least_square()
{}

void 
relpal_least_square::reset(size_t m)
{
  parameters = m;

  observation_count = 0;
  cluster_count     = 0;
  svd_return_code   = 1;

  size_t n = 1;

  AWA.clear();
  AWy.clear();
  yWy.clear();

  AWA.resize_fill(m, m, 0.0);
  AWy.resize_fill(m, n, 0.0);
  yWy.resize_fill(n, n, 0.0);
}

void 
relpal_least_square::add_block(const matrix& y, const trimatrix& W, const matrix& A, trimatrix& Wi)
{
  if(!AWA || !AWy || !yWy)
    return;

  assert( W.size() && W.size() == A.rows() && A.rows() == y.rows() );
  assert( A && A.cols() == parameters );
  assert( y && y.rows() > 0 );

  BA.clear();
  By.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();

  int info = get_tri_inverse(W, Wi);

  if( info != 0 )
  {
    AWA.setstate(matrix::badbit);
    AWy.setstate(matrix::badbit);
    yWy.setstate(matrix::badbit);
    Wi.resize(0);
    return;
  }

  observation_count += y.rows();
  cluster_count     += 1;

  XTtriZ(A, Wi, BA);
  XTtriZ(y, Wi, By);

  multiply(BA, A, temp1);
  multiply(BA, y, temp2);
  multiply(By, y, temp3);

  AWA += temp1;
  AWy += temp2;
  yWy += temp3;

#if 0
  print_trimatrix(Wi, cout, "V^-1");
  print_matrix(AWA, cout, "AWA");
  print_matrix(AWy, cout, "AWy");
  print_matrix(yWy, cout, "yWy");
#endif

  return;
}

void 
relpal_least_square::add_block_kron(const matrix& y, const trimatrix& b, const matrix& A)
{
  if(!AWA || !AWy || !yWy)
    return;

  assert( A && A.cols() == parameters );
  assert( y && y.rows() > 0 );

  BA.clear();
  By.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();

  trimatrix Wi;
  get_tri_kron(b, Wi);

  observation_count += y.rows();
  cluster_count     += 1;

  assert( Wi.size() && Wi.size() == A.rows() && A.rows() == y.rows() );

  XTtriZ(A, Wi, BA);
  XTtriZ(y, Wi, By);

  multiply(BA, A, temp1);
  multiply(BA, y, temp2);
  multiply(By, y, temp3);

  AWA += temp1;
  AWy += temp2;
  yWy += temp3;

#if 0
  print_trimatrix(Wi, cout, "(V*)^1");
  print_matrix(AWA, cout, "AWA");
  print_matrix(AWy, cout, "AWy");
  print_matrix(yWy, cout, "yWy");
#endif

  return;
}

bool 
relpal_least_square::compute()
{
  if( observation_count <= 0 || !AWA || !AWy )
  {
    beta.setstate(matrix::badbit);
    return false;
  }

  Variance.clear();
  information.clear();

  svd.inverse_of(AWA, information);

  if( !information )
  {
    reset();
    information.setstate(matrix::failbit);
    Variance.setstate(matrix::failbit);
    return false;
  }

#if 0
  print_matrix(information, cout, "information");
#endif

  multiply(information, AWy, beta);

  XTZ(AWy, beta, temp2);

  total_variance = yWy / ((double)observation_count);

  if( observation_count - parameters <= 0 )
    return false;

  residual_variance  = ( yWy - temp2 ) / (double)(observation_count-parameters);

  assert(    (information.rows()       == 1 && information.cols()       == 1)
          || (residual_variance.rows() == 1 && residual_variance.cols() == 1) );

  Variance = kron(information, residual_variance);

  return true;
}

void
relpal_least_square::build_residuals(const matrix& y, const matrix& A, matrix& r) const
{
  r = y;
  r -= A*beta;
}  

} // end of namespace RELPAL
} // end of namespace SAGE
