#include "sibpal/gls3.h"

#define CHOLESKY_INVERSE 0

namespace SAGE   {
namespace SIBPAL {

UnivariateGeneralizedLeastSquares3::
UnivariateGeneralizedLeastSquares3(size_t m)
{
  set_inverse_method(SVD);
  set_leverage_adjustment(false);
  reset(m);
}

UnivariateGeneralizedLeastSquares3::~UnivariateGeneralizedLeastSquares3()
{}


size_t
UnivariateGeneralizedLeastSquares3::get_svd_return_code() const
{
  return svd_return_code;
}

void
UnivariateGeneralizedLeastSquares3::set_inverse_method(inverse_type d)
{
  my_inverse_method = d;
}

void
UnivariateGeneralizedLeastSquares3::set_inverse_method(const string& l)
{
  string method = toUpper(l);

  if     ( method == "SVD" )
    set_inverse_method(SVD);
  else if( method == "LU" || method == "LUP")
    set_inverse_method(LUP);
  //else
  //  errors << priority(error) << "Unknown inverse method specified '" << l << "'." << endl;
}

UnivariateGeneralizedLeastSquares3::inverse_type
UnivariateGeneralizedLeastSquares3::inverse_method() const
{
  return my_inverse_method;
}

void
UnivariateGeneralizedLeastSquares3::set_leverage_adjustment(bool adjust)
{
  my_leverage_adjustment = adjust;
}

bool
UnivariateGeneralizedLeastSquares3::leverage_adjustment() const
{
  return my_leverage_adjustment;
}

void
UnivariateGeneralizedLeastSquares3::reset_results()
{
  const double qNaN = std::numeric_limits<double>::quiet_NaN();

  beta.resize_fill(parameters, variates, qNaN);

  information.resize_fill(parameters, parameters, std::numeric_limits<double>::quiet_NaN());
  Variance.resize_fill(parameters, parameters, std::numeric_limits<double>::quiet_NaN());
  SandwichVariance.resize_fill(parameters, parameters, std::numeric_limits<double>::quiet_NaN());
  total_variance.resize_fill(variates, variates, qNaN);
  residual_variance.resize_fill(variates, variates, qNaN);

  beta.setstate(matrix::failbit);
  information.setstate(matrix::failbit);
  Variance.setstate(matrix::failbit);
  SandwichVariance.setstate(matrix::failbit);
  total_variance.setstate(matrix::failbit);
  residual_variance.setstate(matrix::failbit);

  observation_count = 0;
  cluster_count     = 0;
  svd_return_code   = 1;
}  

void 
UnivariateGeneralizedLeastSquares3::reset()
{
  reset(parameters);
}

void 
UnivariateGeneralizedLeastSquares3::reset(size_t m, size_t n)
{
  parameters = m;
  variates   = n;

  reset_results();

  yWy.resize_fill(n,n, 0.0);

  AWA.clear();
  AWy.clear();
  OPG.clear();

  AWA.resize_fill(m, m, 0.0);
  AWy.resize_fill(m, n, 0.0);
  OPG.resize_fill(m, m, 0.0);

  return;
}

// A way to avoid different floating point operation to close-to-0 value on different system
//
void
touch_up_matrix(matrix& m)
{
  for( size_t i = 0; i < m.rows(); ++i )
    for( size_t j = 0; j < m.cols(); ++j )
      if( fabs(m(i, j)) < 1.0e-10 )
        m(i, j) = 0.0;

  return;
}

void 
UnivariateGeneralizedLeastSquares3::add_block(matrix& y, const matrix& W, const matrix& A, 
                                              const weight_status_type& status)
{
#if 0
  print_matrix(y, cout, "y");
  print_matrix(W, cout, "W");
  print_matrix(A, cout, "A");
#endif

  if(!AWA || !AWy || !yWy)
    return;

  // These assertions will only trigger for programming errors, not
  // data errors.  Thus, they can stay unchanged.
  assert( W && W.rows() && W.rows() == W.cols() && A.rows() == W.rows() && A.rows() == y.rows() );
  assert( A && A.cols() == parameters );
  assert( y && y.rows() > 0 );
  assert( status == INVERSEW || status == NORMALW );

  BA.clear();
  By.clear();
  B.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();

  observation_count += y.rows();
  cluster_count     += 1;

#if CHOLESKY_INVERSE
  Cholesky(W, B, false);
  if( status == NORMALW )
    TRIinverse(B, B, false);
#else

  svd_return_code = svd.compute(W);
  //cout << "det(W)  = " << svd.determinant()    << endl;
  //cout << "cond(W) = " << svd.well_condition() << endl;

  if( status == NORMALW )
    svd.inverse_sqrt(B);
  else
    svd.sqrt(B);

#endif

  if(!B)
  {
    AWA.setstate(matrix::badbit);
    AWy.setstate(matrix::badbit);
    yWy.setstate(matrix::badbit);
    return;
  }

  multiply(B, A, BA);
  multiply(B, y, By);

  XTX(BA,     temp1);
  XTZ(BA, By, temp2);
  XTX(By,     temp3);

  AWA += temp1;
  AWy += temp2;
  yWy += temp3;

#if 0
  print_matrix_first10(B*B, cout, "B*B");
  print_matrix(AWA, cout, "AWA");
  print_matrix(AWy, cout, "AWy");
  print_matrix(yWy, cout, "yWy");
#endif

  return;
}

void 
UnivariateGeneralizedLeastSquares3::add_block(matrix& y, const matrix& A)
{
#if 0
  print_matrix(y, cout, "y");
  print_matrix(A, cout, "A");
#endif

  if(!AWA || !AWy || !yWy)
    return;

  // These assertions will only trigger for programming errors, not
  // data errors.  Thus, they can stay unchanged.
  assert( A && A.rows() == y.rows() );
  assert( A && A.cols() == parameters );
  assert( y && y.rows() > 0 && y.cols() == variates);

  temp1.clear();
  temp2.clear();
  temp3.clear();

  observation_count += y.rows();
  cluster_count     += 1;

  XTX(A,     temp1);
  XTZ(A,  y, temp2);
  XTX(y,     temp3);

  AWA += temp1;
  AWy += temp2;
  yWy += temp3;

  return;
}

void 
UnivariateGeneralizedLeastSquares3::compute_covariance(const matrix& X, matrix& C)
{
  size_t n = X.rows();

  if(!X)
  {
    C.resize(n, n, std::numeric_limits<double>::quiet_NaN());
    C.setstate(matrix::failbit);
    return;
  }

  matrix R;

  Cholesky(X, R);

  C.clear();
  C.resize(n, n);

  // See Bjorck (1996) section 2.8.3 for details of this algorithm

  // Compute lower right corner
  C(n-1,n-1) = 1/(R(n-1,n-1)*R(n-1,n-1));

  // Compute right and bottom edges
  for(int i = n - 2; i >= 0; --i)
  {
    double sum = 0.0;
    for(size_t j = i + 1; j < n; ++j)
      sum += R(i,j)*C(j,n-1);
    C(n-1,i) = C(i,n-1) = -sum/R(i,i);
  }

  // Compute rows and columns from right to left
  for(int k = n - 2; k >= 0; --k)
  {
    // Compute diagonal element
    double sum = 0.0;
    for(size_t j = k + 1; j < n; ++j)
      sum += R(k,j)*C(k,j);
    C(k,k) = (1.0/R(k,k) - sum)/R(k,k);

    // Compute rest of row and column
    for(int i = k - 1; i >= 0; --i)
    {
      sum = 0.0;
      for(size_t j = i + 1; j <= (size_t)k; ++j)
        sum += R(i,j)*C(j,k);
      for(size_t j = k + 1; j < (size_t)n; ++j)
        sum += R(i,j)*C(k,j);
      C(k,i) = C(i,k) = -sum/R(i,i);
    }
  }

  return;
}

bool 
UnivariateGeneralizedLeastSquares3::compute()
{
  if(observation_count <= 0 || !AWA || !AWy)
  {
    beta.setstate(matrix::badbit);
    return false;
  }

  Variance.clear();
  information.clear();

#if 1
  if( inverse_method() == SVD)
  {
    svd.inverse_of(AWA, information);
  //cout << "det(AWA)  = " << svd.determinant()    << endl;
  //cout << "cond(AWA) = " << svd.well_condition() << endl;
  }
  else
    lup.inverse_of(AWA, information);
#else
  compute_covariance(AWA, information);
#endif

  if(!information)
  {
    reset_results();
    information.setstate(matrix::failbit);
    Variance.setstate(matrix::failbit);
    return false;
  }

  multiply(information, AWy, beta);

  XTZ(AWy, beta, temp2);

  total_variance     = yWy / ((double)observation_count);

  if(observation_count - parameters <= 0)
    return false;

  residual_variance  = ( yWy - temp2 ) / (double)(observation_count-parameters);

  // NOTICE: Check the order of the Kronecker product arguments!  This has
  //         only been tested for cases where order does not matter (i.e.
  //         one or the other argument is a 1x1 matrix.  Thankfully, this is
  //         true for all of the current analyses, but may not always be.

  assert(    (information.rows()       == 1 && information.cols()       == 1)
          || (residual_variance.rows() == 1 && residual_variance.cols() == 1) );

  Variance = kron(information, residual_variance);

  return true;
}

void 
UnivariateGeneralizedLeastSquares3::update_robust_variance_block(matrix& y, const matrix& A)
{
  // These assertions will only trigger for programming errors, not
  // data errors.  Thus, they can stay unchanged.
  assert( A.rows() == y.rows() );
  assert( A && A.cols() == parameters );

  BA.clear();
  By.clear();
  B.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();

  matrix& residuals = temp1;

  // Compute normal residuals
  residuals = y - A*beta;

  if( leverage_adjustment() )
  {
    // Compute leverage adjustment
    hat = A*information*transpose(A);
    hat = eye<double>(hat.rows()) - hat;
    svd.compute(hat);
    svd.inverse_sqrt(hat);

    // Compute leverage adjusted residuals
    residuals = hat * residuals;
  }

  // Compute outer product gradient
  XTZ(A, residuals, temp2);
  XXT(temp2,temp3);
  OPG += temp3;

  return;
}

void 
UnivariateGeneralizedLeastSquares3::update_robust_variance_block(matrix& y, const matrix& W, const matrix& A, 
                                                                 const weight_status_type& status)
{
  if(!OPG)
    return;

  // These assertions will only trigger for programming errors, not
  // data errors.  Thus, they can stay unchanged.
  assert( W && W.rows() && W.rows() == W.cols() && A.rows() == W.rows() && A.rows() == y.rows() );
  assert( A && A.cols() == parameters );
  assert( y && y.rows() > 0);
  assert( status == INVERSEW || status == NORMALW );

  BA.clear();
  By.clear();
  B.clear();
  temp1.clear();
  temp2.clear();
  temp3.clear();

#if CHOLESKY_INVERSE
  Cholesky(W, B, false);
  if( status == NORMALW )
    TRIinverse(B, B, false);
#else

  svd.compute(W);

  if( status == NORMALW )
    svd.inverse_sqrt(B);
  else
    svd.sqrt(B);

#endif

  if(!B)
  {
    OPG.setstate(matrix::badbit);
    return;
  }
  
  /*
    Robust variance estimator (Sandwich estimator)

    V = working covariance matrix (block diagonal)

            nxn            nxn          nxn
    W = (A' V^-1 A)^-1 [A' V^-1 (y - A beta)]* (A' V^-1 A)^-1

    where [X]* = X X'.

    Leverage adjusted residuals can also be used.

  */

  matrix& residuals = temp1;

  // Compute normal residuals
  residuals = y - A*beta;

  if( leverage_adjustment() )
  {
    // Compute leverage adjustment
    hat = A*information*transpose(A)*transpose(B)*B;
    hat = eye<double>(hat.rows()) - hat;
    svd.compute(hat);
    svd.inverse_sqrt(hat);

    // Compute leverage adjusted residuals
    residuals = hat * residuals;
  }

  // Compute outer product gradient
  multiply(B, A,         BA);
  multiply(B, residuals, By);
  XTZ(BA, By, temp2);
  XXT(temp2,temp3);
  OPG += temp3;

#if 0
  cout << "By =" << endl;
  print_matrix(By, cout);
  cout << "BA =" << endl;
  print_matrix(BA, cout);
  cout << "temp2=" << endl;
  print_matrix(temp2, cout);
  cout << "temp3=" << endl;
  print_matrix(temp3, cout);
  cout << "OPG=" << endl;
  print_matrix(OPG, cout);
#endif

  return;
}

bool 
UnivariateGeneralizedLeastSquares3::compute_robust_variance()
{
  if(observation_count - parameters <= 0)
  {
    //errors << priority(error) << "Too few valid pairs." << std::endl;
    beta.setstate(matrix::badbit);
    return false;
  }

  SandwichVariance.clear();

  if(!information)
  {
    reset_results();
    SandwichVariance.setstate(matrix::failbit);
    return false;
  }

  SandwichVariance  = information * OPG * information;

#if 0
  cout << "Residual variance = " << residual_variance(0,0) << endl;
  svd.compute(OPG);
  print_matrix(Variance, cout, "Variance");
  cout << "det(OPG) = "    << svd.determinant() << endl;
  cout << "1/cond(OPG) = " << svd.well_condition() << endl;
  print_matrix(SVDinverse(OPG), cout, "inv(OPG)");
  print_matrix(SandwichVariance, cout, "Sandwich");
#endif

  return true;
}

void
UnivariateGeneralizedLeastSquares3::build_residuals(const matrix& y, const matrix& A, matrix& r) const
{
  r = y;
  r -= A*beta;
}  

} // end of namespace SIBPAL
} // end of namespace SAGE
