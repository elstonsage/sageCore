////////////////////////////////////////////////////////////////////////////
//             Implementation of lodpal_result(Inline)                    //
////////////////////////////////////////////////////////////////////////////

inline
void
lodpal_result::clear()
{
  my_method                = 'c';
  my_lod_score_mm          = std::numeric_limits<double>::quiet_NaN();
  my_lod_score_mf          = std::numeric_limits<double>::quiet_NaN();
  my_lod_score_ff          = std::numeric_limits<double>::quiet_NaN();
}

inline
void
lodpal_result::set_maxfun_result(const maxfun_results& re)
{
  my_maxfun_result = re;

  my_vc_matrix.resize(0, 0);
  my_vc_matrix.resize(re.getCovarianceMatrix().getSize(), re.getCovarianceMatrix().getSize());

  for( int t1 = 0; t1 < re.getCovarianceMatrix().getSize(); ++t1 )
    for( int t2 = 0; t2 < re.getCovarianceMatrix().getSize(); ++t2 )
      my_vc_matrix(t1, t2) = re.getCovarianceMatrix().getCovariance(t1, t2);
}

inline
void
lodpal_result::set_method(char c)
{
  my_method = c;
}

inline
void
lodpal_result::set_lod_score(const pair<double, double>& ld)
{
  my_lod_score_with_cap    = ld.first;
  my_lod_score_without_cap = ld.second;
}

inline
void
lodpal_result::set_lod_score_mm(double mm)
{
  my_lod_score_mm = mm;
}

inline
void
lodpal_result::set_lod_score_mf(double mf)
{
  my_lod_score_mf = mf;
}

inline
void
lodpal_result::set_lod_score_ff(double ff)
{
  my_lod_score_ff = ff;
}

inline
double
lodpal_result::lod_score() const
{
//  return my_maxfun_result.getFinalFunctionValue();
  return my_lod_score_without_cap;
}

inline
int
lodpal_result::function_evaluations() const
{
  return my_maxfun_result.getIterations();
}

inline
int
lodpal_result::last_error() const
{
  return my_maxfun_result.getExitFlag();
}

inline
int
lodpal_result::ivfl() const
{
  return my_maxfun_result.getCovMatrixStatus();
}

inline
int
lodpal_result::igage() const
{
  return my_maxfun_result.getAgeOfGradientVector();
}

inline
char
lodpal_result::method() const
{
  return my_method;
}

inline
const Matrix2D<double>&
lodpal_result::var_cov_matrix() const
{
  return my_vc_matrix;
}

inline
double
lodpal_result::lod_score_mm() const
{
  return my_lod_score_mm;
}

inline
double
lodpal_result::lod_score_mf() const
{
  return my_lod_score_mf;
}

inline
double
lodpal_result::lod_score_ff() const
{
  return my_lod_score_ff;
}
