////////////////////////////////////////////////////////////////////////////
//             Implementation of lodpal_analysis(Inline)                 //
////////////////////////////////////////////////////////////////////////////

inline
void
lodpal_analysis::set_lod_score(double l)
{
  my_lod_score = l;
}

inline
void
lodpal_analysis::set_function_evaluations(int nfe)
{
  my_function_evaluations = nfe;
}

inline
void
lodpal_analysis::set_last_error(int e)
{
  my_last_error = e;
}

inline
void
lodpal_analysis::set_ivfl(int i)
{
  my_ivfl = i;
}

inline
void
lodpal_analysis::set_igage(int i)
{
  my_igage = i;
}

inline
void
lodpal_analysis::set_method(char c)
{
  my_method = c;
}

inline
void
lodpal_analysis::set_var_cov_matrix(const Matrix2D<double>& m)
{
  my_vc_matrix.resize(m.rows(), m.cols());
  for( size_t i = 0; i < m.rows(); ++i )
    for( size_t j = 0; j < m.cols(); ++j )
      my_vc_matrix(i, j) = m(i, j);
}

inline
void
lodpal_analysis::set_parameters(const lodpal_parameters& p)
{
  invalidate();
  my_parameters = p;
}

inline
void
lodpal_analysis::set_previous_parameters(const lodpal_parameters& pp)
{
  my_previous_parameters = pp;
}

inline
lodpal_parameters&
lodpal_analysis::parameters()
{
  return my_parameters;
}

inline
const lodpal_parameters&
lodpal_analysis::parameters() const
{
  return my_parameters;
}

inline
lodpal_parameters&
lodpal_analysis::previous_parameters()
{
  return my_previous_parameters;
}

inline
const lodpal_parameters&
lodpal_analysis::previous_parameters() const
{
  return my_previous_parameters;
}

inline
lodpal_pairs&
lodpal_analysis::pairs_info()
{
  return my_pairs;
}

inline
const lodpal_pairs&
lodpal_analysis::pairs_info() const
{
  return my_pairs;
}

inline
size_t
lodpal_analysis::valid_parameter_count() const
{
   return my_valid_parameter_count;
}

inline
bool
lodpal_analysis::built() const
{
   if(!parameters().valid())
     return false;
   return my_built;
}

inline
double
lodpal_analysis::lod_score() const
{
  return my_lod_score;
}

inline
int
lodpal_analysis::function_evaluations() const
{
  return my_function_evaluations;
}

inline
int
lodpal_analysis::last_error() const
{
  return my_last_error;
}

inline
int
lodpal_analysis::ivfl() const
{
  return my_ivfl;
}

inline
int
lodpal_analysis::igage() const
{
  return my_igage;
}

inline
char
lodpal_analysis::method() const
{
  return my_method;
}

inline
const Matrix2D<double>&
lodpal_analysis::var_cov_matrix() const
{
  return my_vc_matrix;
}

inline
const RelativePairs&
lodpal_analysis::relative_pairs() const
{
  return my_pairs.relative_pairs();
}

inline
const marker_type& 
lodpal_analysis::current_marker() const
{
  return my_current_marker;
}

inline
void
lodpal_analysis::set_current_marker(const marker_type& m)
{
  my_current_marker = m;
}

inline
double
lodpal_analysis::lod_score_mm() const
{
  return my_lod_score_mm;
}

inline
double
lodpal_analysis::lod_score_mf() const
{
  return my_lod_score_mf;
}

inline
double
lodpal_analysis::lod_score_ff() const
{
  return my_lod_score_ff;
}

inline
bool
lodpal_analysis::is_valid_for_emp_p_value() const
{
  if( relative_pairs().is_x_linked(parameters().marker_parameters(0).marker) )
    return false;

  if(    pairs_info().parameters().autosomal_model().model != autosomal_model_type::one_parameter
      || pairs_info().parameters().autosomal_model().parent_of_origin )
    return false;

  if( pairs_info().parameters().trait_parameters(0).pair_select != trait_parameter::conaff )
    return false;

  if( pairs_info().pair_count() < 20 || pairs_info().pair_count() > 350 )
    return false;

  if( pairs_info().parameters().covariate_count() > 4 )
    return false;

  if( pairs_info().parameters().multipoint() &&
      (relative_pairs().get_average_marker_distance() < 1.0 ||
       relative_pairs().get_average_marker_distance() > 20.0) )
    return false;

  return true;
}
