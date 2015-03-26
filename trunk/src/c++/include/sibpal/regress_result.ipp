////////////////////////////////////////////////////////////////////////////
//             Implementation of regress_result.h (Inline)                //
////////////////////////////////////////////////////////////////////////////

//
//------------------------------------------------------------------------
//

inline void
parameter_estimate::set_test_direction(direction_type d)
{
  my_direction = d;
}

inline void
parameter_estimate::set_value(double v)
{
  my_value = v;
}

inline void
parameter_estimate::set_variance(double ss)
{
  my_variance = ss;
}

inline void
parameter_estimate::set_expected_value(double e)
{
  my_expected_value = e;
}

inline void
parameter_estimate::set_degrees_of_freedom(size_t df)
{
  my_degrees_of_freedom = df;
}

inline void
parameter_estimate::set_adjusted_tvalue(double t)
{
  my_adjusted_tvalue = t;
}

inline void
parameter_estimate::set_adjusted_degrees_of_freedom(size_t df)
{
  my_adjusted_degrees_of_freedom = df;
}

inline void
parameter_estimate::set_adjusted_pvalue(double p)
{
  my_adjusted_pvalue = p;
}

inline parameter_estimate::direction_type
parameter_estimate::test_direction() const
{
  return my_direction;
}

inline double
parameter_estimate::value() const
{
  return my_value;
}

inline double
parameter_estimate::variance() const
{
  return my_variance;
}

inline double
parameter_estimate::standard_error() const
{
  if( finite(my_variance) && my_variance >= 0.0 )
    return sqrt(my_variance);

  return numeric_limits<double>::quiet_NaN();
}

inline double
parameter_estimate::expected_value() const
{
  return my_expected_value;
}

inline size_t
parameter_estimate::degrees_of_freedom() const
{
  return my_degrees_of_freedom;
}

inline double
parameter_estimate::tvalue() const
{
  double p = adjusted_tvalue();
  if( !SAGE::isnan(p) )
    return p;

  return raw_tvalue();
}

inline double
parameter_estimate::raw_tvalue() const
{
  const double v  = value() - expected_value();
  const double ss = variance();

  if( ss < 0 )
    return std::numeric_limits<double>::quiet_NaN();

  // The following code relies heavily on the compiler following IEEE 754
  // Bad things(tm) will happen if it does not.

  return v/sqrt(ss);
}

inline double
parameter_estimate::adjusted_tvalue() const
{
  return my_adjusted_tvalue;
}

inline size_t
parameter_estimate::adjusted_degrees_of_freedom() const
{
  return my_adjusted_degrees_of_freedom;
}

inline double
parameter_estimate::raw_adjusted_pvalue() const
{
  return my_adjusted_pvalue;
}

inline double
parameter_estimate::empirical_pvalue() const
{
  if( my_significant_replicates + my_same_replicates + my_total_replicates == 0 )
    return std::numeric_limits<double>::quiet_NaN();

  double num   = (double)my_significant_replicates+1.0+(0.5*(double)my_same_replicates);
  double denom = (double)my_total_replicates+1.0;

  return num / denom;
}

inline size_t
parameter_estimate::total_replicates() const
{
  return my_total_replicates;
}

inline size_t
parameter_estimate::significant_replicates() const
{
  return my_significant_replicates;
}

inline size_t
parameter_estimate::same_replicates() const
{
  return my_same_replicates;
}

//
//------------------------------------------------------------
//

inline void
regression_results::set_residual_variance(double rv)
{
  my_variances.residual = rv;
}

inline void
regression_results::set_total_variance(double tv)
{
  my_variances.total = tv;
}

inline void
regression_results::set_sum_residual_variance(double rv)
{
  my_variances.sum_residual = rv;
}

inline void
regression_results::set_diff_residual_variance(double rv)
{
  my_variances.diff_residual = rv;
}

inline void
regression_results::set_full_residual_variance(double rv)
{
  my_full_variances.residual = rv;
}

inline void
regression_results::set_full_total_variance(double tv)
{
  my_full_variances.total = tv;
}

inline void
regression_results::set_full_sum_residual_variance(double rv)
{
  my_full_variances.sum_residual = rv;
}

inline void
regression_results::set_full_diff_residual_variance(double rv)
{
  my_full_variances.diff_residual = rv;
}

inline void
regression_results::set_half_residual_variance(double rv)
{
  my_half_variances.residual = rv;
}

inline void
regression_results::set_half_total_variance(double tv)
{
  my_half_variances.total = tv;
}

inline void
regression_results::set_half_sum_residual_variance(double rv)
{
  my_half_variances.sum_residual = rv;
}

inline void
regression_results::set_half_diff_residual_variance(double rv)
{
  my_half_variances.diff_residual = rv;
}

inline void
regression_results::set_full_w(const w_map& ws)
{
  my_full_w = ws;
}

inline void
regression_results::set_half_w(const w_map& ws)
{
  my_half_w = ws;
}

inline void
regression_results::set_residual_info(const SampleInfo& s)
{
  my_residual_info = s;
}

inline void
regression_results::set_residual_square_info(const SampleInfo& s)
{
  my_residual_square_info = s;
}

inline void
regression_results::set_residual_square_info_reduced(const SampleInfo& s)
{
  my_residual_square_info_reduced = s;
}

inline void
regression_results::set_full_residual_info(const SampleInfo& s)
{
  my_full_residual_info = s;
}

inline void
regression_results::set_half_residual_info(const SampleInfo& s)
{
  my_half_residual_info = s;
}

inline void
regression_results::add_result(const reg_result& r)
{
  my_results.push_back(r);
  if( r.type() == reg_result::intercept )
    ++my_intercept_count;
}

inline void
regression_results::set_results(const result_vector& r)
{
  my_results = r;
}

inline void
regression_results::set_F_result(const F_result_type& f)
{
  my_F_result = f;
}

inline void
regression_results::add_full_result(const reg_result& r)
{
  my_full_results.push_back(r);
}

inline void
regression_results::add_half_result(const reg_result& r)
{
  my_half_results.push_back(r);
}

//------------------------------------------------------------

inline double
regression_results::get_residual_variance() const
{
  return my_variances.residual;
}

inline double
regression_results::get_total_variance() const
{
  return my_variances.total;
}

inline double
regression_results::get_sum_residual_variance() const
{
  return my_variances.sum_residual;
}

inline double
regression_results::get_diff_residual_variance() const
{
  return my_variances.diff_residual;
}

inline double
regression_results::get_full_residual_variance() const
{
  return my_full_variances.residual;
}

inline double
regression_results::get_full_total_variance() const
{
  return my_full_variances.total;
}

inline double
regression_results::get_full_sum_residual_variance() const
{
  return my_full_variances.sum_residual;
}

inline double
regression_results::get_full_diff_residual_variance() const
{
  return my_full_variances.diff_residual;
}

inline double
regression_results::get_half_residual_variance() const
{
  return my_half_variances.residual;
}

inline double
regression_results::get_half_total_variance() const
{
  return my_half_variances.total;
}

inline double
regression_results::get_half_sum_residual_variance() const
{
  return my_half_variances.sum_residual;
}

inline double
regression_results::get_half_diff_residual_variance() const
{
  return my_half_variances.diff_residual;
}

inline const map<size_t, pair<double, double> >&
regression_results::get_full_w() const
{
  return my_full_w;
}

inline const map<size_t, pair<double, double> >&
regression_results::get_half_w() const
{
  return my_half_w;
}

inline const SampleInfo&
regression_results::get_residual_info() const
{
  return my_residual_info;
}

inline const SampleInfo&
regression_results::get_residual_square_info() const
{
  return my_residual_square_info;
}

inline const SampleInfo&
regression_results::get_residual_square_info_reduced() const
{
  return my_residual_square_info_reduced;
}

inline const SampleInfo&
regression_results::get_full_residual_info() const
{
  return my_full_residual_info;
}

inline const SampleInfo&
regression_results::get_half_residual_info() const
{
  return my_half_residual_info;
}

inline reg_result&
regression_results::get_result(size_t i)
{
  return my_results[i];
}

inline const reg_result&
regression_results::get_result(size_t i) const
{
  return my_results[i];
}

inline reg_result&
regression_results::get_full_result(size_t i)
{
  return my_full_results[i];
}

inline const reg_result&
regression_results::get_full_result(size_t i) const
{
  return my_full_results[i];
}

inline reg_result&
regression_results::get_half_result(size_t i)
{
  return my_half_results[i];
}

inline const reg_result&
regression_results::get_half_result(size_t i) const
{
  return my_half_results[i];
}

inline size_t
regression_results::get_result_count() const
{
  return my_results.size();
}

inline size_t
regression_results::get_intercept_count() const
{
  return my_intercept_count;
}

inline result_vector&
regression_results::get_results()
{
  return my_results;
}

inline const result_vector&
regression_results::get_results() const
{
  return my_results;
}

inline F_result_type&
regression_results::get_F_result()
{
  return my_F_result;
}

inline const F_result_type&
regression_results::get_F_result() const
{
  return my_F_result;
}
