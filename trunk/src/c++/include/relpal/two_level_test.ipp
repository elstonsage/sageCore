////////////////////////////////////////////////////////////////////////////
//             Implementation of two_level.h (Inline)                     //
////////////////////////////////////////////////////////////////////////////

inline void
two_level_base::set_pairs(const relative_pairs* r)
{
  my_pairs = r;
}

inline void
two_level_base::set_regression_type(const regression_type& r)
{
  my_reg_type = r;
}

inline void
two_level_base::set_model(const regression_model& m)
{
  my_model = m;
}

inline void
two_level_base::set_data(const analysis_data_type& d)
{
  my_data = d;
}

inline void
two_level_base::set_data_out(ostream* o)
{
  my_data_out = o;
}

inline void
two_level_base::set_debug_out(ostream* o)
{
  my_debug_out = o;
}

inline const relative_pairs*
two_level_base::get_pairs() const
{
  return my_pairs;
}

inline const regression_type&
two_level_base::get_regression_type() const
{
  return my_reg_type;
}

inline const regression_model&
two_level_base::get_model() const
{
  return my_model;
}

inline const relpal_least_square&
two_level_base::get_gls_1() const
{
  return my_gls_1;
}

inline ostream&
two_level_base::data_out()
{
  return *my_data_out;
}

inline ostream&
two_level_base::debug_out()
{
  return *my_debug_out;
}

//
//-----------------------------------------------------
//

inline const relpal_least_square&
two_level_regression::get_gls_2() const
{
  return my_gls_2;
}

//
//-----------------------------------------------------
//

inline const relpal_score&
two_level_score_test::get_score_2() const
{
  return my_score_2;
}

inline const vector<matrix>&
two_level_score_test::get_residuals() const
{
  return my_residuals;
}

inline double
two_level_score_test::get_correction(size_t type) const
{
  if( type == 0 )
    return my_correction_na;
  else if( type == 2 )
    return my_correction_al;

  return my_correction_sw;
}

inline double
two_level_score_test::get_emp_pvalue(size_t type) const
{
  if( type == 0 )
    return my_emp_pvalue_na;
  else if( type == 2 )
    return my_emp_pvalue_al;

  return my_emp_pvalue_sw;
}

inline size_t
two_level_score_test::get_rep_count(size_t type) const
{
  if( type == 0 )
    return my_rep_count_na;
  else if( type == 2 )
    return my_rep_count_al;

  return my_rep_count_sw;
}

inline bool
two_level_score_test::has_valid_null() const
{
  return my_valid_null;
}

inline bool
two_level_score_test::has_reliable_score() const
{
  return my_reliable_score;
}
