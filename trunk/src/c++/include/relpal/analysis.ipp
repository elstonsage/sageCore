////////////////////////////////////////////////////////////////////////////
//             Implementation of relpal_analysis.h (Inline)               //
////////////////////////////////////////////////////////////////////////////

inline const relative_pairs*
relpal_analysis::get_pairs() const
{
  return my_pairs;
}

inline const regression_type&
relpal_analysis::get_regression_type() const
{
  return my_reg_type;
}

inline const regression_model&
relpal_analysis::get_current_model() const
{
  return my_current_model;
}

inline const regression_model&
relpal_analysis::get_current_null_model() const
{
  return my_current_null_model;
}

inline const analysis_result&
relpal_analysis::get_current_result() const
{
  return my_current_result;
}

inline size_t
relpal_analysis::get_member_count() const
{
  return my_member_count;
}

inline size_t
relpal_analysis::get_pair_count() const
{
  return my_pair_count;
}

inline size_t
relpal_analysis::get_fsib_pair_count() const
{
  return my_fsib_pair_count;
}

inline size_t
relpal_analysis::get_hsib_pair_count() const
{
  return my_hsib_pair_count;
}

inline const data_options&
relpal_analysis::get_data_options() const
{
  return my_data_opt;
}

inline const output_options&
relpal_analysis::get_output_options() const
{
  return my_output_opt;
}
inline void
relpal_analysis::set_data_options(const relpal_parser& r_parser)
{
  my_data_opt = r_parser.get_data_options();

  return;
}
