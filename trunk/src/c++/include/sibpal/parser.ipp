// ---------------------------------------------------------------------------
// Inline Implementation of meantest_parser
// ---------------------------------------------------------------------------

inline const meantest_parameters&
meantest_parser::parameters() const
{
  return my_parameters;
}

inline bool
meantest_parser::wide_output() const
{
  return my_wide_output;
}

inline bool
meantest_parser::csv_output() const
{
  return my_csv_output;
}

inline void
meantest_parser::set_wide_output(bool b)
{
  my_wide_output = b;
}

inline void
meantest_parser::set_csv_output(bool b)
{
  my_csv_output = b;
}

// ---------------------------------------------------------------------------
// Inline Implementation of regression_parser
// ---------------------------------------------------------------------------

inline void
regression_parser::set_regression_type(regression_type rt)
{
  my_reg_type = rt;
}

inline const vector<LSF_ptr<LSFBase> >&
regression_parser::get_pair_info_file() const
{
  return my_pair_info_file;
}

inline const vector<trait_type>&
regression_parser::get_traits() const
{
  return my_traits;
}

inline const vector<trait_type>&
regression_parser::get_subsets() const
{
  return my_subsets;
}

inline const vector<covariate_type>&
regression_parser::get_covariates() const
{
  return my_covariates;
}

inline const vector<marker_type>&
regression_parser::get_markers() const
{
  return my_markers;
}

inline const vector<independent_variable>&
regression_parser::get_interactions() const
{
  return my_interactions;
}

inline const vector<independent_variable>&
regression_parser::get_batch_interactions() const
{
  return my_batch_interactions;
}

inline const regression_type&
regression_parser::get_regression_type() const
{
  return my_reg_type;
}

inline const string&
regression_parser::get_regression_method_name() const
{
  return my_reg_method_name;
}

inline const data_options&
regression_parser::get_data_options() const
{
  return my_data_opt;
}

inline const analysis_options&
regression_parser::get_analysis_options() const
{
  return my_analysis_opt;
}

inline const pvalue_options&
regression_parser::get_pvalue_options() const
{
  return my_pvalue_opt;
}

inline const output_options&
regression_parser::get_output_options() const
{
  return my_output_opt;
}
