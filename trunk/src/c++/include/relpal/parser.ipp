// ---------------------------------------------------------------------------
// Inline Implementation of relpal_parser
// ---------------------------------------------------------------------------

inline void
relpal_parser::set_default_use_members(bool mibd)
{
  if( my_data_opt.use_members != DEFAULT )
    return;

  if( do_first_level_test() )
    my_data_opt.use_members = EVERY;
  else
    if( mibd )
      my_data_opt.use_members = INFORMATIVE_REGION;
    else
      my_data_opt.use_members = INFORMATIVE_LOCAL;

  return;
}

inline bool
relpal_parser::do_ind_batch_test() const
{
  return my_ind_batch_test;
}

inline bool
relpal_parser::do_first_level_test() const
{
  if( do_ind_batch_test() )
    return true;

  if( get_regression_type() == STZM || get_regression_type() == MTZM )
    return true;

  return false;
}

inline const regression_type&
relpal_parser::get_regression_type() const
{
  return my_reg_type;
}

inline const vector<dependent_variable>&
relpal_parser::get_traits() const
{
  return my_traits;
}

inline const vector<covariate_type>&
relpal_parser::get_ind_covariates() const
{
  return my_ind_covariates;
}

inline const vector<covariate_type>&
relpal_parser::get_ind_batch_covariates() const
{
  return my_ind_batch_covariates;
}

inline const vector<covariate_type>&
relpal_parser::get_ped_null_covariates() const
{
  return my_ped_null_covariates;
}

inline const vector<covariate_type>&
relpal_parser::get_ped_test_covariates() const
{
  return my_ped_test_covariates;
}

inline const vector<marker_type>&
relpal_parser::get_null_markers() const
{
  return my_null_markers;
}

inline const vector<marker_type>&
relpal_parser::get_test_markers() const
{
  return my_test_markers;
}

inline const vector<interaction_type>&
relpal_parser::get_ind_interactions() const
{
  return my_ind_interactions;
}

inline const vector<interaction_type>&
relpal_parser::get_ped_null_interactions() const
{
  return my_ped_null_interactions;
}

inline const vector<interaction_type>&
relpal_parser::get_ped_test_interactions() const
{
  return my_ped_test_interactions;
}

inline const data_options&
relpal_parser::get_data_options() const
{
  return my_data_opt;
}

inline const analysis_options&
relpal_parser::get_analysis_options() const
{
  return my_analysis_opt;
}

inline const pvalue_options&
relpal_parser::get_pvalue_options() const
{
  return my_pvalue_opt;
}

inline const output_options&
relpal_parser::get_output_options() const
{
  return my_output_opt;
}
