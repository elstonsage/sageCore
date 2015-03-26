////////////////////////////////////////////////////////////////////////////
//  Implementation of regression_model.h (Inline)
////////////////////////////////////////////////////////////////////////////

inline bool
regression_model::valid() const
{
  return my_valid;
}

inline size_t
regression_model::get_trait_count() const
{
  return my_traits.size();
}

inline size_t
regression_model::get_ind_parameter_count() const
{
  return my_ind_parameters.size();
}

inline size_t
regression_model::get_ped_parameter_count() const
{
  return my_ped_parameters.size();
}

inline string
regression_model::get_name() const
{
  return my_name;
}

inline const analysis_options&
regression_model::get_analysis_options() const
{
  return my_analysis_opt;
}

inline const pvalue_options&
regression_model::get_pvalue_options() const
{
  return my_pvalue_opt;
}

inline dependent_variable&
regression_model::get_trait(size_t t)
{
  return my_traits[t];
}

inline independent_variable&
regression_model::get_ind_parameter(size_t i)
{
  return my_ind_parameters[i];
}

inline independent_variable&
regression_model::get_ped_parameter(size_t i)
{
  return my_ped_parameters[i];
}

inline const dependent_variable&
regression_model::get_trait(size_t t) const
{
  return my_traits[t];
}

inline const independent_variable&
regression_model::get_ind_parameter(size_t i) const
{
  return my_ind_parameters[i];
}

inline const independent_variable&
regression_model::get_ped_parameter(size_t i) const
{
  return my_ped_parameters[i];
}

inline trait_iterator
regression_model::trait_begin()
{
   return my_traits.begin();
}

inline trait_iterator
regression_model::trait_end()
{
   return my_traits.end();
}

inline parameter_iterator
regression_model::ind_parameter_begin()
{
   return my_ind_parameters.begin();
}

inline parameter_iterator
regression_model::ind_parameter_end()
{
   return my_ind_parameters.end();
}

inline parameter_iterator
regression_model::ped_parameter_begin()
{
   return my_ped_parameters.begin();
}

inline parameter_iterator
regression_model::ped_parameter_end()
{
   return my_ped_parameters.end();
}

inline trait_const_iterator
regression_model::trait_begin() const
{
   return my_traits.begin();
}

inline trait_const_iterator
regression_model::trait_end() const
{
   return my_traits.end();
}

inline parameter_const_iterator
regression_model::ind_parameter_begin() const
{
   return my_ind_parameters.begin();
}

inline parameter_const_iterator
regression_model::ind_parameter_end() const
{
   return my_ind_parameters.end();
}

inline parameter_const_iterator
regression_model::ped_parameter_begin() const
{
   return my_ped_parameters.begin();
}

inline parameter_const_iterator
regression_model::ped_parameter_end() const
{
   return my_ped_parameters.end();
}

inline void
regression_model::clear_traits()
{
  my_traits.clear();
  invalidate();
}

inline void
regression_model::clear_ind_parameters()
{
  invalidate();
  my_ind_parameters.clear();
}

inline void
regression_model::clear_ped_parameters()
{
  invalidate();
  my_ped_parameters.clear();
}

inline pib_value
regression_model::add_ind_covariate(const covariate_type& c)
{
  independent_variable param;

  param.type = independent_variable::COVARIATE;
  param.covariates.push_back(c);

  return add_ind_parameter(param);
}

inline pib_value
regression_model::add_ped_covariate(const covariate_type& c)
{
  independent_variable param;

  param.type = independent_variable::COVARIATE;
  param.covariates.push_back(c);

  return add_ped_parameter(param);
}

inline void
regression_model::add_intercept(size_t t)
{
  independent_variable param;

  param.type = independent_variable::INTERCEPT;
  param.t1   = param.t2 = t;

  my_ind_parameters.push_back(param);

  return;
}

inline void
regression_model::add_random_err_variance(size_t t1, size_t t2)
{
  independent_variable param;

  param.type = independent_variable::RANDOM_ERR;
  param.t1   = t1;
  param.t2   = t2;

  my_ped_parameters.push_back(param);

  return;
}

inline void
regression_model::add_common_env_variance(size_t t1, size_t t2)
{
  independent_variable param;

  param.type = independent_variable::COMMON_ENV;
  param.t1   = t1;
  param.t2   = t2;

  my_ped_parameters.push_back(param);

  return;
}

inline void
regression_model::add_polygenic_variance(size_t t1, size_t t2)
{
  independent_variable param;

  param.type = independent_variable::POLYGENIC;
  param.t1   = t1;
  param.t2   = t2;

  my_ped_parameters.push_back(param);

  return;
}

inline void
regression_model::set_analysis_options(const analysis_options& op)
{
  my_analysis_opt = op;

  return;
}

inline void
regression_model::set_pvalue_options(const pvalue_options& op)
{
  my_pvalue_opt = op;

  return;
}

inline void
regression_model::set_name(const string& s)
{
  my_name = s;

  return;
}

inline bool
regression_model::operator==(const regression_model& r) const
{
  if(    get_trait_count() != r.get_trait_count()
      || get_ind_parameter_count() != r.get_ind_parameter_count()
      || get_ped_parameter_count() != r.get_ped_parameter_count() )
    return false;


  for( size_t t = 0; t < get_trait_count(); ++t )
    if( get_trait(t) != r.get_trait(t) )
      return false;

  for( size_t i = 0; i < get_ind_parameter_count(); ++i )
    if( get_ind_parameter(i) != r.get_ind_parameter(i) )
      return false;

  for( size_t p = 0; p < get_ped_parameter_count(); ++p )
    if( get_ped_parameter(p) != r.get_ped_parameter(p) )
      return false;

  return true;
}
  
inline bool
regression_model::operator!=(const regression_model& r) const
{
  return !( (*this) == r );
}
