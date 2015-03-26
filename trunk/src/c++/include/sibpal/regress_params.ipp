////////////////////////////////////////////////////////////////////////////
//             Implementation of regress_params.h (Inline)                //
////////////////////////////////////////////////////////////////////////////

//
//------------------------------------------------------------------------
//

inline bool
trait_type::operator==(const trait_type& t) const
{
  if( trait_index == (size_t)-1 || t.trait_index == (size_t)-1 )
    return false;

  return trait_index == t.trait_index;
}

inline bool
trait_type::operator!=(const trait_type& t) const
{
  return !( (*this) == t );
}


//
//------------------------------------------------------------------------
//

inline double
dependent_variable::get_mean() const
{
  double mean = fixed_mean;

  if( !finite(mean) )
    mean = trait_all_sibs_info.mean();

  return mean;
}

inline bool
dependent_variable::operator==(const dependent_variable& t) const
{
  if( trait_index == (size_t)-1 || t.trait_index == (size_t)-1 )
    return false;

  return trait_index == t.trait_index;
}

inline bool
dependent_variable::operator!=(const dependent_variable& t) const
{
  return !( (*this) == t );
}

//
//------------------------------------------------------------------------
//

inline bool
covariate_type::operator==(const covariate_type& c) const
{
  return    covariate_index == c.covariate_index
         && power == c.power
         && operation == c.operation;
}

inline bool
covariate_type::operator!=(const covariate_type& c) const
{
  return !( (*this) == c );
}

inline bool
covariate_type::operator<(const covariate_type& c) const
{
  if(covariate_index < c.covariate_index)
    return true;
  if(covariate_index > c.covariate_index)
    return false;
  if(power < c.power)
    return true;
  if(power > c.power)
    return false;
  return operation < c.operation;
}

//
//------------------------------------------------------------------------
//

inline bool
marker_type::operator==(const marker_type& m) const
{
  return marker_index == m.marker_index && effect == m.effect;
}

inline bool
marker_type::operator!=(const marker_type& m) const
{
  return !( (*this) == m);
}

inline bool
marker_type::operator<(const marker_type& m) const
{
  if(marker_index < m.marker_index)
    return true;
  if(marker_index > m.marker_index)
    return false;
  return (effect < m.effect);
}

//
//------------------------------------------------------------------------
//

inline bool
independent_variable::operator==(const independent_variable& c) const
{
    return covariates == c.covariates && markers == c.markers && type == c.type;
}

inline bool
independent_variable::operator!=(const independent_variable& c) const
{
  return !( (*this) == c );
}

inline bool
independent_variable::operator<(const independent_variable& c) const
{
  if(covariates < c.covariates)
    return true;
  if(covariates > c.covariates)
    return false;
  return (markers < c.markers);
}

//
//------------------------------------------------------------------------
//

inline bool
regression_model::valid() const
{
  return my_valid;
}

inline bool
regression_model::is_x_linked() const
{
  return my_x_linked;
}

inline size_t
regression_model::get_parameter_count() const
{
  return my_parameters.size();
}

inline size_t
regression_model::get_subset_count() const
{
  return my_subsets.size();
}

inline dependent_variable&
regression_model::get_trait()
{
  return my_trait;
}

inline independent_variable&
regression_model::get_parameter(size_t i)
{
  return my_parameters[i];
}

inline trait_type&
regression_model::get_subset(size_t t)
{
  return my_subsets[t];
}

inline const dependent_variable&
regression_model::get_trait() const
{
  return my_trait;
}

inline const independent_variable&
regression_model::get_parameter(size_t i) const
{
  return my_parameters[i];
}

inline const trait_type&
regression_model::get_subset(size_t t) const
{
  return my_subsets[t];
}

inline regression_type
regression_model::get_regression_type() const
{
  return my_reg_type;
}

inline string
regression_model::get_regression_method_name() const
{
  return my_reg_method_name;
}

inline const data_options&
regression_model::get_data_options() const
{
  return my_data_opt;
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

inline const output_options&
regression_model::get_output_options() const
{
  return my_output_opt;
}

//------------------------------------------------------------

inline regression_model::parameter_iterator
regression_model::parameter_begin()
{
   return my_parameters.begin();
}

inline regression_model::parameter_iterator
regression_model::parameter_end()
{
   return my_parameters.end();
}

inline regression_model::parameter_const_iterator
regression_model::parameter_begin() const
{
   return my_parameters.begin();
}

inline regression_model::parameter_const_iterator
regression_model::parameter_end() const
{
   return my_parameters.end();
}

//------------------------------------------------------------

inline void
regression_model::set_trait(size_t t, bool binary)
{
  my_trait = dependent_variable(t, binary);
  invalidate();
}

inline void
regression_model::set_trait(const dependent_variable& t)
{
  my_trait = t;
  invalidate();
}

//------------------------------------------------------------

inline void
regression_model::set_regression_type(regression_type r)
{
  my_reg_type = r;
}

inline void
regression_model::set_regression_method_name(string n)
{
  my_reg_method_name = n;
}

inline void
regression_model::set_data_options(const data_options& op)
{
  my_data_opt = op;
}

inline void
regression_model::set_analysis_options(const analysis_options& op)
{
  my_analysis_opt = op;
}

inline void
regression_model::set_pvalue_options(const pvalue_options& op)
{
  my_pvalue_opt = op;
}

inline void
regression_model::set_output_options(const output_options& op)
{
  my_output_opt = op;
}

//------------------------------------------------------------

inline void
regression_model::clear_subsets()
{
  if( my_subsets.size() )
  {
    my_subsets.clear();
    invalidate();
  }
}

inline pair<regression_model::trait_iterator, bool>
regression_model::add_subset(size_t t)
{
  return add_subset( trait_type(t) );
}

inline pair<regression_model::trait_iterator, bool>
regression_model::set_subset(size_t t)
{
  invalidate();
  clear_subsets();

  return add_subset(t);
}

//------------------------------------------------------------

inline void
regression_model::clear_parameters()
{
  invalidate();
  my_parameters.clear();
}

inline regression_model::pib_value
regression_model::add_marker(size_t m, marker_type::effect_type effect, bool x)
{
  independent_variable param;
  param.type = independent_variable::MARKER;
  param.markers.push_back( marker_type(m, effect, x) );

  return add_parameter(param);
}

inline regression_model::pib_value
regression_model::add_covariate(size_t c, covariate_type::cov_op op, double power, double fixed_mean)
{
  independent_variable param;
  param.type = independent_variable::COVARIATE;
  param.covariates.push_back( covariate_type(c, op, power, fixed_mean) );

  return add_parameter(param);
}

