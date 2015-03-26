//============================================================================
// File:      mle_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/4/2 - created.                                djb
//                                                                          
// Notes:     inlines for mle_sub_model class.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



inline std::ostream&
operator<<(std::ostream& out, const mle_sub_model& mle)
{
  out << std::boolalpha
      << "\n" << mle.name() << " values - \n"
      << "Sex-specific:           " << mle.sex_specific << "\n"
      << "Uses alpha:             " << mle.use_alpha    << "\n" 
      << std::endl;
      
  out << "INPUT:" << endl;
  vector<MAXFUN::ParameterInput>::const_iterator  iter;
  for(iter = mle.my_parameters.begin(); iter != mle.my_parameters.end(); ++iter)
  {
    out << "\ngroup_name        " << iter->group_name
        << "\nparam_name        " << iter->param_name
        << "\ninitial_estimate  " << iter->initial_estimate
        << "\ninitial_type      " << MAXFUN::ParamTypeEnum2str(iter->initial_type)
        << "\nlower_bound       " << iter->lower_bound
        << "\nupper_bound       " << iter->upper_bound
        << endl;
  }
  
  const MAXFUN::ParameterMgr*  parameter_manager = mle.getParameterMgr();
  if(parameter_manager)
  {
    out << "\n\nFINAL:" << endl;
    int  param_count = parameter_manager->getParamCount();
    for(int i = 0; i < param_count; ++i)
    {
      const MAXFUN::Parameter&  parameter = parameter_manager->getParameter(i);
      out << "name " << parameter.getName() << endl;
      out << "type " << MAXFUN::ParamTypeEnum2str(parameter.getFinalType()) << endl;
      out << "lower bound " << parameter.getLowerBound() << endl;
      out << "upper bound " << parameter.getUpperBound() << endl;
      out << "value " << parameter.getFinalEstimate() << endl;
    }
  }
  
  out << noboolalpha << std::endl;
  
  return out;
}

//============================================================================
// IMPLEMENTATION:  mle_sub_model
//============================================================================
//
inline  
mle_sub_model::mle_sub_model(bool ss, bool ua, cerrorstream& errors)
    : MAXFUN::Submodel(errors)
{
  set(ss, ua);
}

inline string  
mle_sub_model::option_description() const
{
  return  "N/A";
}

inline string  
mle_sub_model::name() const
{
  return  MLE_NAME;
}

inline const vector<MAXFUN::ParameterInput>&
mle_sub_model::parameters() const
{
  return  my_parameters;
}

inline double
mle_sub_model::average_theta() const
{
  return  my_average_theta;
}

inline double
mle_sub_model::male_theta() const
{
  return  my_male_theta;
}

inline double
mle_sub_model::female_theta() const
{
  return  my_female_theta;
}

inline double
mle_sub_model::average_theta_ub() const
{
  if(sex_specific)
  {
    return  QNAN;
  }
  else
  {
    assert(AVERAGE < my_parameters.size());
    return  my_parameters[AVERAGE].upper_bound;
  }
}

inline double
mle_sub_model::male_theta_ub() const
{
  if(sex_specific)
  {
    assert(MALE < my_parameters.size());
    return  my_parameters[MALE].upper_bound;
  }
  else
  {
    return  QNAN;
  }
}

inline double
mle_sub_model::female_theta_ub() const
{
  if(sex_specific)
  {
    assert(FEMALE < my_parameters.size());
    return  my_parameters[FEMALE].upper_bound;
  }
  else
  {
    return  QNAN;
  }
}

inline double
mle_sub_model::theta(LODLINK::sex s) const
{
  double  t;
  
  if(sex_specific)
  {
    if(s == LODLINK::male)
    {
      t = my_male_theta;
    }
    else
    {
      t = my_female_theta;
    }
  }
  else
  {
    t = my_average_theta;
  }
  
  return t;
}

inline double
mle_sub_model::alpha() const
{
  return  my_alpha;
}

inline bool
mle_sub_model::is_sex_specific() const
{
  return  sex_specific;
}

inline bool
mle_sub_model::uses_alpha() const
{
  return  use_alpha;
}

inline const MAXFUN::ParameterMgr* 
mle_sub_model::parameter_mgr() const
{
  return  getParameterMgr();
}

inline void
mle_sub_model::reset()
{
  set();
}

// - Restrict recombination fractions to being no greater than .5
//
inline void
mle_sub_model::set_strict_limits()
{
  assert(! isLinked());

  if(sex_specific)
  {
    my_parameters[MALE].upper_bound = THETA_STRICT_UPPER_BOUND;
    my_parameters[FEMALE].upper_bound = THETA_STRICT_UPPER_BOUND;
  }
  else
  {
    my_parameters[AVERAGE].upper_bound = THETA_STRICT_UPPER_BOUND;
  }
}

// - Allow recombination fractions to be as large as 1.
//
inline void
mle_sub_model::set_relaxed_limits()
{
  assert(! isLinked());

  if(sex_specific)
  {
    my_parameters[MALE].upper_bound = THETA_UPPER_BOUND;
    my_parameters[FEMALE].upper_bound = THETA_UPPER_BOUND;
  }
  else
  {
    my_parameters[AVERAGE].upper_bound = THETA_UPPER_BOUND;
  }
}

// - Fix alpha at 1.
//
inline void
mle_sub_model::fix_alpha()
{
  assert(! isLinked());
  assert(use_alpha);
  
  int  alpha_idx = sex_specific ? static_cast<int>(ALPHA_TWO) : 
                                  static_cast<int>(ALPHA_ONE);
  
  my_parameters[alpha_idx].initial_estimate = 1;
  my_parameters[alpha_idx].initial_type = MAXFUN::Parameter::FIXED;
  
  init();
}

// - Remove fixed parameter type from alpha.
//
inline void
mle_sub_model::unfix_alpha()
{
  assert(! isLinked());
  assert(use_alpha);

  int  alpha_idx = sex_specific ? static_cast<int>(ALPHA_TWO) : 
                                  static_cast<int>(ALPHA_ONE);
  
  my_parameters[alpha_idx].initial_estimate = ALPHA_INIT_VALUE;
  my_parameters[alpha_idx].initial_type = ALPHA_DEFAULT_STATUS;
    
  init();
}

// - Impose the constraint that the sex specific recombinations sum to 1.
//
inline void
mle_sub_model::constrain_thetas()
{
  assert(! isLinked());
  assert(sex_specific);
  assert(my_parameters[MALE].upper_bound == THETA_UPPER_BOUND);
  assert(my_parameters[FEMALE].upper_bound == THETA_UPPER_BOUND);
  
  my_parameters[MALE].initial_type = MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
  my_parameters[FEMALE].initial_type = MAXFUN::Parameter::DEPENDENT;
  my_parameters[MALE].initial_estimate = .5;
  my_parameters[FEMALE].initial_estimate = .5;
  
  init();
}

// - Remove constraints on sex specific recombination fractions.
//
inline void
mle_sub_model::unconstrain_thetas()
{
  assert(! isLinked());
  assert(sex_specific);
  
  my_parameters[MALE].initial_type = MAXFUN::Parameter::INDEPENDENT;
  my_parameters[FEMALE].initial_type = MAXFUN::Parameter::INDEPENDENT;
}

// - Set value of non sex specific recombination fraction.
//
inline void
mle_sub_model::set_average_theta(double value)
{
  assert(! isLinked());
  assert(! sex_specific);
  assert(my_parameters[AVERAGE].lower_bound <= value &&
         value <= my_parameters[AVERAGE].upper_bound   );
         
  my_parameters[AVERAGE].initial_estimate = value;
  
  init();
}

// - Set value of male recombination fraction.
//
inline void
mle_sub_model::set_male_theta(double value)
{
  assert(! isLinked());
  assert(sex_specific);
  assert(my_parameters[FEMALE].initial_type != MAXFUN::Parameter::DEPENDENT);
  assert(my_parameters[MALE].lower_bound <= value &&
         value <= my_parameters[MALE].upper_bound   );
         
  my_parameters[MALE].initial_estimate = value;
  
  init();
}

// - Set value of female recombination fraction.
//
inline void
mle_sub_model::set_female_theta(double value)
{
  assert(! isLinked());
  assert(sex_specific);
  assert(my_parameters[FEMALE].initial_type != MAXFUN::Parameter::DEPENDENT);
  assert(my_parameters[FEMALE].lower_bound <= value &&
         value <= my_parameters[FEMALE].upper_bound   );
         
  my_parameters[FEMALE].initial_estimate = value;
  
  init();
}

// - Set proportions of families w. linkage.
//
inline void
mle_sub_model::set_alpha(double value)
{
  assert(! isLinked());
  assert(use_alpha);
  
  size_t  index = sex_specific ? static_cast<int>(ALPHA_TWO) : 
                                 static_cast<int>(ALPHA_ONE);
  
  assert(my_parameters[index].lower_bound <= value &&
         value <= my_parameters[index].upper_bound   );
         
  my_parameters[index].initial_estimate = value;
  
  init();
}





