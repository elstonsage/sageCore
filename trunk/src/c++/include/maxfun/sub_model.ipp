//============================================================================
// File:      sub_model.ipp
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   gcw      Initial Implementation                    Apr 2001 
//            Baechle  reformatted.  Added various functions.    Jun 2001
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  maxfun_parameter
//============================================================================
//
inline
maxfun_parameter::maxfun_parameter(string ident, pstatus s, double val, 
                                   double l_bound, double u_bound )
    : value(val), lower_bound(l_bound), upper_bound(u_bound), status(s), identifier(ident)
{}

inline
maxfun_parameter::maxfun_parameter(const maxfun_parameter& other)
{
  value = other.value;
  lower_bound = other.lower_bound;
  upper_bound = other.upper_bound;
  status = other.status;
  identifier = other.identifier;
}

inline maxfun_parameter&
maxfun_parameter::operator=(const maxfun_parameter& other)
{
  if(this != &other)
  {
    value = other.value;
    lower_bound = other.lower_bound;
    upper_bound = other.upper_bound;
    status = other.status;
    identifier = other.identifier;
  }
  
  return *this;
}


//============================================================================
// IMPLEMENTATION:  sub_model
//============================================================================
//
inline
sub_model::sub_model(cerrorstream& errors)
    : my_in_use(false), my_errors(errors)
{}

inline
sub_model::sub_model(const sub_model& other)
{
  my_in_use = false;
  my_errors = other.my_errors;
  my_parameters = other.my_parameters;
}

inline sub_model&
sub_model::operator=(const sub_model& other)
{
  if(this != &other)
  {
    my_in_use = false;
    my_errors = other.my_errors;
    my_parameters = other.my_parameters;
  }
  
  return *this;
}

inline
sub_model::~sub_model()
{}

inline bool
sub_model::in_use() const
{
  return my_in_use;
}

inline string
sub_model::option_description() const
{
  return "";
}

inline string
sub_model::name() const
{
  return "";
}

inline size_t
sub_model::parameter_count() const
{
  return my_parameters.size();
} 

inline double
sub_model::parameter(size_t i) const
{
  if(i < my_parameters.size())
  {
    return my_parameters[i].value;
  }
  else
  {
    return QNAN;
  }
}

inline string
sub_model::parameter_identifier(size_t i) const
{
  if(i < my_parameters.size())
  {
    return my_parameters[i].identifier;
  }
  else
  {
    return "";
  }
}

// - Overridden for those sub-models that require parameter
//   rebuilding.
//
inline const std::vector<maxfun_parameter>&   
sub_model::get_parameters()
{
  return my_parameters;
}

// - Note: sub_model and sequencer are tightly coupled.
//   sub_model 'trusts' the sequencer to insure that it won't
//   'run off the end' of maxfun's parameter container.
//
inline int 
sub_model::synchronize(parameter_iterator start)
{
  for(size_t i = 0; i < my_parameters.size(); ++i)
  {
    my_parameters[i].value = start[i];
  }

  return 0;
}


//============================================================================
// IMPLEMENTATION:  sub_model_sequencer
//============================================================================
//
inline
sub_model_sequencer::sub_model_reference::sub_model_reference(sub_model* m, size_t i)
  : mod(m), pbegin(i)
{}

inline
sub_model_sequencer::sub_model_sequencer_function::sub_model_sequencer_function
  (MaxFunction& mf, sub_model_list& l)
    : my_maxfunction(mf), my_sub_models(l)
{ }
      
inline
sub_model_sequencer::sub_model_sequencer_function::~sub_model_sequencer_function()
{ }

inline double 
sub_model_sequencer::sub_model_sequencer_function::evaluate(parameter_vector& theta)
{
  return my_maxfunction(theta);
}

inline int    
sub_model_sequencer::sub_model_sequencer_function::update_bounds(parameter_vector& theta)
{
  sub_model_list::iterator i = my_sub_models.begin();

  for( ; i != my_sub_models.end(); ++i)
  {
    int v = i->mod->synchronize(theta.begin() + i->pbegin);

    if(v) return v;
  }

  return my_maxfunction.depar(theta);
}

inline
sub_model_sequencer::sub_model_sequencer(MaxFunction& mf)
  : my_sub_models(), my_function(mf, my_sub_models), my_maxfun(my_function)
{}

inline
sub_model_sequencer::~sub_model_sequencer()
{
  sub_model_list::iterator  iter;
  for(iter = my_sub_models.begin(); iter != my_sub_models.end(); ++iter)
  {
    iter->mod->my_in_use = false;
  }
}

inline Maxfun&
sub_model_sequencer::maxfun()
{
  return my_maxfun;
}

inline const Maxfun&
sub_model_sequencer::maxfun() const
{
  return my_maxfun;
}

inline void 
sub_model_sequencer::add_sub_model(sub_model* m)
{
  if(m->in_use())
  {
    return;
  }

  m->my_in_use = true;

  int params = my_maxfun.nt();

  my_sub_models.push_back(sub_model_reference(m, params));

  const vector<maxfun_parameter>& p = m->get_parameters();

  my_maxfun.nt() = params + p.size();

  // Copy the initial values into maxfun
  for(size_t i = 0; i < p.size(); ++i)
  {
    my_maxfun.istin(i + params) = p[i].status;
    my_maxfun.thin (i + params) = p[i].value;
    my_maxfun.thl  (i + params) = p[i].lower_bound;
    my_maxfun.thu  (i + params) = p[i].upper_bound;

    my_maxfun.set_label(i + params, p[i].identifier);
  } 
}


// - For debugging.
//
inline std::string
pstatus_2_string(pstatus status)
{
  switch(status)
  {
    case indep_func:
      return "independent function";
    case indep_non_func:
      return "independent non-function";
    case dependent:
      return "dependent";
    case fixed:
      return "fixed";
    default:
      return "";
  }
}

// - For debugging.
//
inline std::ostream&  
operator<<(std::ostream& out, const maxfun_parameter& max_par)
{
  int  old_precision = out.precision();
  out.precision(12);

  out << "maxfun parameter: " << max_par.identifier << "\n"
      << "value:            " << max_par.value   << "\n"
      << "status:           " << pstatus_2_string(max_par.status) << "\n"
      << "lower bound:      " << max_par.lower_bound << "\n"
      << "upper bound:      " << max_par.upper_bound << "\n" << std::endl;
  
  out.precision(old_precision);
  
  return out;
}

inline std::ostream&
operator<<(std::ostream& out, const sub_model& sm)
{
  out << "\n" << sm.name() << " values: \n";
  out << "Option: " << sm.option_description() << std::endl;
  vector<maxfun_parameter>::const_iterator  iter;
  for(iter = sm.my_parameters.begin(); iter != sm.my_parameters.end(); ++iter)
  {
    out << (*iter);
  }
  
  return out;
}

// - Return appropriate message fragment depending on whether specified
//   parameter is fixed or not.
//
inline std::string  
value_phrase(const maxfun_parameter& mp)
{
  std::ostringstream  return_buff; 
  return_buff << (mp.status == fixed ? "fixed value of " : "initial estimate of ");
  return_buff << mp.value;
  
  return return_buff.str();  
}

//============================================================================
// IMPLEMENTATION:  model_input
//============================================================================
//
inline
model_input::model_input(double i, bool f)
    : value(i), fixed(f)
{}


inline
model_input::model_input(const model_input& other)
{
  value = other.value;
  fixed = other.fixed;
}

inline model_input&  
model_input::operator=(const model_input& other)
{
  if(this != &other)
  {
    value = other.value;
    fixed = other.fixed;
  }
    
  return *this;
}

