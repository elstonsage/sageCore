//===========================================================
//
//  File:	ParameterInput.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2004 R. C. Elston
//===========================================================

#include "maxfunapi/ParameterInput.h"

namespace SAGE   {
namespace MAXFUN {

ParameterInput::ParameterInput()
{
  group_name       =  "";
  param_name       =  "";
  initial_estimate =  0.0;
  initial_type     =  Parameter::INDEPENDENT;
  lower_bound      = -MF_INFINITY;
  upper_bound      =  MF_INFINITY;
  index            =  0;
}

ParameterInput::ParameterInput(string                   _group_name,
                               string                   _param_name,
                               Parameter::ParamTypeEnum _initial_type, 
                               double                   _initial_estimate,
                               double                   _lower_bound,
                               double                   _upper_bound)

  : group_name       (_group_name),
    param_name       (_param_name),
    initial_estimate (_initial_estimate),
    initial_type     (_initial_type),
    lower_bound      (_lower_bound),
    upper_bound      (_upper_bound)
{}

ParameterInput::ParameterInput(const ParameterInput& other)
{
  copy(other);
}

ParameterInput& 
ParameterInput::operator=(const ParameterInput& other)
{
  if(&other != this)
  {
    copy(other);
  }

  return *this;
}

void 
ParameterInput::copy(const ParameterInput& other)
{
  group_name       = other.group_name;
  param_name       = other.param_name;
  initial_estimate = other.initial_estimate;
  initial_type     = other.initial_type;
  lower_bound      = other.lower_bound;
  upper_bound      = other.upper_bound;
  index            = other.index;
}

} // End namespace MAXFUN
} // End namespace SAGE
