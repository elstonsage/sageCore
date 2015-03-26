//============================================================================
// File:      transformation_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/2/01 - created.                       -djb
//                                                                          
// Notes:     implementation of transformation_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <cmath>
#include "maxfun/transf_sub_model.h"
#include "numerics/log_double.h"

using namespace std;
namespace SAGE
{

const std::string  TRANSFORMATION_NAME   = "Transformation";
const double       LAMBDA_ONE_DEFAULT_VALUE = 1;
const bool         LAMBDA_ONE_DEFAULT_FIXED = false;
const double       LAMBDA_ONE_DEFAULT_LB = -1;             
const double       LAMBDA_ONE_DEFAULT_UB = std::numeric_limits<double>::infinity();
const double       LAMBDA_TWO_DEFAULT_VALUE = 0;
const bool         LAMBDA_TWO_DEFAULT_FIXED = true;

//============================================================================
// IMPLEMENTATION:  transformation_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
void
transformation_sub_model::set_option(sm_option opt)
{
     my_option = opt;
}

bool
transformation_sub_model::set(sm_option opt, const model_input& lambda1,
                              const model_input& lambda2, double lambda1_lb,
                              double lambda1_ub)
{
  clear_each_sync = true;

  if(my_in_use)
  {
    return false;
  }

  my_option = opt;

  if(opt == no_trans)
  {
    return set_none(lambda1, lambda2, lambda1_lb, lambda1_ub);
  }

  // - Set default values.
  //
  my_parameters.resize(2);
  my_parameters[0] = maxfun_parameter("lambda_one", indep_non_func,
                                      LAMBDA_ONE_DEFAULT_VALUE, 
                                      LAMBDA_ONE_DEFAULT_LB,
                                      LAMBDA_ONE_DEFAULT_UB);

  my_parameters[1] = maxfun_parameter("lambda_two", fixed, LAMBDA_TWO_DEFAULT_VALUE);

  bool  return_value = set_lambda_one_limits(lambda1, lambda1_lb, lambda1_ub) &&
                       set_lambda_one(lambda1)                                &&
                       set_lambda_two(lambda2);

  // - Sync. internally.
  //
  internally_synchronize();

  return return_value;
}

//
// ==============  Ancillary Functions
//
bool  
transformation_sub_model::set_none(const model_input& lambda1, const model_input& lambda2,
                                   double lambda1_lb, double lambda1_ub)
{
  if(! (SAGE::isnan(lambda1.value) && SAGE::isnan(lambda2.value) &&
        SAGE::isnan(lambda1_lb)    && SAGE::isnan(lambda1_ub)      ))
  {
    my_errors << priority(error) << "Unnecessary parameters specified in "
              << "transformation sub-block with 'none' option.  Ignoring ..." << endl;
  }
    
  my_parameters.resize(0);
  my_lambda_one = QNAN;
  my_lambda_two = QNAN;
    
  return true;
}

// - Power parameter.  Assumes power parameter limits are already set.
//
bool  
transformation_sub_model::set_lambda_one(const model_input& lambda1)
{
  bool  return_value = false;

  if(! SAGE::isnan(lambda1.value))
  {
    if(! lambda1.fixed)
    {
      my_parameters[0].value = lambda1.value;
      my_parameters[0].status = indep_non_func;
    
      if(! (my_parameters[0].lower_bound <= lambda1.value &&
            lambda1.value <= my_parameters[0].upper_bound    ))   
      {
        my_errors << priority(error) << "Value specified for transformation "
                  << "power parameter is not within bounds.  Ignoring ..." << endl;
        
        adjust_lambda_one();
      }
      
      return_value = true; 
    }
    else
    {
      my_parameters[0].value = lambda1.value;
      my_parameters[0].status = lambda1.fixed ? fixed : indep_non_func;
      
      return_value = true;
    }
  }
  else
  {
    // - Value must be supplied if fixed.
    //
    assert(! lambda1.fixed);
    adjust_lambda_one();
  
    return_value = true;
  }
  
  return return_value;
}
   
// - Power parameter limits.
//
bool
transformation_sub_model::set_lambda_one_limits(const model_input& lambda1,
                                                double lambda1_lb, double lambda1_ub)
{
  bool  return_value = false;

  // - Limits not relevant if lambda one is fixed.
  //
  if(lambda1.fixed && (! SAGE::isnan(lambda1_lb) || ! SAGE::isnan(lambda1_ub)))
  {
    my_errors << priority(error) << "Limit(s) specified for fixed transformation "
              << "power parameter.  Ignoring ..." << endl;
    return_value = true;
  }
  else
  {
    if(! SAGE::isnan(lambda1_lb))
    {
      my_parameters[0].lower_bound = lambda1_lb;
    }
    
    if(! SAGE::isnan(lambda1_ub))
    {
      my_parameters[0].upper_bound = lambda1_ub;
    }
    
    // - Check the limits.
    //
    if(my_parameters[0].lower_bound >= my_parameters[0].upper_bound)
    {
      my_errors << priority(critical) << "Conflicting values specified for "
                << "lower and upper limits of transformation power parameter.  "
                << "Skipping analysis ..." << endl;
    }
    else
    {
      return_value = true;
    }
  }

  return return_value;    
} 
   
// - If lambda one is not within bounds, change it so that it is.
//
void
transformation_sub_model::adjust_lambda_one()
{
  assert(my_parameters.size());
  assert(my_parameters[0].status != fixed);
  
  if(! (my_parameters[0].lower_bound <= my_parameters[0].value &&
        my_parameters[0].value <= my_parameters[0].upper_bound    ))
  {
    if(my_parameters[0].upper_bound != POSITIVE_INF)
    {
      my_parameters[0].value = (my_parameters[0].lower_bound + my_parameters[0].upper_bound) / 2;
    }
    else
    {
      if(my_parameters[0].lower_bound != 0)
      {
        my_parameters[0].value = my_parameters[0].lower_bound  + .5 * std::fabs(my_parameters[0].lower_bound);
      }
      else
      {
        my_parameters[0].value = 1;
      }
    }
  }
}
   
// - Shift parameter.
//
bool  
transformation_sub_model::set_lambda_two(const model_input& lambda2)
{
  if(! SAGE::isnan(lambda2.value))
  {
    my_parameters[1].value = lambda2.value;
    my_parameters[1].status = lambda2.fixed ? fixed : indep_non_func;
  }
  
  return true;
}

//
// ==============  Synchronization w. Maxfun
//
int
transformation_sub_model::synchronize(parameter_iterator start)
{
  //lint -e{534}
  sub_model::synchronize(start);
  internally_synchronize();

  if(clear_each_sync)
    my_geometric_mean = QNAN;

  return 0;
}

//
// ===============  Internal synchronization
//
void  
transformation_sub_model::internally_synchronize()
{
  if(my_option != no_trans)
  {
    my_lambda_one = my_parameters[0].value;
    my_lambda_two = my_parameters[1].value;
  }
}

//
// ===============  Transformation
//
// - Given a vector of doubles, transform each member of the vector
//   per D21 of Yi Dong's write-up dated 8/3/01.
//
//   For Box-Cox option trait value + lambda two must be greater
//   than 0 for every trait or return value is false and no trans-
//   formations are performed.
//
bool
transformation_sub_model::transform(std::vector<double>& traits) const
{
  if(my_option == no_trans) return true;

  if(SAGE::isnan(my_geometric_mean))
    return false;

  std::vector<double>::iterator  iter;
  for(iter = traits.begin(); iter != traits.end(); ++iter)
  {
    if(!transform(*iter))
    {
      return false;
    }
  }

  return true;
}

//
// ===============  Transformation
//
// - Given a double, transform it per D21 of Yi Dong's write-up dated
// - 8/3/01.
//
//   For Box-Cox option trait value + lambda two must be greater than 0 for
//   every trait or return value is false and no trans- formations are
//   performed.
//
bool transformation_sub_model::transform(double& trait) const
{
  // If the trait is not finite, there is nothing to do.

  if(!finite(trait)) return true;

  switch(my_option)
  {
    case no_trans:

      break;
      
    case box_cox:

      if((trait + my_lambda_two) < 0.0) return false;

      if(my_lambda_one == 0)
      {
        bc_transform_power_zero(my_geometric_mean, trait);
      }
      else
      {
        bc_transform_power_non_zero(my_geometric_mean, trait);
      }
      
      break;
      
    case george_elston:
      
      if(my_lambda_one == 0)
      {
        ge_transform_power_zero(my_geometric_mean, trait);
      }
      else
      {
        ge_transform_power_non_zero(my_geometric_mean, trait);
      }
      
      break;
      
    default:
      SAGE_internal_error();
  }

  // If, after transform, the trait is not finite, then something is wrong,
  // and we return an error.
  //lint -e{734}
  return finite(trait);
}

// - Returns QNAN if trait + shift parameter not greater than 0 for all traits.
//
double  
transformation_sub_model::bc_geom_mean(const std::vector<double>& traits) const
{
  size_t  adj_trait_count = 0;
  log_double  G(1.0);
  
  std::vector<double>::const_iterator  iter;
  for(iter = traits.begin(); iter != traits.end(); ++iter)
  {
    if(SAGE::isnan(*iter))
    {
      continue;
    }

    // If the value is not nan

    ++adj_trait_count;

    double  shifted_t = *iter + my_lambda_two;

    if(shifted_t <= 0)
    {
      return QNAN;
    }
    else
    {
      G *= shifted_t;
    }
  }

  if(adj_trait_count == 0) return QNAN;

  //lint -e{534}
  G.pow(1.0 / adj_trait_count);
  
  return G.get_double();
}

double  
transformation_sub_model::ge_geom_mean(const std::vector<double>& traits) const
{
  size_t  trait_count = traits.size();
  size_t  adj_trait_count = trait_count;
  log_double  G(1.0);
  
  std::vector<double>::const_iterator  iter;
  for(iter = traits.begin(); iter != traits.end(); ++iter)
  {
    if(SAGE::isnan(*iter))
    {
      --adj_trait_count;
      continue;
    }
    
    double  shifted_t = *iter + my_lambda_two;
    G *= fabs(shifted_t) + 1;
  }

  if(adj_trait_count == 0) return QNAN;

  //lint -e{534}
  G.pow(1.0 / adj_trait_count);
  
  return G.get_double();
}

// These local helper functions help parse a transformation sub block.

void  get_lambda_one_bounds
  (double&        lambda_one_lb,
   double&        lambda_one_ub,
   const LSFBase* param,
   const string&  name_phrase,
   bool           fixed,
   cerrorstream&  errors);

bool  get_lambda_one
  (model_input&   lambda_one,
   double&        lambda_one_lb, 
   double&        lambda_one_ub,
   const LSFBase* param,
   cerrorstream&  errors);

bool  get_lambda_two
  (model_input&   lambda_two,
   const LSFBase* param,
   cerrorstream&  errors);

bool
parse_transformation_sub_model
    (transformation_sub_model&           tsm,
     const LSFBase*                      param,
     transformation_sub_model::sm_option option,
     cerrorstream&                       errors)
{
  static char exp_list[][13] = 
    { "VAL", "FIXED", "VALUE", "INTERACTION", "LOWER_BOUND", "UPPER_BOUND" };

  //bool  option_specified = false;                             
  
  model_input  lambda_one(QNAN, LAMBDA_ONE_DEFAULT_FIXED);
  model_input  lambda_two(QNAN, LAMBDA_TWO_DEFAULT_FIXED);
  double       lambda_one_lb = QNAN;
  double       lambda_one_ub = QNAN;

  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());

      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          string option_string = toUpper(a.String());
          if(option_string == "NONE")
          {
            option = transformation_sub_model::no_trans;
            //option_specified = true;
          }
          else if(option_string == "BOX_COX")
          {
            option = transformation_sub_model::box_cox;
            //option_specified = true;
          }
          else if(option_string == "GEORGE_ELSTON")
          {
            option = transformation_sub_model::george_elston;
            //option_specified = true;
          }
          else
          {
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of transformation sub-block not recognized.  "
                   << "Ignoring ... " << endl;
          }
        }
      }
      else if(param_name == "LAMBDA1")
      {
        APP::check_for_incorrect_attributes(*iter, exp_list, 6, errors);
        
        if(! get_lambda_one(lambda_one, lambda_one_lb, lambda_one_ub, *iter, errors))
        {
          return false;
        }
        
      }
      else if(param_name == "LAMBDA2")
      {
        // Note: Checks only the first two attributes, because the others
        //       don't belong for lambda_two
        APP::check_for_incorrect_attributes(*iter, exp_list, 2, errors);

        if(! get_lambda_two(lambda_two, *iter, errors))
        {
          return false;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in transformation sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "transformation sub-block contains no parameters."
           << endl;
  }
  

  if(! tsm.set(option, lambda_one, lambda_two, lambda_one_lb, lambda_one_ub))
  {
    return false;
  }

  return true;

}

void  get_lambda_one_bounds
  (double&        lambda_one_lb,
   double&        lambda_one_ub,
   const LSFBase* param,
   const string&  name_phrase,
   bool           fixed,
   cerrorstream&  errors)
{
  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    // - Lower bound.
    //
    iter = a_list->find("LOWER_BOUND");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          lambda_one_lb = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for 'lower_bound' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
    
    // - Upper bound.
    //
    iter = a_list->find("UPPER_BOUND");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          lambda_one_ub = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for 'upper_bound' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  }
}

bool get_value
    (model_input&   mi,
     double         def,
     const LSFBase* param,
     const string&  name_phrase,
     cerrorstream&  errors)
{
  double value = QNAN;
  double val   = QNAN;

  APP::LSFConvert::error_t err1 = APP::LSFConvert::to_real(param, value);
  APP::LSFConvert::error_t err2 = APP::LSFConvert::to_real(param, "val", val);

  if(err1 == APP::LSFConvert::GOOD && err2 == APP::LSFConvert::GOOD)
  {
    errors << priority(error) << "Both value and val attribute set for "
           << name_phrase << ".  Will ignore value." << endl;
    err1 = APP::LSFConvert::INVALID;
  }

  switch(err1)
  {
    case APP::LSFConvert::GOOD          :
      mi.value = value;
      break;

    case APP::LSFConvert::ATTR_INVALID  :
      errors << priority(error) << "Value of " << name_phrase
             << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }

  switch(err2)
  {
    case APP::LSFConvert::GOOD          :
      mi.value = val;
      break;

    case APP::LSFConvert::ATTR_INVALID  :
      errors << priority(error) << "'val' of " << name_phrase
             << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }

  // If we haven't got a value by now, we can use the default.
  // If the default is nan, there's no real change.   
 
  if(SAGE::isnan(mi.value))
    mi.value = def;

  // Parse fixed attribute

  bool fixed;

  APP::LSFConvert::error_t fixed_err = APP::LSFConvert::to_boolean(param, "FIXED", fixed);

  switch(fixed_err)
  {
    case APP::LSFConvert::GOOD :
      mi.fixed = fixed;

      // - User specified value for fixed, but not for val.  Use default
      //   for val.  If left as QNAN, sub-model will ignore entire model_input.
      if(fixed && SAGE::isnan(mi.value))
      {
        errors << priority(critical) << "No value given for " << name_phrase
               << " with attribute, fixed, equal to 'true'.  Skipping analysis ..." << endl;
        return false;
      }
      break;

    case APP::LSFConvert::ATTR_INVALID :
      errors << priority(error) << "Value of 'fixed' attribute of " 
             << name_phrase << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }
  
  return true;
}


bool  get_lambda_one
  (model_input&   lambda_one,
   double&        lambda_one_lb, 
   double&        lambda_one_ub,
   const LSFBase* param,
   cerrorstream&  errors)
{
  string  name_phrase = "lambda_one in transformation sub-block";
  
  if(! get_value(lambda_one, QNAN, param, name_phrase, errors))
  {
    return false;
  }

  get_lambda_one_bounds(lambda_one_lb, lambda_one_ub, param, name_phrase, 
                        lambda_one.fixed, errors);
                        
  return true;
}


bool  get_lambda_two
  (model_input&   lambda_two,
   const LSFBase* param,
   cerrorstream&  errors)
{
  string  name_phrase = "lambda_two in transformation sub-block";
  
  // - fixed default available for this parameter.
  //
  if(! get_value(lambda_two, LAMBDA_TWO_DEFAULT_VALUE, param, name_phrase, errors))
  {
    return false;
  }
  
  return true;
}

}
