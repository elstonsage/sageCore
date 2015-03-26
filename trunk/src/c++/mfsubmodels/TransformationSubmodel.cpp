#include "mfsubmodels/TransformationSubmodel.h"

using namespace std;
namespace SAGE        {
namespace MFSUBMODELS {

const std::string  TRANSFORMATION_NAME      = "Transformation";
const double       LAMBDA_ONE_DEFAULT_VALUE =  1;
const bool         LAMBDA_ONE_DEFAULT_FIXED =  false;
const double       LAMBDA_ONE_DEFAULT_LB    = -1.0;
const double       LAMBDA_ONE_DEFAULT_UB    =  MAXFUN::MF_INFINITY;
const double       LAMBDA_TWO_DEFAULT_VALUE =  0;
const bool         LAMBDA_TWO_DEFAULT_FIXED =  true;

//============================================================================
// IMPLEMENTATION:  TransformationSubmodel
//============================================================================
//
// ===============  Setting the sub-model
//
void
TransformationSubmodel::set_option(sm_option opt)
{
     my_option = opt;
}

bool
TransformationSubmodel::set(sm_option opt, const model_input& lambda1,
                              const model_input& lambda2, double lambda1_lb,
                              double lambda1_ub)
{
  clear_each_sync = true;

  my_option = opt;

  if(opt == no_trans)
  {
    return set_none(lambda1, lambda2, lambda1_lb, lambda1_ub);
  }

  // - Set default values.
  //
  my_parameters.resize(2);
  my_parameters[0] = MAXFUN::ParameterInput("Transformation", "lambda_one", MAXFUN::Parameter::INDEPENDENT,
                                             LAMBDA_ONE_DEFAULT_VALUE, 
                                             LAMBDA_ONE_DEFAULT_LB,
                                             LAMBDA_ONE_DEFAULT_UB);

  my_parameters[1] = MAXFUN::ParameterInput("Transformation", "lambda_two", MAXFUN::Parameter::FIXED, 
                                             LAMBDA_TWO_DEFAULT_VALUE, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);

  bool  return_value = set_lambda_one_limits(lambda1, lambda1_lb, lambda1_ub) &&
                       set_lambda_one(lambda1)                                &&
                       set_lambda_two(lambda2);

  my_lambda_one = my_parameters[0].initial_estimate;
  my_lambda_two = my_parameters[1].initial_estimate;
                       
  return return_value;
}

int
TransformationSubmodel::finalizeConfiguration()
{
  if(my_option != no_trans)
  {
    my_parameters[0].initial_estimate = my_lambda_one;
    my_parameters[1].initial_estimate = my_lambda_two;
  }

  return 0;
}

//
// ==============  Ancillary Functions
//
bool  
TransformationSubmodel::set_none(const model_input& lambda1, const model_input& lambda2,
                                   double lambda1_lb, double lambda1_ub)
{
  if(! (SAGE::isnan(lambda1.value) && SAGE::isnan(lambda2.value) &&
        SAGE::isnan(lambda1_lb)    && SAGE::isnan(lambda1_ub)      ))
  {
    my_errors << priority(error) 
              << "Unnecessary parameters specified in "
              << "transformation sub-block with 'none' option.  Ignoring ..." 
              << endl;
  }
    
  my_parameters.resize(0);
    
  return true;
}

// - Power parameter.  Assumes power parameter limits are already set.
//
bool  
TransformationSubmodel::set_lambda_one(const model_input& lambda1)
{
  bool  return_value = false;

  if(! SAGE::isnan(lambda1.value))
  {
    if(! lambda1.fixed)
    {
      my_parameters[0].initial_estimate = lambda1.value;
      my_parameters[0].initial_type     = MAXFUN::Parameter::INDEPENDENT;
    
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
      my_parameters[0].initial_estimate = lambda1.value;
      my_parameters[0].initial_type     = lambda1.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
      
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
TransformationSubmodel::set_lambda_one_limits(const model_input& lambda1,
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
TransformationSubmodel::adjust_lambda_one()
{
  assert(my_parameters.size());
  assert(my_parameters[0].initial_type != MAXFUN::Parameter::FIXED);
  
  if(! (my_parameters[0].lower_bound <= my_parameters[0].initial_estimate &&
        my_parameters[0].initial_estimate <= my_parameters[0].upper_bound    ))
  {
    if(my_parameters[0].upper_bound != POSITIVE_INF)
    {
      my_parameters[0].initial_estimate = (my_parameters[0].lower_bound + my_parameters[0].upper_bound) / 2;
    }
    else
    {
      if(my_parameters[0].lower_bound != 0)
      {
        my_parameters[0].initial_estimate = my_parameters[0].lower_bound  + .5 * std::fabs(my_parameters[0].lower_bound);
      }
      else
      {
        my_parameters[0].initial_estimate = 1;
      }
    }
  }
}
   
// - Shift parameter.
//
bool  
TransformationSubmodel::set_lambda_two(const model_input& lambda2)
{
  if(! SAGE::isnan(lambda2.value))
  {
    my_parameters[1].initial_estimate = lambda2.value;
    my_parameters[1].initial_type     = lambda2.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
  }
  
  return true;
}

int
TransformationSubmodel::update()
{
  if(my_option != no_trans)
  {
    my_lambda_one = getParam(0);
    my_lambda_two = getParam(1);

    if(clear_each_sync)
      my_geometric_mean = QNAN;
  }

  return 0;
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
TransformationSubmodel::transform(std::vector<double>& traits) const
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
bool TransformationSubmodel::transform(double& trait) const
{
  // If the trait is not finite, there is nothing to do.

  if(!finite(trait)) return true;

  switch(my_option)
  {
    case no_trans:

      break;
      
    case box_cox:

      if((trait + lambda_two()) < 0.0) return false;

      if(lambda_one() == 0)
      {
        bc_transform_power_zero(my_geometric_mean, trait);
      }
      else
      {
        bc_transform_power_non_zero(my_geometric_mean, trait);
      }
      
      break;
      
    case george_elston:
      
      if(lambda_one() == 0)
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
TransformationSubmodel::bc_geom_mean(const std::vector<double>& traits) const
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

    double  shifted_t = *iter + lambda_two();

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
TransformationSubmodel::ge_geom_mean(const std::vector<double>& traits) const
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
    
    double  shifted_t = *iter + lambda_two();
    G *= fabs(shifted_t) + 1;
  }

  if(adj_trait_count == 0) return QNAN;

  //lint -e{534}
  G.pow(1.0 / adj_trait_count);
  
  return G.get_double();
}

//============================================================================
// IMPLEMENTATION:  NewTransformationSubmodel
//============================================================================
//
// ===============  Setting the sub-model
//
void
NewTransformationSubmodel::set_option(sm_option opt)
{
     my_option = opt;
}

bool
NewTransformationSubmodel::set(sm_option opt, const model_input& lambda1,
                              const model_input& lambda2, double lambda1_lb,
                              double lambda1_ub)
{
  clear_each_sync = true;

  my_option = opt;

  if(opt == no_trans)
  {
    return set_none(lambda1, lambda2, lambda1_lb, lambda1_ub);
  }

  // - Set default values.
  //
  my_parameters.resize(2);
  my_parameters[0] = MAXFUN::ParameterInput("Transformation", "lambda_one", MAXFUN::Parameter::INDEPENDENT,
                                             LAMBDA_ONE_DEFAULT_VALUE, 
                                             LAMBDA_ONE_DEFAULT_LB,
                                             LAMBDA_ONE_DEFAULT_UB);

  my_parameters[1] = MAXFUN::ParameterInput("Transformation", "lambda_two", MAXFUN::Parameter::FIXED, 
                                             LAMBDA_TWO_DEFAULT_VALUE, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);

  bool  return_value = set_lambda_one_limits(lambda1, lambda1_lb, lambda1_ub) &&
                       set_lambda_one(lambda1)                                &&
                       set_lambda_two(lambda2);

  my_lambda_one = my_parameters[0].initial_estimate;
  my_lambda_two = my_parameters[1].initial_estimate;
                       
  return return_value;
}

int
NewTransformationSubmodel::finalizeConfiguration()
{
  if(my_option != no_trans)
  {
    my_parameters[0].initial_estimate = my_lambda_one;
    my_parameters[1].initial_estimate = my_lambda_two;
  }

  return 0;
}

//
// ==============  Ancillary Functions
//
bool  
NewTransformationSubmodel::set_none(const model_input& lambda1, const model_input& lambda2,
                                   double lambda1_lb, double lambda1_ub)
{
  if(! (SAGE::isnan(lambda1.value) && SAGE::isnan(lambda2.value) &&
        SAGE::isnan(lambda1_lb)    && SAGE::isnan(lambda1_ub)      ))
  {
    my_errors << priority(error) 
              << "Unnecessary parameters specified in "
              << "transformation sub-block with 'none' option.  Ignoring ..." 
              << endl;
  }
    
  my_parameters.resize(0);
    
  return true;
}

// - Power parameter.  Assumes power parameter limits are already set.
//
bool  
NewTransformationSubmodel::set_lambda_one(const model_input& lambda1)
{
  bool  return_value = false;

  if(! SAGE::isnan(lambda1.value))
  {
    if(! lambda1.fixed)
    {
      my_parameters[0].initial_estimate = lambda1.value;
      my_parameters[0].initial_type     = MAXFUN::Parameter::INDEPENDENT;
    
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
      my_parameters[0].initial_estimate = lambda1.value;
      my_parameters[0].initial_type     = lambda1.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
      
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
NewTransformationSubmodel::set_lambda_one_limits(const model_input& lambda1,
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
NewTransformationSubmodel::adjust_lambda_one()
{
  assert(my_parameters.size());
  assert(my_parameters[0].initial_type != MAXFUN::Parameter::FIXED);
  
  if(! (my_parameters[0].lower_bound <= my_parameters[0].initial_estimate &&
        my_parameters[0].initial_estimate <= my_parameters[0].upper_bound    ))
  {
    if(my_parameters[0].upper_bound != POSITIVE_INF)
    {
      my_parameters[0].initial_estimate = (my_parameters[0].lower_bound + my_parameters[0].upper_bound) / 2;
    }
    else
    {
      if(my_parameters[0].lower_bound != 0)
      {
        my_parameters[0].initial_estimate = my_parameters[0].lower_bound  + .5 * std::fabs(my_parameters[0].lower_bound);
      }
      else
      {
        my_parameters[0].initial_estimate = 1;
      }
    }
  }
}
   
// - Shift parameter.
//
bool  
NewTransformationSubmodel::set_lambda_two(const model_input& lambda2)
{
  if(! SAGE::isnan(lambda2.value))
  {
    my_parameters[1].initial_estimate = lambda2.value;
    my_parameters[1].initial_type     = lambda2.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
  }
  
  return true;
}

int
NewTransformationSubmodel::update()
{
  if(my_option != no_trans)
  {
    my_lambda_one = getParam(0);
    my_lambda_two = getParam(1);

    if(clear_each_sync)
      my_geometric_mean = QNAN;
  }

  return 0;
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
NewTransformationSubmodel::transform(std::vector<double>& traits) const
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
bool NewTransformationSubmodel::transform(double& trait) const
{
  // If the trait is not finite, there is nothing to do.

  if(!finite(trait)) return true;

  switch(my_option)
  {
    case no_trans:

      break;
      
    case box_cox:

      if((trait + lambda_two()) < 0.0) return false;

      if(lambda_one() == 0)
      {
        bc_transform_power_zero(my_geometric_mean, trait);
      }
      else
      {
        bc_transform_power_non_zero(my_geometric_mean, trait);
      }
      
      break;
      
    case george_elston:
      
      if(lambda_one() == 0)
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
NewTransformationSubmodel::bc_geom_mean(const std::vector<double>& traits) const
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

    double  shifted_t = *iter + lambda_two();

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
NewTransformationSubmodel::ge_geom_mean(const std::vector<double>& traits) const
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
    
    double  shifted_t = *iter + lambda_two();
    G *= fabs(shifted_t) + 1;
  }

  if(adj_trait_count == 0) return QNAN;

  //lint -e{534}
  G.pow(1.0 / adj_trait_count);
  
  return G.get_double();
}



} // End namespace MFSUBMODELS
} // End namespace SAGE

