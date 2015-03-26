//============================================================================
// File:      mle_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/4/2 - created.                       -djb
//                                                                          
// Notes:     implementation of mle_sub_model class.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/mle_sub_model.h"

using namespace std;

namespace SAGE
{

namespace LODLINK
{

const string GROUP_NAME = "mle sub model";
const string MLE_NAME = "Maximum likelihood estimation";

const string  AVERAGE_THETA_DESC = "average recombination fraction";
const string  MALE_THETA_DESC    = "male recombination fraction";
const string  FEMALE_THETA_DESC  = "female recombination fraction";
const pstatus  THETA_DEFAULT_STATUS = MAXFUN::Parameter::INDEPENDENT;
const double  THETA_INIT_VALUE         = .25;
const double  THETA_LOWER_BOUND        = 0;
const double  THETA_STRICT_UPPER_BOUND = .5;
const double  THETA_UPPER_BOUND        = 1;

const string  ALPHA_DESC = "proportion of families with linkage";
const pstatus  ALPHA_DEFAULT_STATUS = MAXFUN::Parameter::INDEPENDENT;
const double  ALPHA_INIT_VALUE  = .5;
const double  ALPHA_LOWER_BOUND = 0;
const double  ALPHA_UPPER_BOUND = 1;

//============================================================================
// IMPLEMENTATION:  mle_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
void
mle_sub_model::set(bool ss, bool ua)
{
  assert(! isLinked());
  
  sex_specific = ss;
  use_alpha = ua;

  // - Create one or two thetas and zero or one alpha parameters
  //   as function arguments dictate.
  //
  int  param_count = ss ? (ua ? 3 : 2) : (ua ? 2 : 1);
  my_parameters.resize(param_count);
  
  switch(param_count)
  {
    case 1:
      my_parameters[AVERAGE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                      AVERAGE_THETA_DESC,
                                                      THETA_DEFAULT_STATUS,
                                                      THETA_INIT_VALUE,
                                                      THETA_LOWER_BOUND,
                                                      THETA_UPPER_BOUND  );
      break;
      
    case 2:
      if(ss)
      {
        my_parameters[MALE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                     MALE_THETA_DESC,
                                                     THETA_DEFAULT_STATUS,
                                                     THETA_INIT_VALUE,
                                                     THETA_LOWER_BOUND,
                                                     THETA_UPPER_BOUND );
        my_parameters[FEMALE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                       FEMALE_THETA_DESC,
                                                       THETA_DEFAULT_STATUS,
                                                       THETA_INIT_VALUE,
                                                       THETA_LOWER_BOUND,
                                                       THETA_UPPER_BOUND );
      }
      else
      {
        my_parameters[AVERAGE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                        AVERAGE_THETA_DESC,
                                                        THETA_DEFAULT_STATUS,
                                                        THETA_INIT_VALUE,
                                                        THETA_LOWER_BOUND,
                                                        THETA_UPPER_BOUND  );
        my_parameters[ALPHA_ONE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                          ALPHA_DESC,
                                                          ALPHA_DEFAULT_STATUS,
                                                          ALPHA_INIT_VALUE,
                                                          ALPHA_LOWER_BOUND,
                                                          ALPHA_UPPER_BOUND ); 
      }
      break;
      
    case 3:
      my_parameters[MALE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                   MALE_THETA_DESC,
                                                   THETA_DEFAULT_STATUS,
                                                   THETA_INIT_VALUE,
                                                   THETA_LOWER_BOUND,
                                                   THETA_UPPER_BOUND );
      my_parameters[FEMALE] = MAXFUN::ParameterInput(GROUP_NAME,
                                                     FEMALE_THETA_DESC,
                                                     THETA_DEFAULT_STATUS,
                                                     THETA_INIT_VALUE,
                                                     THETA_LOWER_BOUND,
                                                     THETA_UPPER_BOUND );
      my_parameters[ALPHA_TWO] = MAXFUN::ParameterInput(GROUP_NAME,
                                                        ALPHA_DESC,
                                                        ALPHA_DEFAULT_STATUS,
                                                        ALPHA_INIT_VALUE,
                                                        ALPHA_LOWER_BOUND,
                                                        ALPHA_UPPER_BOUND ); 
      break;
      
    default:
      assert(false);
  }

  init();
}

// - This function is responsible for:
//    1. Copying Maxfun parameter values to sub_model variables
//       used by the sub_model interface, internal synchronization.
//    2. Supplying maxfun w. values for dependent variables.
//    3. Checking that any constraints between two or more indep_func
//       variables are met.  Return positive integer if a constraint
//       is violated, 0 otherwise.  There are no cases in this sub_model
//       where this is an issue.
//
int
mle_sub_model::update()
{
  // - Special case: sex_specific recombination fractions constrained
  //   to sum to 1.
  //
  if(sex_specific && my_parameters[FEMALE].initial_type == MAXFUN::Parameter::DEPENDENT)
  {
    getParam(FEMALE) = 1 - getParam(MALE);
  }
  
  internally_synchronize();
  
  return  0;
}

void  
mle_sub_model::init()
{
  switch(my_parameters.size())
  {
    case 1:
      my_average_theta = my_parameters[AVERAGE].initial_estimate;
      my_male_theta = QNAN;
      my_female_theta = QNAN;
      my_alpha = QNAN;
      break;
      
    case 2:
      if(sex_specific)
      {
        my_average_theta = QNAN;
        my_male_theta = my_parameters[MALE].initial_estimate;
        my_female_theta = my_parameters[FEMALE].initial_estimate;
        my_alpha = QNAN;
      }
      else
      {
        my_average_theta = my_parameters[AVERAGE].initial_estimate;
        my_male_theta = QNAN;
        my_female_theta = QNAN;
        my_alpha = my_parameters[ALPHA_ONE].initial_estimate;
      }      
      break;
      
    case 3:
      my_average_theta = QNAN;
      my_male_theta = my_parameters[MALE].initial_estimate;
      my_female_theta = my_parameters[FEMALE].initial_estimate;
      my_alpha = my_parameters[ALPHA_TWO].initial_estimate;
      break;
    
    default:
      assert(false);
  }
}

//
// ===============  Internal synchronization
//
void  
mle_sub_model::internally_synchronize()
{
  assert(getParameterMgr());

  switch(my_parameters.size())
  {
    case 1:
      my_average_theta = getParam(AVERAGE);
      my_male_theta = QNAN;
      my_female_theta = QNAN;
      my_alpha = QNAN;
      break;
      
    case 2:
      if(sex_specific)
      {
        my_average_theta = QNAN;
        my_male_theta = getParam(MALE);
        my_female_theta = getParam(FEMALE);
        my_alpha = QNAN;
      }
      else
      {
        my_average_theta = getParam(AVERAGE);
        my_male_theta = QNAN;
        my_female_theta = QNAN;
        my_alpha = getParam(ALPHA_ONE);
      }      
      break;
      
    case 3:
      my_average_theta = QNAN;
      my_male_theta = getParam(MALE);
      my_female_theta = getParam(FEMALE);
      my_alpha = getParam(ALPHA_TWO);
      break;
    
    default:
      assert(false);
  }
}

}
}
