//============================================================================
// File:      resid_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/13/01 - created.                            djb
//                                                                          
// Notes:     implementation of residual_correlation_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/resid_sub_model.h"

using namespace std;
namespace SAGE
{

namespace SEGREG
{

const std::string  RESID_NAME            = "Residual correlations";
const double       RESID_SP_DEFAULT_VALUE = 0;
const double       RESID_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // 0;
const double       RESID_EPSILON = .00000001;
const bool         RESID_DEFAULT_FIXED = false;           
const double       RESID_LB = -1;                       // Per rce 8/14/01.
const double       RESID_UB = 1;

//============================================================================
// IMPLEMENTATION:  residual_correlation_sub_model
//============================================================================
//
bool  
residual_correlation_sub_model::set_as_default(model_class m_class)
{
  set(residual_correlation_sub_model::equal_po_ss,
      model_input(),model_input(),model_input(),
      model_input(),m_class);
      
  my_is_in_default_mode = true;
  
  return false;
}

bool  
residual_correlation_sub_model::set
    (sm_option   opt, 
     model_input fm_val,
     model_input mo, 
     model_input fo,
     model_input sib_sib,
     model_class m_class)
{
  if(isLinked())
  {
    return false;
  }

  my_model_class = m_class;
  
  my_is_in_default_mode = false;

  // if the model is continuous, the elements are limited to the (-1, 1)
  // range.

  if(my_model_class != model_MLM)
  {
    if(! (input_meets_constraints(fm_val, RESID_EPSILON) &&
          input_meets_constraints(mo)                &&
          input_meets_constraints(fo)                &&
          input_meets_constraints(sib_sib)              ))
    {
      return false;
    }
  }

  bool valid = false;

  switch(opt)
  {
    case equal_po_ss:
      valid = set_equal_po_ss(fm_val, mo, fo, sib_sib);
      break;
    case equal_po:
      valid = set_equal_po(fm_val, mo, fo, sib_sib);
      break;
    case arb:
      valid = set_arb(fm_val, mo, fo, sib_sib);
      break;
    default:
      assert(false);  
  }

  supply_missing_values();
    
  initialize();

  if(valid)
  {

    calculate_alpha_and_delta();
  }
  
  return valid;
}

int
residual_correlation_sub_model::finalizeConfiguration()
{
  if(has_residuals())
    initialize();

  return 0;
}

void 
residual_correlation_sub_model::initialize()
{
  switch(my_option)
  {
    case equal_po_ss:
      initialize_equal_po_ss();
      break;
    case equal_po:
      initialize_equal_po();
      break;
    case arb:
      initialize_arb();
      break;
    default:
      assert(false);  
  }
}


// - Fill in any parameter values that are QNAN.
//
void
residual_correlation_sub_model::supply_missing_values()
{
    switch(my_model_class)
    {
      case model_A:
        supply_missing_values(0.5, 0.5);
        break;
      
      case model_D:
        supply_missing_values(0.5, 0.98);
        break;

      case model_MLM:
        supply_missing_values(0.3, 0.3);
        break;

      case model_FPMM:
      case model_INVALID:
      default:
        ;
    }
}

inline string get_resid_name(const model_class& mc)
{
  if(mc == model_MLM) return "assoc";
  else                return "corr";
}

bool  
residual_correlation_sub_model::set_equal_po_ss(const model_input& fm_val, const model_input& mo, 
                                              const model_input& fo, const model_input& sib_sib)
{
  my_option = equal_po_ss;

  // Set our two correlation parameters to default arguments
  
  my_corrs       [fm] = RESID_SP_DEFAULT_VALUE;
  my_corrs_fixed [fm] = true;

  my_corrs       [fs] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [fs] = false;
                                           
  bool  return_value = false;
  
  // If fm correlation defined, reset it.
  
  if(! SAGE::isnan(fm_val.value))
  {
    my_corrs      [fm] = fm_val.value;
    my_corrs_fixed[fm] = fm_val.fixed;
  }

  // Define the po, ss correlation.

  return_value = set_po_ss(mo, fo, sib_sib);
  
  if(return_value)
  {
    // If we're valid, copy our state information to the other corrs and build our
    // initialization parameters
    
    initialize_equal_po_ss();
  }
  
  return return_value;
}

void  
residual_correlation_sub_model::initialize_equal_po_ss()
{
  // First, copy our correlations
  
  my_corrs[ms] = 
  my_corrs[md] =
  my_corrs[fd] =
  my_corrs[bb] =
  my_corrs[ss] =
  my_corrs[bs] = my_corrs[fs];

  my_corrs_fixed[ms] = 
  my_corrs_fixed[md] =
  my_corrs_fixed[fd] =
  my_corrs_fixed[bb] =
  my_corrs_fixed[ss] =
  my_corrs_fixed[bs] = my_corrs_fixed[fs];

  // Now setup maxfun parameters.

  string rn = get_resid_name(my_model_class);

  my_parameters.resize(2);
  
  my_parameters[eps_fm] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "marital resid " + rn,
       my_corrs_fixed[fm] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[fm],
       get_lower_bound() + RESID_EPSILON,
       get_upper_bound() - RESID_EPSILON);
  my_parameters[eps_po_ss] =
    MAXFUN::ParameterInput(
      "RESIDUALS",
      "po = ss resid " + rn,
      my_corrs_fixed[fs] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_corrs[fs],
      get_lower_bound(),
      get_upper_bound());
}


bool  
residual_correlation_sub_model::set_equal_po(const model_input& fm_val, const model_input& mo, 
                                              const model_input& fo, const model_input& sib_sib)
{
  my_option = equal_po;

  // Set our three correlation parameters to default arguments

  my_corrs       [fm] = RESID_SP_DEFAULT_VALUE;
  my_corrs_fixed [fm] = true;

  my_corrs       [fs] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [fs] = false;

  my_corrs       [ss] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [ss] = false;

  bool  return_value = false;
  
  // Fix fm correlation, if given.

  if(! SAGE::isnan(fm_val.value))
  {
    my_corrs       [fm] = fm_val.value;
    my_corrs_fixed [fm] = fm_val.fixed;
  }
  
  // Fix po correlation, if given.

  bool  po_specified = false;

  // Check and set the mo correlation
  
  if(! SAGE::isnan(mo.value))
  {
    // Note: we use the fs correlation arbitrarily throughout this process,
    //       even though we're looking at mother offspring right now.  It
    //       is copied to all the po options at the end assuming everything is
    //       valid.
    
    my_corrs       [fs] = mo.value;
    my_corrs_fixed [fs] = mo.fixed;

    po_specified = true;
    return_value = true;
  }
  
  // Check the fo correlation.  Note that only one of mo and fo should be set,
  // and if both are, various errors may result.

  if(! SAGE::isnan(fo.value))
  {
    if(! po_specified)
    {
      // Only the fo is set, so we can set the values happily.

      my_corrs       [fs] = fo.value;
      my_corrs_fixed [fs] = fo.fixed;

      po_specified = true;
      return_value = true;
    }
    else
    {
      // We have both the fo and mo set.  This causes problems.
      
      if(fo.value == mo.value && fo.fixed == mo.fixed)
      {
        // If the fo values match the mo values, there is redundant information
        // but it's not actually bad.
        
        /* This message superflous and misleading.  Removed 4/9/2.  -djb
        my_errors << priority(warning) << "Redundant specification.  Parameters "
                  << "mo and fo both specified in resid sub-block "
                  << "with equal_po option." << endl;
        */
        
        ;
      }
      else
      {
        // the fo and mo don't match.  This is a problem.
        
        if(fo.fixed != mo.fixed)
        {
          // The fixed status doesn't match.  This is a critical problem.
          
          my_errors << priority(critical) << "Parameters mo and fo not in agreement "
                    << "in resid sub-block with equal_po option.  Skipping analysis ..." << endl;
          return_value = false;
        }
        else
        {
          // The fixed statuses match, but the values don't.  
          
          if(! fo.fixed)
          {
            // Since they're not fixed, this isn't critical.  We just ignore one
            // of them.
            
            my_errors << priority(error) << "Parameters mo and fo not in agreement "
                      << "in resid sub-block with equal_po option.  Ignoring fo ..." << endl; 
          }
          else
          {
            // If they're both fixed, but don't match, we don't know what to do,
            // so it's a critical error.
            
            my_errors << priority(critical) << "Parameters mo and fo not in agreement "
                      << "in resid sub-block with equal_po option.  Skipping analysis ..." << endl;
            return_value = false;
          }
        }
      }
    }
  }
  
  // - No value specifed.  Using defaults.
  
  if(! po_specified)
  {
    return_value = true;
  }
  
  // Fix the ss correlation, if given.
  
  if(! SAGE::isnan(sib_sib.value))
  {
    my_corrs       [ss] = sib_sib.value;
    my_corrs_fixed [ss] = sib_sib.fixed;
  }
  
  if(return_value)
  {
    initialize_equal_po();
  }
  
  return return_value;
}

void
residual_correlation_sub_model::initialize_equal_po()
{
  // Copy our state information to the other corrs
  
  my_corrs[ms] = 
  my_corrs[md] =
  my_corrs[fd] = my_corrs[fs];

  my_corrs[bb] =
  my_corrs[bs] = my_corrs[ss];

  my_corrs_fixed[ms] = 
  my_corrs_fixed[md] =
  my_corrs_fixed[fd] = my_corrs_fixed[fs];

  my_corrs_fixed[bb] =
  my_corrs_fixed[bs] = my_corrs_fixed[ss];

  // Build our initialization parameters

  string rn = get_resid_name(my_model_class);

  my_parameters.resize(3);
  
  my_parameters[ep_fm] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "marital resid " + rn,
       my_corrs_fixed[fm] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[fm],
       get_lower_bound() + RESID_EPSILON,
       get_upper_bound() - RESID_EPSILON);

  my_parameters[ep_po] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "po resid " + rn,
       my_corrs_fixed[fs] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[fs],
       get_lower_bound(),
       get_upper_bound());

  my_parameters[ep_ss] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "ss resid " + rn,
       my_corrs_fixed[ss] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[ss],
       get_lower_bound(),
       get_upper_bound());
}

bool  
residual_correlation_sub_model::set_arb(const model_input& fm_val, const model_input& mo, 
                                        const model_input& fo, const model_input& sib_sib)
{
  my_option = arb;

  // Set our initial four correlations.
  
  my_corrs       [fm] = RESID_SP_DEFAULT_VALUE;
  my_corrs_fixed [fm] = true;

  my_corrs       [fs] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [fs] = false;

  my_corrs       [ms] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [ms] = false;

  my_corrs       [ss] = RESID_DEFAULT_VALUE;
  my_corrs_fixed [ss] = true;

  // Update fm correlation if fixed.
                                           
  if(! SAGE::isnan(fm_val.value))
  {
    my_corrs       [fm] = fm_val.value;
    my_corrs_fixed [fm] = fm_val.fixed;
  }
  
  // Update fo correlation if fixed.
                                           
  if(! SAGE::isnan(fo.value))
  {
    my_corrs       [fs] = fo.value;
    my_corrs_fixed [fs] = fo.fixed;
  }
  
  // Update mo correlation if fixed.
                                           
  if(! SAGE::isnan(mo.value))
  {
    my_corrs       [ms] = mo.value;
    my_corrs_fixed [ms] = mo.fixed;
  }
  
  // - Sib - sib correlation.
  //
  if(! SAGE::isnan(sib_sib.value))
  {
    my_corrs       [ss] = sib_sib.value;
    my_corrs_fixed [ss] = sib_sib.fixed;
  }
  
  // Initialize our maxfun parameters
  
  initialize_arb();
  
  return true;
}
    
void
residual_correlation_sub_model::initialize_arb()
{
  // Set all the other correlations

  my_corrs[fd] = my_corrs[fs];

  my_corrs[md] = my_corrs[ms]; 

  my_corrs[bb] =
  my_corrs[bs] = my_corrs[ss];

  my_corrs_fixed[fd] = my_corrs_fixed[fs];

  my_corrs_fixed[md] = my_corrs_fixed[ms]; 

  my_corrs_fixed[bb] =
  my_corrs_fixed[bs] = my_corrs_fixed[ss];
  string rn = get_resid_name(my_model_class);

  // Initial our maximization parameters
  
  my_parameters.resize(4);
  
  my_parameters[arb_fm] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "marital resid " + rn,
       my_corrs_fixed[fm] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[fm],
       get_lower_bound() + RESID_EPSILON,
       get_upper_bound() - RESID_EPSILON);

  my_parameters[arb_fo] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "fo resid " + rn,
       my_corrs_fixed[fs] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[fs],
       get_lower_bound(),
       get_upper_bound());

  my_parameters[arb_mo] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "mo resid " + rn,
       my_corrs_fixed[ms] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[ms],
       get_lower_bound(),
       get_upper_bound());

  my_parameters[arb_ss] =
     MAXFUN::ParameterInput(
       "RESIDUALS",
       "ss resid " + rn,
       my_corrs_fixed[ss] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
       my_corrs[ss],
       get_lower_bound(),
       get_upper_bound());
}
//
// ===============  Ancillary Functions
//
// - If a fixed input value is not w/i bounds (or QNAN), return false.
//   If a non-fixed input value is not w/i bounds (or QNAN), set value to
//   QNAN.  In either case, write an appropriate message.  An optional epsilon
//   argument is applied to the bounds.
//
bool
residual_correlation_sub_model::input_meets_constraints(model_input& input, double epsilon)
{
  bool  return_value = false;
  if(SAGE::isnan(input.value))
  {
    return_value = true;
  }
  else
  {
    if((get_lower_bound() + epsilon) <= input.value && input.value <= (get_upper_bound() - epsilon))
    {
      return_value = true;
    }
    else
    {
      if(input.fixed == true)
      {
        my_errors << priority(critical) << "A correlation specified "
                  << "in resid sub-block is not within allowable limits.  "
                  << "Skipping analysis ..." << endl;
      }
      else
      {
        input.value = QNAN;
        my_errors << priority(error) << "A correlation specified "
                  << "in resid sub-block is not within allowable limits.  "
                  << "Ignoring ..." << endl;
        return_value = true;
      }
    }
  }
  
  return return_value;  
}

bool
residual_correlation_sub_model::set_po_ss(const model_input& mo, const model_input& fo,
                                          const model_input& sib_sib)
{
  bool  po_ss_specified = false;
  bool  return_value = false;

  // Fill in the blanks based on mo, if specified.
  
  if(! SAGE::isnan(mo.value))
  {
    my_corrs       [fs] = mo.value;
    my_corrs_fixed [fs] = mo.fixed;

    po_ss_specified = true;
    return_value    = true;
  }

  // Fill in based on fo.  Note that if both mo and fo are specified, there are
  // problems.
  
  if(! SAGE::isnan(fo.value))
  {
    if(! po_ss_specified)
    {
      my_corrs       [fs] = fo.value;
      my_corrs_fixed [fs] = fo.fixed;
  
      po_ss_specified = true;
      return_value    = true;
    }
    else
    {
      //lint --e{777}
      if(fo.value == mo.value && fo.fixed == mo.fixed)
      {
        // We're ok, since fo matches mo.
        
        /* This message superflous and misleading.  Removed 4/9/2.  -djb
        my_errors << priority(warning) << "Redundant specification.  Parameters "
                  << "mo and fo both specified in resid sub-block "
                  << "with equal_po_ss option." << endl;
        */
        
        ;
      }
      else
      {
        if(fo.fixed != mo.fixed)
        {
          my_errors << priority(critical) << "Parameters mo and fo not in agreement "
                    << "in resid sub-block with equal_po_ss option.  Skipping analysis ..." 
                    << endl;
          return_value = false;
        }
        else
        {
          if(! fo.fixed)
          {
            my_errors << priority(error) << "Parameters mo and fo not in agreement "
                      << "in resid sub-block with equal_po_ss option.  Ignoring fo ..." << endl; 
          }
          else
          {
            my_errors << priority(critical) << "Parameters mo and fo not in agreement "
                      << "in resid sub-block with equal_po_ss option.  Skipping analysis ..." 
                      << endl;
            return_value = false;
          }
        }
      }
    }
  }
  
  if(! SAGE::isnan(sib_sib.value))
  {
    if(! po_ss_specified)
    {
      my_corrs       [fs] = sib_sib.value;
      my_corrs_fixed [fs] = sib_sib.fixed;
  
      po_ss_specified = true;
      return_value    = true;
    }
    else
    {
      // We have a potential problem.  sib_sib is specified when mo or fo has
      // also been specified (or both).
      
      //lint --e{777}
      if(sib_sib.value == my_corrs[fs] && 
         sib_sib.fixed == my_corrs_fixed[fs])
      {
        // We're ok, because everything matches.
        
        /* This message superflous and misleading.  Removed 4/9/2.  -djb
        my_errors << priority(warning) << "Redundant specification.  More "
                  << "than one of the parameters, mo, fo and ss specified in resid sub-block "
                  << "with equal_po_ss option." << endl;
        */
        
        ;
      }
      else
      {
        // Something doesn't match
        
        if(sib_sib.fixed != my_corrs_fixed[fs])
        {
          // Critical problem.  Fixed flags are different.
          my_errors << priority(critical) << "Parameter ss  not in agreement with mo or fo "
                    << "in resid sub-block with equal_po_ss option.  Skipping analysis ..." 
                    << endl;
          return_value = false;
        }
        else
        {
          // Fixed flags are the same, so values must be different.  We can
          // ignore, but only if we're not fixed.  IF we *are* fixed, this
          // is a problem we can't deal with.
          if(! sib_sib.fixed)
          {
            my_errors << priority(error) << "Parameter ss not in agreement with mo or fo "
                      << "in resid sub-block with equal_po_ss option.  Ignoring ss ..." << endl;  
          }
          else
          {
            my_errors << priority(critical) << "Parameter ss  not in agreement with mo or fo "
                      << "in resid sub-block with equal_po_ss option.  Skipping analysis ..." 
                      << endl;
            return_value = false;
          }
        }
      }
    }
  }
  
  // - No value specifed.  Using defaults.
  //
  else if(! po_ss_specified)
  {
    return_value = true;
  }
  
  return return_value;
}

void  
residual_correlation_sub_model::supply_missing_values(double pov, double ssv)
{
  switch(my_option)
  {
    case equal_po_ss:
      if(SAGE::isnan(my_corrs[fm]))
      {    
        my_corrs[fm] = 0;
      }
      if(SAGE::isnan(my_corrs[fs]))
      {    
        my_corrs[fs] = pov;
      }
      break;
    
    case equal_po:
      if(SAGE::isnan(my_corrs[fm]))
      {    
        my_corrs[fm] = 0;
      }
      if(SAGE::isnan(my_corrs[fs]))
      {    
        my_corrs[fs] = pov;
      }
      if(SAGE::isnan(my_corrs[ss]))
      {    
        my_corrs[ss] = ssv;
      }
      break;
    
    case arb:
      if(SAGE::isnan(my_corrs[fm]))
      {    
        my_corrs[fm] = 0;
      }
      if(SAGE::isnan(my_corrs[fs]))
      {    
        my_corrs[fs] = pov;
      }
      if(SAGE::isnan(my_corrs[ms]))
      {    
        my_corrs[ms] = pov;
      }
      if(SAGE::isnan(my_corrs[ss]))
      {    
        my_corrs[ss] = ssv;
      }
      break;
    
    default:
      assert(false);
  }
}

//
// ==============  Synchronization w. Maxfun
//    
int 
residual_correlation_sub_model::update()
{
  int  return_value = 1;

  switch(my_option)
  {
    case equal_po_ss:
      return_value = synchronize_equal_po_ss();
      break;
    case equal_po:
      return_value = synchronize_equal_po();
      break;
    case arb:
      return_value = synchronize_arb();
      break;
    default:
      assert(false);
  }
  
  calculate_alpha_and_delta();

  return return_value;
}

// - For this option, spouse corr. is fixed at 0 and the other corr. is fixed
//   or independent non_func.
//
bool
residual_correlation_sub_model::synchronize_equal_po_ss()
{
  internally_synchronize_equal_po_ss();
  return 0;
}

// - For this option there are three parameters wh. may be fixed or independent
//   non_func.
//
bool
residual_correlation_sub_model::synchronize_equal_po()
{
  internally_synchronize_equal_po();
  return 0;
}

// - For this option there are four parameters wh. may be fixed or independent
//   non_func.
//
bool
residual_correlation_sub_model::synchronize_arb()
{
  internally_synchronize_arb();
  return 0;
}

//
// ===============  Internal synchronization
//
void
residual_correlation_sub_model::internally_synchronize()
{
  switch(my_option)
  {
    case equal_po_ss:
      internally_synchronize_equal_po_ss();
      break;
      
    case equal_po:
      internally_synchronize_equal_po();
      break;
      
    case arb:
      internally_synchronize_arb();
      break;
      
    default:
      assert(false);
  }
}
    
void  
residual_correlation_sub_model::internally_synchronize_equal_po_ss()
{
  my_corrs[fm] = getParam(eps_fm);
  my_corrs[ms] = getParam(eps_po_ss);
  my_corrs[fs] = getParam(eps_po_ss);
  my_corrs[md] = getParam(eps_po_ss);
  my_corrs[fd] = getParam(eps_po_ss);
  my_corrs[bb] = getParam(eps_po_ss);
  my_corrs[ss] = getParam(eps_po_ss);
  my_corrs[bs] = getParam(eps_po_ss);
}

void  
residual_correlation_sub_model::internally_synchronize_equal_po()
{
  my_corrs[fm] = getParam(ep_fm);
  my_corrs[ms] = getParam(ep_po);
  my_corrs[fs] = getParam(ep_po);
  my_corrs[md] = getParam(ep_po);
  my_corrs[fd] = getParam(ep_po);
  my_corrs[bb] = getParam(ep_ss);
  my_corrs[ss] = getParam(ep_ss);
  my_corrs[bs] = getParam(ep_ss);
}

void  
residual_correlation_sub_model::internally_synchronize_arb()
{
  my_corrs[fm] = getParam(arb_fm);
  my_corrs[ms] = getParam(arb_mo);
  my_corrs[fs] = getParam(arb_fo);
  my_corrs[md] = getParam(arb_mo);
  my_corrs[fd] = getParam(arb_fo);
  my_corrs[bb] = getParam(arb_ss);
  my_corrs[ss] = getParam(arb_ss);
  my_corrs[bs] = getParam(arb_ss);
}

}
}
