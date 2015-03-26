//============================================================================
// File:      freq_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/9/01 - created.                            djb
//                                                                          
// Notes:     implementation of genotype_frequency_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/freq_sub_model.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{

const std::string  GENO_FREQ_NAME        = "Type frequency";
const double       PROB_AA_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // .25;
const double       PROB_AB_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // .5;
const double       PROB_BB_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // .25;
const double       PROB_AA_DEFAULT_VALUE_NONE = .25;
const double       PROB_AB_DEFAULT_VALUE_NONE = .5; 
const double       PROB_BB_DEFAULT_VALUE_NONE = .25;
const double       PROB_TOTAL = 1;
const bool         PROBS_FIXED_DEFAULT = false;
const double       FREQ_EPSILON = .00001;           
const double       FREQ_A_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN(); // .5;
const double       FREQ_A_DEFAULT_VALUE_NONE = .5;  
const double       FREQ_A_LB = 0;        
const double       FREQ_A_UB = 1;
const double       FREQ_A_EPSILON = numeric_limits<double>::epsilon();      
const double       PROB_LB = 0;        
const double       PROB_UB = 1;        
const double       CORR_DEFAULT_VALUE = 0;
const bool         CORR_DEFAULT_FIXED = false;           
const double       CORR_LB = -1;
const double       CORR_UB = 1;

//============================================================================
// IMPLEMENTATION:  genotype_frequency_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
// - Call the appropriate sub-model setting function for the 
//   specified sub-model option.  A return value of false means
//   that sub-model set failed.
//
//   Does not take the state of the mean sub-model into account.
//
bool
genotype_frequency_sub_model::set
      (sm_option opt, double freq_a,
       double prob_AA, double prob_AB,
       double prob_BB, const model_input& gcorr, bool probs_fixed)
{
  if(isLinked())
  {
    return false;
  }
 
  bool return_value = false;
  switch(opt)
  {
    case hwe:
    
      if(! freq_A_input_meets_constraints(freq_a, probs_fixed))
      {
        return false;
      }
    
      return_value = set_hwe(freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
      break;
    case nhwe:
      return_value = set_nhwe(freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
      break;
    case NONE:
      return_value = set_none();
      break;
    default:
      assert(false);  
  }
  
  return return_value;
}

// - Call the appropriate sub-model setting function for the 
//   specified sub-model option.  A return value of false means
//   that sub-model set failed.
//
//   Assumes mean sub-model is already set.
//
bool
genotype_frequency_sub_model::set
      (sm_option opt, double freq_a,
       double prob_AA, double prob_AB,
       double prob_BB, const model_input& gcorr, 
       const genotype_specific_mean_susc_sub_model& type_sm, 
       bool type_missing, bool trans_missing, bool probs_fixed)
{
  if(isLinked())
  {
    return false;
  }

  // - Not specifying mean sub-model is a way of specifying commingling
  //   analysis, but certain constraints apply.  See also set functions
  //   in variance and transmission frequency sub-models.
  //
  if(type_missing)
  {
    if(opt == hwe)
    {
      return set(opt, freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
    }
    else
    {
      my_errors << priority(error) << "nhwe option specified in geno_freq sub-block "
                << "for commingling analysis.  Ignoring ..." << endl;
      return true;
    }
  }
  
  // - Specifying mean (or suscept) sub-model and not specifying transmission model 
  //   is a way of specifying transmission analysis, but certain constraints 
  //   apply.  See also set function in the type_var sub-model.
  //
  if(trans_missing && type_sm.option () != genotype_specific_mean_sub_model::one && opt != hwe)
  {
    my_errors << priority(error) << "geno_freq sub-block specified with "
              << "mean or suscept sub-block specified and transmission sub-block not "
              << "specified.  Ignoring geno_freq sub-block ..." << endl;
    return true;
  }
  
  if(type_sm.option() != genotype_specific_mean_sub_model::one)
  {
    // - SCR 190. nhwe option not allowed with two types.
    //
    if(opt == nhwe && type_sm.is_two_option())
    {
      my_errors << priority(error) << "nhwe option not allowed in geno_freq "
                << "sub-block with two types in type_mean or type_suscept sub-block.  Using "
                << "hwe option ..." << endl;
      return set(hwe, freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
    }
    else
    {
      return set(opt, freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
    }
  }
  else
  {
    my_errors << priority(error) << "geno_freq sub-block not "
                << "relevant for one type.  Ignoring ..."  << endl; 
    return set_none();
  }
}

bool
genotype_frequency_sub_model::set_none()
{
  my_option = NONE;
  my_freq_A = FREQ_A_DEFAULT_VALUE_NONE;
  my_probs[index_AA] = PROB_AA_DEFAULT_VALUE_NONE;    
  my_probs[index_AB] = PROB_AB_DEFAULT_VALUE_NONE;    
  my_probs[index_BB] = PROB_BB_DEFAULT_VALUE_NONE;    
  
  initialize_none();
  
  return true;
}

int
genotype_frequency_sub_model::finalizeConfiguration()
{
  switch(my_option)
  {
    case NONE :
      initialize_none();
      break;
    case hwe  :
      initialize_hwe(my_parameters[index_AA].initial_type == MAXFUN::Parameter::FIXED);
      break;
    case nhwe  :
      initialize_nhwe(my_parameters[index_AA].initial_type == MAXFUN::Parameter::FIXED);
      break;
    default:
      SAGE_internal_error();
      break;
  }

  return 0;
}

void
genotype_frequency_sub_model::initialize_none()
{
  my_parameters.resize(0);
}

// - Only the frequency of allele A is considered for this option.  All
//   probabilities are dependent on it.
//
bool  
genotype_frequency_sub_model::set_hwe(double freq_a, double prob_AA, 
                                      double prob_AB, double prob_BB,
                                      const model_input& gcorr, bool probs_fixed)
{
  // Initialize our parameters to QNAN

  my_freq_A = QNAN;

  my_probs[index_AA] = QNAN;
  my_probs[index_BB] = QNAN;
  my_probs[index_AB] = QNAN;

  // - No fixed default values provided per rce request of 8-28-01.
  //
  if(probs_fixed && SAGE::isnan(freq_a))
  {
    my_errors << priority(critical) << "freq_A not specified "
              << "in geno_freq sub-block with probs_fixed equal to true "
              << "and hwe option.  Skipping analysis ..." << endl;
    return false;
  }

  // - Set default values for hwe option.
  //
  my_option = hwe;
  
  bool probs_ok = false;

  // - Genotype probabilities.
  //
  if(! SAGE::isnan(freq_a))
  {
    if(PROB_LB <= freq_a && freq_a <= PROB_UB)
    {
      // We have a freq_A within bounds, so we know we're ok and can compute 
      // our initial values.
      
      my_default = false;
      probs_ok   = true; 

      // Compute initial values
      
      my_freq_A = freq_a;

      my_probs[index_AA] = freq_a * freq_a;
      my_probs[index_BB] = (1 - freq_a) * (1 - freq_a);
      my_probs[index_AB] = PROB_TOTAL - (my_probs[index_AA] + my_probs[index_BB]);
    }
    else
    {
      if(probs_fixed)
      {
        my_errors << priority(critical) << "Value specified for frequency "
                  << "of allele A outside of allowable range.  "
                  << "Skipping analysis ..." << endl;
      }
      else
      {
        my_errors << priority(error) << "Value specified for frequency "
                  << "of allele A outside of allowable range.  "
                  << "Ignoring ..." << endl;
        probs_ok = true;
      }
    }
  }
  else
  {
    probs_ok = true;
  }

  if(! (SAGE::isnan(prob_AA) && SAGE::isnan(prob_AB) && SAGE::isnan(prob_BB)))
  {
    my_errors << priority(warning) << "type_prob parameter(s) specified in geno_freq sub-block "
              << "with option, hwe.  Ignoring ..." << endl;
  }
  
  if(probs_ok)
    initialize_hwe(probs_fixed);
  
  return probs_ok;
}
void  
genotype_frequency_sub_model::initialize_hwe(bool probs_fixed)
{
  // - Set default values for hwe option.
  //
  my_parameters.resize(NUM_OF_PROBS + 1);
  
  my_parameters[index_AA] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_AA",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT, 
            my_probs[index_AA],
            PROB_LB,
            PROB_UB);
  my_parameters[index_AB] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_AB",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT, 
            my_probs[index_AB],
            PROB_LB,
            PROB_UB);
  my_parameters[index_BB] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_BB",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT, 
            my_probs[index_BB],
            PROB_LB,
            PROB_UB);
  my_parameters[index_freq_A] = 
          MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "freq_A",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            my_freq_A,
            PROB_LB + FREQ_A_EPSILON,
            PROB_UB - FREQ_A_EPSILON);
}

// - Only genotype probabilities are considered for this option.  Probabilities
//   AA and BB are indepedent.  Probability AB is dependent on them.  Either none
//   or all probabilities are fixed.
//

bool  
genotype_frequency_sub_model::set_nhwe(double freq_a, double prob_AA, 
                                       double prob_AB, double prob_BB,
                                       const model_input& gcorr, bool probs_fixed)
{
  // - Set default values for nhwe option.
  //
  my_option = nhwe;
  
  bool probs_ok = false;
  
  // - Genotype probabilities.
  //
  genotype_info  info = total_info(prob_AA, prob_AB, prob_BB);
  double  probs_total = prob_AA + prob_AB + prob_BB;

  // - If we have no probabilities, we may use freq_A to determine them.
  //
  if(info == no_geno)
  {
  
    // - Frequency of A is available?
    //
    if(! SAGE::isnan(freq_a))
    {
      if(! freq_A_input_meets_constraints(freq_a, probs_fixed))
      {
        return false;
      }
      
      // - freq_A_input_meets_constraints() may set freq_a to QNAN!
      //   If so, try again.  The 'right' thing to do is probably to break
      //   freq_A_input_meets_constraints() into two or more pieces, but
      //   this works as is.
      //
      if(SAGE::isnan(freq_a))
      {
        return set_nhwe(freq_a, prob_AA, prob_AB, prob_BB, gcorr, probs_fixed);
      }
    
      // - We set everything canonically
      //
      double freq_B = 1.0 - freq_a;

      prob_AA = freq_a * freq_a;
      prob_AB = 2.0 * freq_a * freq_B;
      prob_BB = freq_B * freq_B;

      probs_total = 1.0;

      info = all;
    }
  }
  else
  {
    if(! SAGE::isnan(freq_a))
    {
      my_errors << priority(error) << "Parameters freq_A and prob(s) "
                << "specified in geno_freq sub-block "
                << "with nhwe option.  Ignoring parameter, freq_A ..." << endl;
    }
  }  
  
  switch(info)
  {
    case no_geno:
      
      // - No fixed default values provided per rce request of 8-28-01.
      //
      if(probs_fixed)
      {
        my_errors << priority(critical) << "No genotype probabilities "
                  << "specified in geno_freq sub-block with probs_fixed equal to true "
                  << "and nhwe option.  Skipping analysis ..." << endl;
        return false;
      }
    
      probs_ok = true;
      break;
      
    case AA:      
    case AB:
    case BB:
      probs_ok = set_one_prob(probs_fixed);
      break;
      
    case AA_AB:   
      probs_ok = set_two_probs(prob_AA, index_AA, 
                               prob_AB, index_AB, probs_fixed);
      break;
      
    case AA_BB:
      probs_ok = set_two_probs(prob_AA, index_AA, 
                               prob_BB, index_BB, probs_fixed);
      break;
      
    case AB_BB:
      probs_ok = set_two_probs(prob_AB, index_AB, 
                               prob_BB, index_BB, probs_fixed);
      break;
      
    case all:
      if((PROB_LB <= prob_AA && prob_AA <= PROB_UB) &&
         (PROB_LB <= prob_AB && prob_AB <= PROB_UB) &&
         (PROB_LB <= prob_BB && prob_BB <= PROB_UB) &&
         (PROB_TOTAL - FREQ_EPSILON < probs_total  && 
          probs_total < PROB_TOTAL + FREQ_EPSILON))
      {
        my_default = false;

        my_probs[index_AA] = prob_AA;
        my_probs[index_AB] = PROB_TOTAL - (prob_AA + prob_BB);
        my_probs[index_BB] = prob_BB; 
        my_freq_A = my_probs[index_AA] + .5 * my_probs[index_AB];

        probs_ok = true; 
      }
      else
      {
        if(probs_fixed)
        {
          my_errors << priority(critical) << "Sum of specified genotype probabilities "
                    << "not equal to 1 or individual probability out of range.  "
                    << "Skipping analysis ..." << endl;
        }
        else
        {
          my_errors << priority(error) << "Sum of specified genotype probabilities not equal to 1.  "
                    << "or individual probability out of range.  "
                    << "Ignoring all specified genotype probabilities..." << endl;
          probs_ok = true;            
        }
      }
      
    break;
    
    default:
      assert(false);
  }
  
  initialize_nhwe(probs_fixed);
  
  return probs_ok;
}
 
void  
genotype_frequency_sub_model::initialize_nhwe(bool probs_fixed)
{
  my_parameters.resize(NUM_OF_PROBS + 1);
  
  my_parameters[index_AA] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_AA",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            my_probs[index_AA],
            PROB_LB,
            PROB_UB);
  my_parameters[index_AB] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_AB",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT, 
            my_probs[index_AB],
            PROB_LB,
            PROB_UB);
  my_parameters[index_BB] =
        MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "prob_BB",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 
            my_probs[index_BB],
            PROB_LB,
            PROB_UB);
  my_parameters[index_freq_A] = 
          MAXFUN::ParameterInput
          ( "TYPE FREQUENCIES",
            "freq_A",
            probs_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT, 
            my_freq_A,
            FREQ_A_LB + FREQ_A_EPSILON,
            FREQ_A_UB - FREQ_A_EPSILON);
}

//
// ===============  Ancillary functions
//
// - If fixed is true and input is not w/i bounds (or QNAN), return false.
//   If fixed is false and input value is not w/i bounds (or QNAN), set input to
//   QNAN.  In either case, write an appropriate message.  
//
bool
genotype_frequency_sub_model::freq_A_input_meets_constraints(double& input,
                                                             bool fa_fixed)
{
  bool  return_value = false;
  if(SAGE::isnan(input))
  {
    return_value = true;
  }
  else
  {
    if(FREQ_A_LB + FREQ_A_EPSILON <= input && input <= FREQ_A_UB - FREQ_A_EPSILON)
    {
      return_value = true;
    }
    else
    {
      if(fa_fixed)
      {
        my_errors << priority(critical) << "freq_A specified "
                  << "in geno_freq sub-block is not within allowable limits.  "
                  << "Skipping analysis ..." << endl;
      }
      else
      {
        input = QNAN;
        my_errors << priority(error) << "freq_A specified "
                  << "in geno_freq sub-block is not within allowable limits.  "
                  << "Ignoring ..." << endl;
        return_value = true;
      }
    }
  }
  
  return return_value;  
}

// - Set genotype probabilities when user has specified only one.
//
bool  
genotype_frequency_sub_model::set_one_prob(bool probs_fixed)
{
  bool  probs_ok = false;

  if(probs_fixed)
  {
    my_errors << priority(critical) << "Only one genotype probability "
              << "specified in geno_freq sub-block with probs_fixed equal to true "
              << "and nhwe option.  Skipping analysis ..." << endl;
  }
  else
  {
    my_errors << priority(error) << "Only one genotype probability specified.  "
              << "Ignoring ..." << endl;
    probs_ok = true;
  }
  
  return probs_ok;
}

// - Set genotype probabilities when user has specified only two.
//
bool
genotype_frequency_sub_model::set_two_probs
      (double prob_one, genotype_index index_one,
       double prob_two, genotype_index index_two, bool probs_fixed)
{
  bool  probs_ok = false;
  genotype_index  index_three = third_index(index_one, index_two);

  if( (PROB_LB <= prob_one && prob_one <= PROB_UB) &&
      (PROB_LB <= prob_two && prob_two <= PROB_UB)    )
  {
    if(prob_one + prob_two <= PROB_TOTAL)
    {
      my_default = false;

      my_probs[index_one]   = prob_one;
      my_probs[index_two]   = prob_two;
      my_probs[index_three] = PROB_TOTAL - (prob_one + prob_two); 
      my_freq_A = my_probs[index_AA] + .5 * my_probs[index_AB];

      probs_ok = true; 
    }
    else
    {
      if(probs_fixed)  
      {
        my_errors << priority(critical) << "Sum of specified genotype probabilities "
                  << "exceeds 1.  Skipping analysis ..." << endl;
      }
      else
      {
        my_errors << priority(error) << "Sum of specified genotype probabilities exceeds 1.  "
                  << "Ignoring all specified genotype probabilities ..." << endl;
        probs_ok = true;
      }
    }
  }
  else
  {
    if(probs_fixed)
    {
      my_errors << priority(critical) << "Specified genotype probability "
                << "outside of allowable range.  "
                << "Skipping analysis ..." << endl;
    }
    else
    {
      my_errors << priority(error) << "Specified genotype probability "
                << "outside of allowable range.  "
                << "Ignoring ..." << endl;
      probs_ok = true;
    }
  }
  
  return probs_ok;
}

//
// ==============  Synchronization w. Maxfun
//
int
genotype_frequency_sub_model::update()
{
  int  return_value = 1;

  switch(my_option)
  {
    case hwe:
      return_value = synchronize_hwe();
      break;
    case nhwe:
      return_value = synchronize_nhwe();
      break;
    case NONE:
      return_value = 0;
      break;
    default:
      assert(false);
  }
  
  return return_value;
}

int
genotype_frequency_sub_model::synchronize_hwe()
{
  // - Probabilities of AA, AB and BB are all dependent on frequency of A.
  //   Either all are fixed or none are fixed.  If we're not fixed, we must
  //   compute the values from the freq_A value.
  
  if(my_parameters[index_freq_A].initial_type != MAXFUN::Parameter::FIXED)
  {
    // Copy our independent variables from maxfun.
  
    my_freq_A          = getParam(index_freq_A);

    // Compute our new frequencies
    
    double one_minus_fa = 1.0 - my_freq_A;
    
    my_probs[index_AA] = my_freq_A * my_freq_A;
    my_probs[index_BB] = one_minus_fa * one_minus_fa;
    my_probs[index_AB] = PROB_TOTAL - (my_probs[index_AA] + my_probs[index_BB]);

    // Copy the newly computed values back to maxfun
    
    getParam(index_AA) = my_probs[index_AA];
    getParam(index_AB) = my_probs[index_AB];
    getParam(index_BB) = my_probs[index_BB];
  }
  
  return 0;
}

int
genotype_frequency_sub_model::synchronize_nhwe()
{
  // Assume there is an error unless told otherwise.
  
  bool  return_value = 1;

  // - Probabilities AA and BB are either fixed or indep_func.  Probability AB
  //   fixed or dependent.  Frequency of allele A either fixed or dependent.
  //   Either all are fixed or none are fixed.  IF we're not fixed, we must compute
  //   the freq_A and the index_AB frequency

  if(my_parameters[index_AA].initial_type != MAXFUN::Parameter::FIXED)
  {
    // Copy the values from Maxfun which are independent
    
    my_probs[index_AA] = getParam(index_AA);
    my_probs[index_BB] = getParam(index_BB);

    // Compute the total of the independent (AA and AB) frequencies.
    
    double  indep_probs_total = my_probs[index_AA] + my_probs[index_BB];

    // If this total is less than the maximum, we can do our computation.
    // It is an error otherwise.
                                
    if(indep_probs_total <= PROB_TOTAL)
    {
      my_probs[index_AB] = PROB_TOTAL - (my_probs[index_AA] + my_probs[index_BB]);

      my_freq_A = my_probs[index_AA] + 0.5 * my_probs[index_AB];
          
      // Copy the newly computed parameters back into maxfun
          
      getParam(index_AB)     = my_probs[index_AB];
      getParam(index_freq_A) = my_freq_A;

      return_value = 0;
    }
  }
  else
  {
    return_value = 0;
  }
  
  return return_value;
}

}
}
