//============================================================================
// File:      fpmm_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/22/01 - created.                            djb
//                                                                          
// Notes:     implementation of finite_polygenic_mixed_model_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/fpmm_sub_model.h"

using namespace std;

namespace SAGE   
{
namespace SEGREG 
{

const std::string  FPMM_NAME             = "Finite polygenic mixed model";
const double       FPMM_DEFAULT_FREQ = .5;
const size_t       FPMM_DEFAULT_LOCI = 3;
const size_t       FPMM_MAX_LOCI = 6;
const double       FPMM_DEFAULT_VALUE =  numeric_limits<double>::quiet_NaN();  // 1;
const double       FPMM_EPSILON = .00001;
const double       FPMM_LB = 0;
const double       FPMM_UB =  numeric_limits<double>::infinity();
const bool         FPMM_DEFAULT_FIXED = false;

//============================================================================
// IMPLEMENTATION:  finite_polygenic_mixed_model_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
bool
finite_polygenic_mixed_model_sub_model::set(model_input var, double freq, size_t nloci)
{
// go around for no polygenic loci;
bool no_loci_found;

 if (nloci == 0) 
       {
         no_loci_found = true;
         var = FPMM_LB + FPMM_EPSILON;
       } else 
            {
            no_loci_found = false;
            }
                   
  if (no_loci_found) nloci = 1;
  
  if(isLinked())
  {
    return false;
  }

// Original location 
  if(! input_meets_constraints(var))
  {
    return false;
  }

  // - Variance.
  if (! no_loci_found) 
    { // modified from previous code by JA for zero polygenic loci
      if(! SAGE::isnan(var.value))
      {
        my_variance       = var.value;
        my_variance_fixed = var.fixed;
      }
      else
      {
        my_variance       = QNAN;
        my_variance_fixed = false;
      }
    } 
   else 
    { // activated if no polygenic loci present
         my_variance = FPMM_LB + FPMM_EPSILON;
         my_variance_fixed = true;
//         cout << "In fpmm checking variance " << my_variance << endl;
     }      
  // - Frequency.
  //

  my_frequency = .5;
  
  bool  freq_ok = false;

  if (! no_loci_found) 
   { // added by JA , block activates with polygenic loci > 0
  if(! SAGE::isnan(freq)) 
   {
    if(0 < freq && freq < 1)
    {
      my_frequency = freq;
      freq_ok = true;
    }
    else
    {
      my_errors << priority(critical) << "Parameter freq in fpmm sub-block outside of "
                << "allowable range.  Skipping analysis ..." << endl;
    }
  }
  else
  {
    freq_ok = true;
  }
      } else  {freq_ok = true;}
  // - Number of polygenic loci.
  //
  bool  loci_ok = false;
  
  if(freq_ok)
  {
       if (no_loci_found){ // block due to JA
          my_max_pgt = nloci*2 + 1; 
          my_means.clear();
          my_means.resize(max_pgt(), QNAN);
          loci_ok = true;  
         } 
    if (! no_loci_found){
    if(nloci != (size_t)(-1))
    {
      if(0 < nloci && nloci <= FPMM_MAX_LOCI)
      {
        my_max_pgt = nloci * 2 + 1;
        loci_ok = true;

        my_means.clear();
        my_means.resize(max_pgt(), QNAN);
      }
      else
      {
        my_errors << priority(critical) << "Parameter loci in fpmm sub-block outside of "
                << "allowable range.  Skipping analysis ..." << endl;
      }
    }
        else
        {
          loci_ok = true;
         } 
   }
  }
                     
  
  if(freq_ok && loci_ok)
  {
    calculate_means();
  }
  
return freq_ok && loci_ok;
}

int
finite_polygenic_mixed_model_sub_model::finalizeConfiguration()
{
  my_parameters.resize(1);
  my_parameters[0] = MAXFUN::ParameterInput
    ( "FPMM",
      "polyg variance",
      my_variance_fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_variance,
      FPMM_LB + FPMM_EPSILON,
      FPMM_UB);

  return 0;
}

//
// ===============  Ancillary functions
//
// - If a fixed input value is not w/i bounds (or QNAN), return false.
//   If a non-fixed input value is not w/i bounds (or QNAN), set value to
//   QNAN.  In either case, write an appropriate message.
//
bool
finite_polygenic_mixed_model_sub_model::input_meets_constraints(model_input& input)
{
  bool  return_value = false;
  if(SAGE::isnan(input.value))
  {
    return_value = true;
  }
  else
  {
    if(FPMM_LB + FPMM_EPSILON <= input.value)
    {
      return_value = true;
    }
    else
    {
      if(input.fixed == true)
      {
        my_errors << priority(critical) << "Parameter var specified "
                  << "in fpmm sub-block is not within allowable limits.  "
                  << "Skipping analysis ..." << endl;
      }
      else
      {
        input.value = QNAN;
        my_errors << priority(error) << "Parameter var specified "
                  << "in fpmm sub-block is not within allowable limits.  "
                  << "Ignoring ..." << endl;
        return_value = true;
      }
    }
  }
  return return_value;  
}


//
// ==============  Synchronization w. Maxfun
//
int
finite_polygenic_mixed_model_sub_model::update()
{
  // If we're not a fixed submodel, get the variance from maxfun and
  // calculate our means.

  if(!my_variance_fixed)
  {
    my_variance = getParam(0);

    calculate_means();
  }
  
  return 0;
}

void  
finite_polygenic_mixed_model_sub_model::calculate_means()
{
  if(finite(my_variance))
  {
    double p = frequency();
    double l = loci();
    double s = variance();

    for(size_t v = 0; v < my_means.size(); ++v)
    {
      double temp1 = 2.0 * p * l;

      double first  = (v - temp1) / (1.0 - p);
      double second = sqrt(s * (1.0 - p) / temp1);

      my_means[v] = first * second;
    }
  }
}

}
}
