#include "segreg/FPMMSubmodel.h"

using namespace std;

namespace SAGE   {
namespace SEGREG {

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
// IMPLEMENTATION:  FPMMSubmodel
//============================================================================
//
// ===============  Setting the sub-model
//
bool
FPMMSubmodel::set(model_input var, double freq, size_t nloci)
{
  if(my_in_use)
  {
    return false;
  }

  if(! input_meets_constraints(var))
  {
    return false;
  }

  // - Defaults.
  //
  my_parameters.resize(1);
                                       
  my_parameters[0] = ParameterInput("FPMM", 
                                    "polygenic variance",
                                    MAXFUN::Parameter::INDEPENDENT, 
                                    FPMM_DEFAULT_VALUE, 
                                    FPMM_LB + FPMM_EPSILON,
                                    FPMM_UB, true);

  my_frequency = .5;
  
  // - Variance.
  //
  if(! SAGE::isnan(var.value))
  {
    my_parameters[0].initial_estimate = var.value;
    my_parameters[0].initial_type     = var.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
  }
  
  // - Frequency.
  //
  bool  freq_ok = false;
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
  
  // - Number of polygenic loci.
  //
  bool  loci_ok = false;
  if(freq_ok)
  {
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
  
  if(freq_ok && loci_ok)
  {
    internally_synchronize();
  }
  
  return freq_ok && loci_ok;
}

//
// ===============  Ancillary functions
//
// - If a fixed input value is not w/i bounds (or QNAN), return false.
//   If a non-fixed input value is not w/i bounds (or QNAN), set value to
//   QNAN.  In either case, write an appropriate message.
//
bool
FPMMSubmodel::input_meets_constraints(model_input& input)
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


int
FPMMSubmodel::update()
{
  my_variance = getParam(0);

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

} // End namespace SEGREG
} // End namespace SAGE
