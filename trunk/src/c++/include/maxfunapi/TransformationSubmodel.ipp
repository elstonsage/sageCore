//============================================================================
// File:      TransformationSubmodel.ipp
//                                                                          
// Author:    Dan Baechle
//            Stephen Gross
//                                                                          
// History:   7/2/01 - created.                                djb
//            5/25/04 - Adapted to use the new Maxfun API      sag
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef TRANSFORMATIONSUBMODEL_H
#include "maxfunapi/TransformationSubmodel.h"
#endif

namespace SAGE   {
namespace MAXFUN {

//============================================================================
// Constructor
//============================================================================
inline  
TransformationSubmodel::TransformationSubmodel(cerrorstream & errors)
  : Submodel(errors)
{
  set(box_cox, model_input(QNAN), model_input(QNAN), QNAN, QNAN);
}

//============================================================================
// Copy constructor
//============================================================================
inline
TransformationSubmodel::TransformationSubmodel(const TransformationSubmodel & other)
  : Submodel(other)
{
  copy(other);
}

//============================================================================
// operator=
//============================================================================
inline TransformationSubmodel &
TransformationSubmodel::operator=(const TransformationSubmodel& other)
{
  if(this != &other)
    copy(other);

  return *this;
}

//============================================================================
// copy
//============================================================================
inline void
TransformationSubmodel::copy(const TransformationSubmodel& other)
{
  Submodel::operator=(other);

  my_option               = other.my_option;
  my_geometric_mean       = other.my_geometric_mean;
  clear_each_sync         = other.clear_each_sync;
  
  my_lambda_one = other.my_lambda_one;
  my_lambda_two = other.my_lambda_two;
}

//============================================================================
// Destructor
//============================================================================
inline
TransformationSubmodel::~TransformationSubmodel()
{}

//============================================================================
// option
//============================================================================
inline TransformationSubmodel::sm_option
TransformationSubmodel::option() const
{
  return my_option;
}

//============================================================================
// option_description
//============================================================================
inline string  
TransformationSubmodel::option_description() const
{
  return option_2_description(my_option);
}

//============================================================================
// name
//============================================================================
inline string  
TransformationSubmodel::name() const
{
  return TRANSFORMATION_NAME;
}

//============================================================================
// option_2_description
//============================================================================
inline string
TransformationSubmodel::option_2_description(sm_option opt)
{
  switch(opt)
  {
    case box_cox:       return "Box-Cox";
    case george_elston: return "George-Elston";
    case no_trans:      return "none";
    default:            return "";
  }
}

//============================================================================
// option_2_parameter
//============================================================================
inline string
TransformationSubmodel::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case box_cox:       return "box_cox";
    case george_elston: return "george_elston";
    case no_trans:      return "none";
    default:            return "";
  }
}

//============================================================================
// lambda_one
//============================================================================
inline double
TransformationSubmodel::lambda_one() const
{
  return my_lambda_one;
}

//============================================================================
// lambda_two
//============================================================================
inline double
TransformationSubmodel::lambda_two() const
{
  return my_lambda_two;
}

//============================================================================
// dump
//============================================================================
inline void
TransformationSubmodel::dump(ostream & out) const
{
  int old_precision = out.precision();

  // *** XXX ***
  // The 12 (below) was originally DUMP_PRECISION in segreg, but had to be
  // removed when moving tranformation out of segreg.  It should probably be
  // a functional arguement. -- GCW 2003-01-27.

  out.precision(12);

  bool lambda_one_fixed = my_parameters[0].initial_type == Parameter::FIXED;
  bool lambda_two_fixed = my_parameters[1].initial_type == Parameter::FIXED;

  assert(!SAGE::isnan(lambda_one()));

  out << "# " << name() << "\n"
      << "transformation\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n"
      << "  lambda1, val=" << lambda_one() << ", fixed=" << std::boolalpha << lambda_one_fixed 
      << ", lower_bound=" << my_parameters[0].lower_bound;
      
  if(finite(my_parameters[0].upper_bound))
    out << ", upper_bound=" << my_parameters[0].upper_bound;
      
  assert(!SAGE::isnan(lambda_two()));
      
  out << "\n"
      << "  lambda2, val=" << lambda_two() << ", fixed=" << lambda_two_fixed << "\n"
      << "}" << std::noboolalpha << endl;

  out.precision(old_precision);
}
  
//=============================================================================
// sign
//=============================================================================
inline double
TransformationSubmodel::sign(double value)
{
  return value < 0.0 ? -1.0 : (value == 0.0 ? 0.0 : 1.0);
}

//=============================================================================
// calculate_geom_mean
//=============================================================================
inline bool
TransformationSubmodel::calculate_geom_mean(const std::vector<double>& traits) const
{
  switch(my_option)
  {
    case no_trans:      return true;
    case box_cox:       my_geometric_mean = bc_geom_mean(traits); break;
    case george_elston: my_geometric_mean = ge_geom_mean(traits); break;
    default:            SAGE_internal_error();
  }

  return !SAGE::isnan(my_geometric_mean);
}

/// Returns true if lambda1 is close to 0 (less than 1e-9 away).  In this case,
/// we use the 0 formulation of the transform.
inline bool TransformationSubmodel::lambda_one_close_to_zero() const
{
  return (abs(lambda_one()) < 1e-9);
}

//=============================================================================
// bc_transform_power_zero
//=============================================================================
inline void  
TransformationSubmodel::bc_transform_power_zero(double G, double& trait) const
{
  if(!SAGE::isnan(trait))
    trait = G * log(trait + lambda_two());
}

//=============================================================================
// bc_transform_power_non_zero
//=============================================================================
inline void  
TransformationSubmodel::bc_transform_power_non_zero
    (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double numerator   = pow(trait + lambda_two(), lambda_one()) - 1,
           denominator = lambda_one() * pow(G, lambda_one() - 1);

    trait = numerator / denominator;
  }
}

//=============================================================================
// ge_transform_power_zero
//=============================================================================
inline void  
TransformationSubmodel::ge_transform_power_zero
  (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double adj_trait = trait + lambda_two(),
           sn        = sign(adj_trait),
           abs_trait = (adj_trait < 0.0) ? -adj_trait : adj_trait;

    trait = sn * G * log(abs_trait + 1);
  }
}

//=============================================================================
// ge_transform_power_non_zero
//=============================================================================
inline void  
TransformationSubmodel::ge_transform_power_non_zero
    (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double adj_trait   = trait + lambda_two(),
           sn          = sign(adj_trait),
           abs_trait   = (adj_trait < 0.0) ? -adj_trait : adj_trait,
           numerator   = sn * (pow(abs_trait + 1, lambda_one()) - 1),
           denominator = lambda_one() * pow(G, lambda_one() - 1);

    trait = numerator / denominator;
  }
}

} // End namespace MAXFUN
} // End namespace SAGE

