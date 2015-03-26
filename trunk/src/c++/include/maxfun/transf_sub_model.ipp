//============================================================================
// File:      transformation_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   7/2/01 - created.                                djb
//                                                                          
// Notes:     inlines for transformation_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef TRANSF_SUB_MODEL_H
#include "maxfun/transf_sub_model.h"
#endif

namespace SAGE
{

//============================================================================
// IMPLEMENTATION:  transformation_sub_model
//============================================================================
//
inline  
transformation_sub_model::transformation_sub_model
      (cerrorstream& errors)
    : sub_model(errors)
{
  //lint -e{534}
  set(box_cox, model_input(QNAN), model_input(QNAN), QNAN, QNAN);
}

inline
transformation_sub_model::transformation_sub_model
      (const transformation_sub_model& other)
    : sub_model(other)
{
  my_option = other.my_option;
  my_lambda_one = other.my_lambda_one;
  my_lambda_two = other.my_lambda_two;

  my_geometric_mean = QNAN;
}

inline transformation_sub_model&
transformation_sub_model::operator=
        (const transformation_sub_model& other)
{
  if(this != &other)
  {
    sub_model::operator=(other);

    my_option = other.my_option;
    my_lambda_one = other.my_lambda_one;
    my_lambda_two = other.my_lambda_two;

    my_geometric_mean = other.my_geometric_mean;
  }
  
  return *this;
}

inline
transformation_sub_model::~transformation_sub_model()
{}

inline transformation_sub_model::sm_option
transformation_sub_model::option() const
{
  return my_option;
}

inline string  
transformation_sub_model::option_description() const
{
  return option_2_description(my_option);
}

inline string  
transformation_sub_model::name() const
{
  return TRANSFORMATION_NAME;
}

inline string
transformation_sub_model::option_2_description(sm_option opt)
{
  switch(opt)
  {
    case box_cox:
      return "Box-Cox";
    case george_elston:
      return "George-Elston";
    case no_trans:
      return "none";
    default:
      return "";
  }
}

inline string
transformation_sub_model::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case box_cox:
      return "box_cox";
    case george_elston:
      return "george_elston";
    case no_trans:
      return "none";
    default:
      return "";
  }
}

inline double
transformation_sub_model::lambda_one() const
{
  return my_lambda_one;
}

inline double
transformation_sub_model::lambda_two() const
{
  return my_lambda_two;
}

// - Write sub-model values in LSF readable format.
//
inline void
transformation_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();

  // *** XXX ***
  // The 12 (below) was originally DUMP_PRECISION in segreg, but had to be
  // removed when moving tranformation out of segreg.  It should probably be
  // a functional arguement. -- GCW 2003-01-27.

  out.precision(12);

  bool  lambda_one_fixed = my_parameters[0].status == fixed;
  bool  lambda_two_fixed = my_parameters[1].status == fixed;

  assert(! SAGE::isnan(my_lambda_one));

  out << "# " << name() << "\n"
      << "transformation\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n"
      << "  lambda1, val=" << my_lambda_one << ", fixed=" << std::boolalpha << lambda_one_fixed 
      << ", lower_bound=" << my_parameters[0].lower_bound;
      
  if(finite(my_parameters[0].upper_bound))
  {
    out << ", upper_bound=" << my_parameters[0].upper_bound;
  }
      
  assert(! SAGE::isnan(my_lambda_two));
      
  out << "\n"
      << "  lambda2, val=" << my_lambda_two << ", fixed=" << lambda_two_fixed << "\n"
      << "}" << std::noboolalpha << std::endl;

  out.precision(old_precision);
}
  
inline double
transformation_sub_model::sign(double value)
{
  return value < 0.0 ? -1.0 : (value == 0.0 ? 0.0 : 1.0);
}

inline bool
transformation_sub_model::calculate_geom_mean(const std::vector<double>& traits) const
{
  switch(my_option)
  {
    case no_trans:
      return true;
      
    case box_cox:

      my_geometric_mean = bc_geom_mean(traits);

      break;
      
    case george_elston:
      my_geometric_mean = ge_geom_mean(traits);

      break;
      
    default:
      SAGE_internal_error();
  }

  return !SAGE::isnan(my_geometric_mean);
}

inline void  
transformation_sub_model::bc_transform_power_zero(double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    trait = G *log(trait + my_lambda_two);
  }
}

inline void  
transformation_sub_model::bc_transform_power_non_zero
    (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double  num = pow(trait + my_lambda_two, my_lambda_one) - 1;
    double  denum = my_lambda_one * pow(G, my_lambda_one - 1);
    trait = num / denum;
  }
}

inline void  
transformation_sub_model::ge_transform_power_zero
  (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double adj_trait = trait + my_lambda_two;

    double sn        = sign(adj_trait);
  
    double abs_trait = (adj_trait < 0.0) ? -adj_trait : adj_trait;

    trait = sn * G * log(abs_trait + 1);
  }
}

inline void  
transformation_sub_model::ge_transform_power_non_zero
    (double G, double& trait) const
{
  if(!SAGE::isnan(trait))
  {
    double adj_trait = trait + my_lambda_two;
    double sn        = sign(adj_trait);
    double abs_trait = (adj_trait < 0.0) ? -adj_trait : adj_trait;

    double  num = sn * (pow(abs_trait + 1, my_lambda_one) - 1);
    double  denum = my_lambda_one * pow(G, my_lambda_one -1);
    trait = num / denum;
  }
}

}
