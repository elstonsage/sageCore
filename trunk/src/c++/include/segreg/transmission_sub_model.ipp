//============================================================================
// File:      transmission_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/1/01 - created.                             djb
//                                                          
// Notes:     inlines for transmission_sub_model.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_TRANS_SUB_MODEL_H
#include "segreg/transmission_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{


//============================================================================
// IMPLEMENTATION:  transmission_sub_model
//============================================================================
//
inline  
transmission_sub_model::transmission_sub_model
      (const genotype_frequency_sub_model* gf_ptr, cerrorstream& errors)
    : SegregSubmodel(errors), my_gf_ptr(gf_ptr)
{
  //lint -e{534}
  set(homog_no_trans,
      model_input(QNAN, false),
      model_input(QNAN, false),
      model_input(QNAN, false), false, false);
      
  my_default = true;
}


inline
transmission_sub_model::transmission_sub_model
      (const transmission_sub_model& other)
    : SegregSubmodel(other)
{
  my_option = other.my_option;
  //
  my_gf_ptr = other.my_gf_ptr;            //lint !e1554 CAUTION: sub_models will refer to same model.
  my_default = other.my_default;

  for(int i = 0; i < NUM_OF_TYPES; ++i)
  {
    my_taus[i] = other.my_taus[i];
  }
}

inline transmission_sub_model&
transmission_sub_model::operator=
        (const transmission_sub_model& other)
{
  if(this != &other)
  {
    SegregSubmodel::operator=(other);
    
    my_option  = other.my_option;
    my_gf_ptr  = other.my_gf_ptr;           //lint !e1555 CAUTION: sub_models will refer to same model.
    my_default = other.my_default;
    
    for(int i = 0; i < NUM_OF_TYPES; ++i)
    {
      my_taus[i] = other.my_taus[i];
    }
  }
  
  return *this;
}

inline
transmission_sub_model::~transmission_sub_model()
{
  //lint -e{1540}
}

inline transmission_sub_model::sm_option
transmission_sub_model::option() const
{
  return my_option;
}

inline string  
transmission_sub_model::option_description() const
{
  return option_2_description(my_option);
}

inline string  
transmission_sub_model::name() const
{
  return TRANSMISSION_NAME;
}

/// Returns \c true if the transmission probabilities are affected by sex,
/// \c false otherwise.
inline bool
transmission_sub_model::has_sex_effect() const
{
  return option() == mitochondrial;
}

inline string
transmission_sub_model::option_2_description(sm_option opt)
{
  switch(opt)
  {
    case no_trans:
      return "no transmission";
    case homog_no_trans:
      return "homogeneous no transmission";
    case homog_mendelian:
      return "homogeneous mendelian";
    case homog_general:
      return "homogeneous general";
    case general:
      return "general";
    case tau_ab_free:
      return "tau AB free";
    case mitochondrial:
      return "mitochondrial";
    default:
      return "";
  }
}

inline string
transmission_sub_model::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case no_trans:
      return "no_trans";
    case homog_no_trans:
      return "homog_no_trans";
    case homog_mendelian:
      return "homog_mendelian";
    case homog_general:
      return "homog_general";
    case general:
      return "general";
    case tau_ab_free:
      return "tau_ab_free";
    case mitochondrial:
      return "mitochondrial";
    default:
      return "";
  }
}

inline double
transmission_sub_model::tau(genotype_index gt) const
{
  return my_taus[gt];
}

// - Equation 79 Yi Dong's Programmer Documentation (8/3/01).
//   **** Special case 'homogeneity between generations'??
//
inline double          
transmission_sub_model::prob(genotype_index indiv, genotype_index mother,
                             genotype_index father) const
{
  double  return_value = 0.0;

  switch(my_option)
  {
    case no_trans :
    
      assert(my_gf_ptr != 0);
    
      // - This is true for both the hwe and nhwe cases.
      //
      return_value = my_gf_ptr->prob(indiv);
      break;
    case mitochondrial :
      return_value = (indiv == mother);
      break;
    default:
      switch(indiv)
      {
        case index_AA:
          return_value = tau(mother) * tau(father);
          break;
        case index_AB:
          return_value = tau(mother) * (1.0 - tau(father)) + tau(father) * (1.0 - tau(mother));
          break;
        case index_BB:
          return_value = (1.0 - tau(father)) * (1.0 - tau(mother));
          break;
        case index_INVALID:
        default:
          SAGE_internal_error();
      }
  }
  
  return return_value;
}

inline bool
transmission_sub_model::default_() const
{
  return my_default;
}

// - Write sub-model values in LSF readable format.
//
inline void
transmission_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  out << std::boolalpha
      << "# " << name() << "\n"
      << "transmission\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n";

  switch(my_option)
  {
    case homog_no_trans:     // Nothing else to specify.
    case homog_mendelian:
    case no_trans:
      break;
      
    case homog_general:
      assert(! SAGE::isnan(my_taus[index_AA]));
      assert(! SAGE::isnan(my_taus[index_BB]));
    
      out << "  tau=AA, val=" << my_taus[index_AA] << ", fixed=" 
                << (my_parameters[index_AA].initial_type == MAXFUN::Parameter::FIXED) << "\n"
          << "  tau=BB, val=" << my_taus[index_BB] << ", fixed="
                << (my_parameters[index_BB].initial_type == MAXFUN::Parameter::FIXED) << "\n";
      break;
      
    case general:
      assert(! SAGE::isnan(my_taus[index_AA]));
      assert(! SAGE::isnan(my_taus[index_AB]));
      assert(! SAGE::isnan(my_taus[index_BB]));
    
      out << "  tau=AA, val=" << my_taus[index_AA] << ", fixed=" 
                << (my_parameters[index_AA].initial_type == MAXFUN::Parameter::FIXED) << "\n"
          << "  tau=AB, val=" << my_taus[index_AB] << ", fixed="
                << (my_parameters[index_AB].initial_type == MAXFUN::Parameter::FIXED) << "\n"
          << "  tau=BB, val=" << my_taus[index_BB] << ", fixed=" 
                << (my_parameters[index_BB].initial_type == MAXFUN::Parameter::FIXED) << "\n";
      break;
      
    case tau_ab_free:
      assert(! SAGE::isnan(my_taus[index_AB]));
    
      out << "  tau=AB, val=" << my_taus[index_AB] << ", fixed=" 
                << (my_parameters[index_AB].initial_type == MAXFUN::Parameter::FIXED) << "\n";
      break;
    case mitochondrial:
      break;
      
    default:
      SAGE_internal_error();
  }

  out << "}" << std::noboolalpha << std::endl;
  out.precision(old_precision);
}

inline bool transmission_sub_model::is_complete() const
{
  bool complete = false;

  switch(option())
  {
    case homog_no_trans  :
    case no_trans        :
    case mitochondrial   :
    case homog_mendelian : complete = true; break;

    case homog_general   : complete = finite(tau(index_AA)) &&
                                      finite(tau(index_BB));
                           break;

    case tau_ab_free     : complete = finite(tau(index_AB)); //lint !e734
                           break;

    case general         : complete = finite(tau(index_AA)) && 
                                      finite(tau(index_AB)) &&
                                      finite(tau(index_BB));
                           break;
  }

  return complete;
}

inline void  
transmission_sub_model::copy_freq_a_into_taus()
{
  assert(my_gf_ptr != 0);
  
  double  qA = my_gf_ptr->freq_A();
  my_taus[index_AA] = qA;
  my_taus[index_AB] = qA;
  my_taus[index_BB] = qA;
}

}}
