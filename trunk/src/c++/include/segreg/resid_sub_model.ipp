//============================================================================
// File:      resid_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   8/2/01 - created.                                djb
//                                                                          
// Notes:     inlines for residual_correlation_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_RESID_SUB_MODEL_H
#include "segreg/resid_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{

//============================================================================
// IMPLEMENTATION:  residual_correlation_sub_model
//============================================================================
//
inline  
residual_correlation_sub_model::residual_correlation_sub_model
      (cerrorstream& errors)
    : SegregSubmodel(errors),
      my_model_class(model_A),
      my_is_in_default_mode(true)
{
  //lint -e{534}
  set_equal_po_ss(model_input(), model_input(), model_input(), model_input());  

  initialize();
  
  my_is_in_default_mode = true;
}

inline
residual_correlation_sub_model::residual_correlation_sub_model
      (const residual_correlation_sub_model& other)
    : SegregSubmodel(other)
{
  my_option = other.my_option;
  my_is_in_default_mode = other.my_is_in_default_mode;
  my_model_class = other.my_model_class;
  
  for(int i = 0; i < NUM_OF_CORRS; ++i)
  {
    my_corrs[i] = other.my_corrs[i];
    my_corrs_fixed[i] = other.my_corrs_fixed[i];
  }

  for(int i = 0; i < 4; ++i)
  {
    my_alpha_mother[i] = other.my_alpha_mother[i];
    my_alpha_father[i] = other.my_alpha_father[i];
    my_delta       [i] = other.my_delta       [i];
  }
}

inline residual_correlation_sub_model&
residual_correlation_sub_model::operator=
        (const residual_correlation_sub_model& other)
{
  if(this != &other)
  {
    SegregSubmodel::operator=(other);

    my_option = other.my_option;
    my_is_in_default_mode = other.my_is_in_default_mode;
    my_model_class = other.my_model_class;
    
    for(int i = 0; i < NUM_OF_CORRS; ++i)
    {
      my_corrs[i] = other.my_corrs[i];
      my_corrs_fixed[i] = other.my_corrs_fixed[i];
    }

    for(int i = 0; i < 4; ++i)
    {
      my_alpha_mother[i] = other.my_alpha_mother[i];
      my_alpha_father[i] = other.my_alpha_father[i];
      my_delta       [i] = other.my_delta       [i];
    }
  }
  
  return *this;
}

inline
residual_correlation_sub_model::~residual_correlation_sub_model()
{}

inline residual_correlation_sub_model::sm_option
residual_correlation_sub_model::option() const
{
  return my_option;
}

inline string  
residual_correlation_sub_model::option_description() const
{
  return option_2_description(my_option);
}

inline string  
residual_correlation_sub_model::name() const
{
  return RESID_NAME;
}

inline bool
residual_correlation_sub_model::is_in_default_mode() const
{
  return my_is_in_default_mode;
}

inline string
residual_correlation_sub_model::option_2_description(sm_option opt)
{
  switch(opt)
  {
    case equal_po_ss:
      return "No spouse correlation, parent-offspring and sib-sib correlations equal";
    case equal_po:
      return "Spouse and sib-sib correlations independent, mother-"
             "offspring and father-offspring correlations equal";
    case arb:
      return "Arbitrary";
    default:
      return "";
  }
}

inline string
residual_correlation_sub_model::option_2_parameter(sm_option opt)
{
  switch(opt)
  {
    case equal_po_ss:
      return "equal_po_ss";
    case equal_po:
      return "equal_po";
    case arb:
      return "arb";
    default:
      return "";
  }
}

inline double          
residual_correlation_sub_model::father_mother_correlation() const
{
  return my_corrs[fm];
}
  
inline double          
residual_correlation_sub_model::mother_son_correlation() const
{
  return my_corrs[ms];
}
    
inline double          
residual_correlation_sub_model::father_son_correlation() const
{
  return my_corrs[fs];
}

inline double          
residual_correlation_sub_model::mother_daughter_correlation() const
{
  return my_corrs[md];
}
 
inline double         
residual_correlation_sub_model::father_daughter_correlation() const
{
  return my_corrs[fd];
}
    
inline double          
residual_correlation_sub_model::brother_brother_correlation() const
{
  return my_corrs[bb];
}
    
inline double          
residual_correlation_sub_model::sister_sister_correlation() const
{
  return my_corrs[ss];
}
 
inline double          
residual_correlation_sub_model::brother_sister_correlation() const
{
  return my_corrs[bs];
}

inline double          
residual_correlation_sub_model::correlation(corr rs) const
{
  return my_corrs[rs];
}
    
inline bool
residual_correlation_sub_model::is_correlation_fixed(corr rs) const
{
  return my_corrs_fixed[rs];
}

inline bool 
residual_correlation_sub_model::has_residuals() const
{
  if(my_model_class == model_FPMM) return false;
  
  for(int i = 0; i != NUM_OF_CORRS; ++i)
  {
    // If the correlations aren't fixed, or a correlation is fixed, but not
    // to 0.0, we have residuals.
    
    if(!my_corrs_fixed[i] || my_corrs[i] != 0.0)
      return true;
  }

  return false;
}

/// Returns \c true if there is a sex effect, false otherwise
inline
bool
  residual_correlation_sub_model::has_sex_effect() const
{
  return my_model_class != model_FPMM && option() == arb;
}

inline double          
residual_correlation_sub_model::alpha_mother(double t_f, double t_m) const
{
  return alpha_mother(!SAGE::isnan(t_f), !SAGE::isnan(t_m));
}

inline double          
residual_correlation_sub_model::alpha_father(double t_f, double t_m) const
{
  return alpha_father(!SAGE::isnan(t_f), !SAGE::isnan(t_m));
}

inline double          
residual_correlation_sub_model::delta(double t_f, double t_m) const
{
  return delta(!SAGE::isnan(t_f), !SAGE::isnan(t_m));
}

inline double          
residual_correlation_sub_model::alpha_mother(bool t_f, bool t_m) const
{
  return my_alpha_mother[2 * t_f + t_m];
}

inline double          
residual_correlation_sub_model::alpha_father(bool t_f, bool t_m) const
{
  return my_alpha_father[2 * t_f + t_m];
}

inline double          
residual_correlation_sub_model::delta(bool t_f, bool t_m) const
{
  return my_delta[2 * t_f + t_m];
}

// - Write sub-model values in LSF readable format.
//
inline void
residual_correlation_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  bool  corr_fm_fixed = false;
  bool  corr_mo_fixed = false;
  bool  corr_fo_fixed = false;
  bool  corr_ss_fixed = false;

  switch(my_option)
  {
    case equal_po_ss:
      corr_fm_fixed  = my_corrs_fixed[fm];
      corr_mo_fixed  = my_corrs_fixed[ms];
      corr_fo_fixed  = my_corrs_fixed[fs];
      corr_ss_fixed  = my_corrs_fixed[bb];
      break;
      
    case equal_po:
      corr_fm_fixed  = my_corrs_fixed[fm];
      corr_mo_fixed  = my_corrs_fixed[ms];
      corr_fo_fixed  = my_corrs_fixed[fs];
      corr_ss_fixed  = my_corrs_fixed[bb];
      break;
      
    case arb:
      corr_fm_fixed  = my_corrs_fixed[fm];
      corr_mo_fixed  = my_corrs_fixed[fs];
      corr_fo_fixed  = my_corrs_fixed[ms];
      corr_ss_fixed  = my_corrs_fixed[bb];
      break;
      
    default:
      SAGE_internal_error();
  }

  assert(! SAGE::isnan(my_corrs[fm]));
  assert(! SAGE::isnan(my_corrs[ms]));
  assert(! SAGE::isnan(my_corrs[fs]));
  assert(! SAGE::isnan(my_corrs[bb]));

  out << "# " << name() << "\n"
      << "resid\n" 
      << "{\n"
      << "  # " << option_2_description(my_option) << "\n"
      << "  option=" << option_2_parameter(my_option) << "\n"
      << "  fm, val=" << my_corrs[fm] << ", fixed=" << std::boolalpha << corr_fm_fixed << "\n" 
      << "  mo, val=" << my_corrs[ms] << ", fixed=" << corr_mo_fixed << "\n" 
      << "  fo, val=" << my_corrs[fs] << ", fixed=" << corr_fo_fixed << "\n" 
      << "  ss, val=" << my_corrs[bb] << ", fixed=" << corr_ss_fixed << "\n" 
      << "}" << std::noboolalpha << std::endl;
      
  out.precision(old_precision);
}

inline void residual_correlation_sub_model::calculate_alpha_and_delta()
{
  // Syntactic simplicity declarations

  double Pms = my_corrs[ms];
  double Pfs = my_corrs[fs];
  double Pfm = my_corrs[fm];

  // Make sure the values won't break

  double denom = (1 - Pfm * Pfm);

  if(denom < 1.e-6) denom = 1.e-6;          // added 9-25-1 djb
                                            // modified 2003-05-05 gcw

  // For each of alpha_m, alpha_f, and delta, there are 4 options, depending
  // on which of the parents are available.  These are organized bitwise,
  // such that:
  //
  // Case 0: neither parent available
  // Case 1: mother available, father isn't
  // Case 2: father available, mother isn't
  // Case 3: both parents available

  // Case 0:

  my_alpha_mother[0] = 0.0;
  my_alpha_father[0] = 0.0;
  my_delta       [0] = 0.0;

  // Case 1:

  my_alpha_mother[1] = Pms;
  my_alpha_father[1] = 0.0;
  my_delta       [1] = Pms * Pms;
  
  // Case 2:

  my_alpha_mother[2] = 0.0;
  my_alpha_father[2] = Pfs;
  my_delta       [2] = Pfs * Pfs;
  
  // Case 3:

  my_alpha_mother[3] = (Pms - Pfm * Pfs) / denom;
  my_alpha_father[3] = (Pfs - Pfm * Pms) / denom;
  my_delta       [3] = (Pms * Pms + Pfs * Pfs - 2.0 * Pfm * Pms * Pfs) / denom;
}

inline double
residual_correlation_sub_model::get_lower_bound() const
{
  if(my_model_class == model_MLM) return NEGATIVE_INF;
  else                            return RESID_LB;
}

inline double
residual_correlation_sub_model::get_upper_bound() const
{
  if(my_model_class == model_MLM) return POSITIVE_INF;
  else                            return RESID_UB;
}

}}
