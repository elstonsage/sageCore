//============================================================================
// File:      ascertainment_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/18/01  - created.                        djb
//                                                                          
// Notes:     inlines for ascertainment_sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_ASCER_SUB_MODEL_H
#include "segreg/ascertainment_sub_model.h"
#endif

namespace SAGE
{

namespace SEGREG
{


//============================================================================
// IMPLEMENTATION:  ascertainment_sub_model
//============================================================================
//
inline  
ascertainment_sub_model::ascertainment_sub_model
      (cerrorstream& errors)
    : my_errors(errors)
{
  //lint -e{534}
  set(none, not_specified, 0, string(), std::vector<double>(), string(), 
      QNAN, QNAN, QNAN, QNAN);
}

inline
ascertainment_sub_model::ascertainment_sub_model
      (const ascertainment_sub_model& other)
{
  my_errors   = other.my_errors;
  
  my_s_option = other.my_s_option;
  my_v_option = other.my_v_option;
  my_psf_indicator = other.my_psf_indicator;
  my_thresh_indicator = other.my_thresh_indicator;
  my_includes = other.my_includes;
  my_thresh = other.my_thresh;               
  my_thresh_high = other.my_thresh_high;             
  my_thresh_low = other.my_thresh_low;
  my_indic_thresh = other.my_indic_thresh;            
  my_t_specified = other.my_t_specified;       
  my_th_specified = other.my_th_specified;        
  my_tl_specified = other.my_tl_specified;
}

inline ascertainment_sub_model&
ascertainment_sub_model::operator=
        (const ascertainment_sub_model& other)
{
  if(this != &other)
  {
    my_errors = other.my_errors;
  
    my_s_option = other.my_s_option;
    my_v_option = other.my_v_option;
    my_psf_indicator = other.my_psf_indicator;
    my_thresh_indicator = other.my_thresh_indicator;
    my_includes = other.my_includes;
    my_thresh = other.my_thresh;               
    my_thresh_high = other.my_thresh_high;             
    my_thresh_low = other.my_thresh_low;            
    my_t_specified = other.my_t_specified;       
    my_th_specified = other.my_th_specified;        
    my_tl_specified = other.my_tl_specified;    
    my_indic_thresh = other.my_indic_thresh;
  }
  
  return *this;
}

inline
ascertainment_sub_model::~ascertainment_sub_model()
{}
    
inline ascertainment_sub_model::s_sm_option     
ascertainment_sub_model::s_option() const
{
  return my_s_option;
}

inline ascertainment_sub_model::v_sm_option     
ascertainment_sub_model::v_option() const
{
  return my_v_option;
}

inline string  
ascertainment_sub_model::option_description() const
{
  return s_option_description() + ", " + v_option_description();
}

inline string  
ascertainment_sub_model::s_option_description() const
{
  return s_option_2_description(my_s_option);
}

inline string  
ascertainment_sub_model::v_option_description() const
{
  return v_option_2_description(my_v_option);
}


inline string  
ascertainment_sub_model::s_option_2_description(s_sm_option option)
{
  switch(option)
  {
    case none:
      return "none";
    case founders:
      return "founders";
    case psf:
      return "proband sampling frame";
    case founders_plus_psf:
      return "founders plus proband sampling frame";
    default:
      return "";
  }
}

inline string  
ascertainment_sub_model::s_option_2_parameter(s_sm_option option)
{
  switch(option)
  {
    case none:
      return "none";
    case founders:
      return "founders";
    case psf:
      return "psf";
    case founders_plus_psf:
      return "founders_plus_psf";
    default:
      return "";
  }
}

inline string  
ascertainment_sub_model::v_option_2_description(v_sm_option option)
{
  switch(option)
  {
    case actual:
      return "actual";
    case gte_thresh:
      return "greater than or equal to threshold";
    case lte_thresh:
      return "less than or equal to threshold";
    case thresh_indic:
      return "threshold indicator";
    case onset:
      return "by onset";
    case not_specified:
    default:
      return "";
  }
}

inline string  
ascertainment_sub_model::v_option_2_parameter(v_sm_option option)
{
  switch(option)
  {
    case actual:
      return "actual";
    case gte_thresh:
      return "gte_thresh";
    case lte_thresh:
      return "lte_thresh";
    case thresh_indic:
      return "thresh_indic";
    case onset:
      return "by_onset";
    case not_specified:
    default:
      return "";
  }
}
    
inline string  
ascertainment_sub_model::name() const
{
  return ASCERTAINMENT_NAME;
}

inline string 
ascertainment_sub_model::psf_indicator() const
{
  return my_psf_indicator;
}

inline const std::vector<double>&
ascertainment_sub_model::includes() const
{
  return my_includes;
}

inline string          
ascertainment_sub_model::thresh_indicator() const
{
  return my_thresh_indicator;
}

inline double          
ascertainment_sub_model::thresh() const
{
  return my_thresh;
}
 
inline double          
ascertainment_sub_model::thresh_high() const
{
  return my_thresh_high;
}

inline double          
ascertainment_sub_model::thresh_low() const
{
  return my_thresh_low;
}

inline double          
ascertainment_sub_model::indic_thresh() const
{
  return my_indic_thresh;
}

    
// - Write sub-model values in LSF readable format.
//
inline void
ascertainment_sub_model::dump(std::ostream& out) const
{
  int  old_precision = out.precision();
  out.precision(DUMP_PRECISION);

  out << "# " << name() << "\n"
      << "ascertainment\n" 
      << "{\n"
      << "  # " << s_option_2_description(my_s_option) << "\n"
      << "  cond_set=" << s_option_2_parameter(my_s_option) << "\n";
      
  if(my_psf_indicator.size())
  {
    out << "  psf_indic=" << my_psf_indicator << "\n";
    
    std::vector<double>::const_iterator  iter;
    for(iter = my_includes.begin(); iter != my_includes.end(); iter++)
    {
      out << "  psf_indic_include=" << *iter << "\n";
    }
  }
  
  if(my_v_option != ascertainment_sub_model::not_specified)
  {
    out << endl 
        << "  # " << v_option_2_description(my_v_option) << endl;
    out << "  cond_val=";
    
    string vopt = v_option_2_parameter(my_v_option); 
    out << vopt;
    if(my_t_specified)
    {
      out << ", thresh=" << my_thresh;
    }
    if(my_th_specified)
    {
      out << ", thresh_indic_high=" << my_thresh_high;
    }
    if(my_tl_specified)
    {
      out << ", thresh_indic_low=" << my_thresh_low;
    }
    
    out << endl;        
  }
  
  if(my_thresh_indicator.size())
  {
    out << "  thresh_indic=" << my_thresh_indicator;
    out << ", thresh=" << my_indic_thresh << "\n";
  }
      
  out << "}" << std::endl;
  out.precision(old_precision);
}

inline bool ascertainment_sub_model::is_ind_in_C
    (FPED::MemberConstPointer m) const
{
  // If the thresh_indic isn't valid, we can't use the person (missing)
  // NOTE that if the v_option is not thresh_indic, this always returns
  // true.
  if(!does_ind_have_valid_thresh(m)) return false;

  // Check to see if the individual is in the conditioned subset
  bool use = false;

  // NOTE: For ascertainment individuals are only founders if they actually
  // belong to a subpedigree.  Unconnecteds are not founders

  switch(s_option())
  {
    case psf :
      use = is_ind_in_PSF(m);
      break;

    case founders:

      use = (m->parent1() == NULL && m->subindex() != (uint) -1);
      break;

    case founders_plus_psf:
      use = (m->parent1() == NULL && m->subindex() != (uint) -1) || is_ind_in_PSF(m);
      break;

    case none:
    default:
      use = false;
  }

  return use;
}

inline ascertainment_sub_model::v_sm_option
    ascertainment_sub_model::get_ind_type
    (FPED::MemberConstPointer m) const
{
  // If the individual is not in C, then it's always missing

  if(!is_ind_in_C(m))
  {
    return not_specified;
  }

  // If the individual is in C through founderness

  if(s_option() == founders || (s_option() == founders_plus_psf && !is_ind_in_PSF(m)))
  {
    return actual;
  }

  // Now we know that the individual is in PSF

  switch(v_option())
  {
    case gte_thresh   : 
      return gte_thresh;

    case lte_thresh   :
      return lte_thresh;

    case actual       : 
      return actual;

    case thresh_indic :
    {
      double thr = get_ind_thresh(m);

      if(thr == indic_thresh())
        return actual;

      return (thr < indic_thresh()) ? lte_thresh : gte_thresh;
    }
    case onset :
      return onset;

    case not_specified :
    default            : return not_specified;
  }
}

inline bool ascertainment_sub_model::is_ind_in_PSF
    (FPED::MemberConstPointer m) const
{
  // Get the member's trait value

  double psf_value = get_ind_trait_value(m, psf_indicator());

  // If the value is missing, then the ind can't be in PSF

  if(SAGE::isnan(psf_value))
    return false;

  // Individual is in psf if psf_value is any of the included elements

  for(size_t j = 0; j < my_includes.size(); ++j)
    if(psf_value == my_includes[j])
      return true; 

  return false;
}

inline bool ascertainment_sub_model::does_ind_have_valid_thresh
    (FPED::MemberConstPointer m) const
{
  // If the v_option isn't thresh_indic, there's nothing to do here
  if(v_option() != thresh_indic)
    return true;

  // Return true if the member has a valid thresh trait (not missing)

  return (!SAGE::isnan(get_ind_thresh(m)));
}

inline double ascertainment_sub_model::get_ind_thresh
    (FPED::MemberConstPointer m) const
{
  // If the v_option isn't thresh_indic, there's nothing to do here
  if(v_option() != thresh_indic)
    return indic_thresh();

  // Get member and pedigree information for later lookups

  return get_ind_trait_value(m, thresh_indicator());
}

inline double ascertainment_sub_model::get_ind_trait_value
    (FPED::MemberConstPointer m, const string& trait) const
{
  // Get member and pedigree information for later lookups

  size_t member_index = m->index();

  const FPED::PedigreeConstPointer ped  = m->pedigree();

  // Get the member's thresh trait
  
  size_t trait_num = m->multipedigree()->info().trait_find(trait);

  double trait_value = ped->info().trait(member_index, trait_num);

  return trait_value;
}

inline bool ascertainment_sub_model::ind_uses_thresh_high
    (FPED::MemberConstPointer m) const
{
  return get_ind_type(m) == gte_thresh;
}

inline bool ascertainment_sub_model::ind_uses_thresh_low
    (FPED::MemberConstPointer m) const
{
  return get_ind_type(m) == lte_thresh;
}

inline ostream& operator<< (ostream& out, const ascertainment_sub_model& sm)
{
  out << "\n" << sm.name() << " values: \n";
  out << "Option: " << sm.option_description() << std::endl;
  
  return out;
}

}}
