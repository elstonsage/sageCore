//============================================================================
// File:      onset_sub_model.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/28/01 - created.                             djb
//                                                          
// Notes:     inlines for onset_sub_model.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef SEGREG_ONSET_SUB_MODEL_H
#define SEGREG_ONSET_SUB_MODEL_H
#endif

#include "onset_sub_model.h"

namespace SAGE
{

namespace SEGREG
{


//============================================================================
// IMPLEMENTATION:  onset_sub_model
//============================================================================
//
inline  
onset_sub_model::onset_sub_model(cerrorstream& errors)
    : my_errors(errors), my_t_option(t_A), my_m_option(m_N)
{
  set_a_option();
}

inline
onset_sub_model::onset_sub_model(const onset_sub_model& other)
    : my_errors(other.my_errors)
{
  my_t_option = other.my_t_option;
  my_m_option = other.my_m_option;
  my_affection_status = other.my_affection_status;
  my_age_of_onset = other.my_age_of_onset;
  my_age_at_exam = other.my_age_at_exam;

  set_a_option();
}

inline onset_sub_model&
onset_sub_model::operator=(const onset_sub_model& other)
{
  if(this != &other)
  {
    my_errors = other.my_errors;
    
    my_t_option = other.my_t_option;
    my_m_option = other.my_m_option;
    my_affection_status = other.my_affection_status;
    my_age_of_onset = other.my_age_of_onset;
    my_age_at_exam = other.my_age_at_exam;

    set_a_option();
  }
  
  return *this;
}

inline
onset_sub_model::~onset_sub_model()
{}

// due to JA
inline bool
onset_sub_model::no_poly_loci()
{
 if (m_option_2_parameter(my_m_option) == "N")
     {return true;} else {return false;} 
}

/*
inline bool
onset_sub_model::static_no_poly_loci()
{

 multi_option m_opt;
   switch(m_opt) 
    {
      case m_N:
       return true;
       break;
      case m_A:
       return false;
       break;
      case m_S:
       return false;
       break;
      default:
       return false; 
     }

 }
// back to regular code
*/
inline onset_sub_model::type_option
onset_sub_model::t_option() const
{
  return my_t_option;
}

inline onset_sub_model::multi_option
onset_sub_model::m_option() const
{
  return my_m_option;
}

inline onset_sub_model::ageon_option
onset_sub_model::a_option() const
{
  return my_a_option;
}

inline string
onset_sub_model::option_description() const
{
  return t_option_description() + ", " + m_option_description();
}

inline string  
onset_sub_model::t_option_description() const
{
  return t_option_2_description(my_t_option);
}

inline string  
onset_sub_model::m_option_description() const
{
  return m_option_2_description(my_m_option);
}

inline string  
onset_sub_model::name() const
{
  return ONSET_NAME;
}

inline string
onset_sub_model::t_option_2_description(type_option opt)
{
  string  return_string;
  switch(opt)
  {
    case t_A:
      return_string = "age of onset depends on type";
      break;
    case t_S:
      return_string = "susceptibility depends on type";
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_string;
}

inline string
onset_sub_model::m_option_2_description(multi_option opt)
{
  string  return_string;
  switch(opt)
  {
    case m_N:
      return_string = "no polygenic component";
      break;
    case m_A:
      return_string = "age of onset has a polygenic component";
      break;
    case m_S:
      return_string = "susceptibility has a polygenic component";
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_string;
}

inline string
onset_sub_model::t_option_2_parameter(type_option opt)
{
  string  return_string;
  switch(opt)
  {
    case t_A:
      return_string = "A";
      break;
    case t_S:
      return_string = "S";
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_string;
}

inline string
onset_sub_model::m_option_2_parameter(multi_option opt)
{
  string  return_string;
  switch(opt)
  {
    case m_N:
      return_string = "N";
      break;
    case m_A:
      return_string = "A";
      break;
    case m_S:
      return_string = "s";
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_string;
}

inline string
onset_sub_model::affection_status() const
{
  return my_affection_status;
}

inline string
onset_sub_model::age_of_onset() const
{
  return my_age_of_onset;
}

inline string
onset_sub_model::age_at_exam() const
{
  return my_age_at_exam;
}

// - Write sub-model values in LSF readable format.
//
inline void
onset_sub_model::dump(std::ostream& out) const
{
  out << "# " << name() << "\n"
      << "onset\n" 
      << "{\n"
      << "  # " << t_option_2_description(my_t_option) << "\n"
      << "  type_dependent=" << t_option_2_parameter(my_t_option) << "\n\n"
      << "  # " << m_option_2_description(my_m_option) << "\n"
      << "  multi_dependent=" << m_option_2_parameter(my_m_option) << "\n\n";
      
  my_affection_status.size() ? out << "  status=" << my_affection_status << "\n" : out << "";
  my_age_of_onset.size() ? out << "  age_onset=" << my_age_of_onset << "\n" : out << "";
  my_age_at_exam.size() ? out << "  age_exam=" << my_age_at_exam << "\n" : out << "";

  out << "}" << std::endl;
}

inline void
onset_sub_model::set_a_option()
{
  my_a_option = (ageon_option) (3 * my_t_option + my_m_option);
}

inline ostream& operator<< (ostream& out, const onset_sub_model& sm)
{
  out << "\n" << sm.name() << " values: \n";
  out << "Option: " << sm.option_description() << std::endl;
  
  return out;
}

}}
