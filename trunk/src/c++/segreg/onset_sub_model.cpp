//============================================================================
// File:      onset_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/9/01 - created.                            djb
//                                                                          
// Notes:     implementation of onset sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/onset_sub_model.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{

const std::string  ONSET_NAME            = "Age of onset";

//============================================================================
// IMPLEMENTATION:  onset_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
bool
onset_sub_model::set(type_option t_opt, multi_option m_opt, const RPED::RefMultiPedigree* mp, 
                     const string& status, const string& onset, const string& exam)
{
//  if(my_in_use)
//  {
//    return false;
//  }

  my_t_option = t_opt;
  my_m_option = m_opt;

  bool  return_value = false;
  
  if(trait_given(status, "status") && trait_given(onset, "age_onset") && 
     trait_given(exam, "age_exam"))
  {
    if(trait_valid(status, RPED::RefTraitInfo::binary_trait, mp)     && 
       trait_valid(onset, RPED::RefTraitInfo::continuous_trait, mp)  && 
       trait_valid(exam, RPED::RefTraitInfo::continuous_trait, mp)      )
    {
      my_affection_status = status;
      my_age_of_onset = onset;
      my_age_at_exam = exam;
      return_value = true;
    }
  }

  set_a_option();
  
  return return_value;
}

//
// ===============  Ancillary functions
//
bool 
onset_sub_model::trait_valid(const string& trait_name, 
                             RPED::RefTraitInfo::trait_t req_type, const RPED::RefMultiPedigree* mp)
{
  bool  return_value = false;

  if(mp->info().trait_exists(trait_name))
  {
    if(req_type == mp->info().get_trait_type(trait_name))
    {
      return_value = true;
    }
    else
    {
      my_errors << priority(critical) << "Trait '" << trait_name
                << "' specified in onset sub-block is the not correct type.  "
                << "Skipping analysis ..." << endl;
    }  
  }
  else
  {
    my_errors << priority(critical) << "Could not find trait '" << trait_name
              << "' specified in onset sub-block.  Skipping analysis ..." << endl; 
  }
  
  return return_value;
}

bool 
onset_sub_model::trait_given(const string& trait_name, const string& param)
{
  if(! trait_name.empty())
  {
    return true;  
  }
  else
  {
    my_errors << priority(critical) << "No trait given for " << param
              << " in onset sub-block.  Skipping analysis ..." << endl; 
    return false;
  }
}

}
}
