//============================================================================
// File:      ascertainment_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/12/01 - created.                            djb
//                                                                          
// Notes:     implementation of ascertainment sub_model class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/ascertainment_sub_model.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{

const std::string  ASCERTAINMENT_NAME    = "Ascertainment";
const double       ASCER_INCLUDE_DEFAULT_VALUE = 1;
const double       ASCER_INDIC_THRESH_DEFAULT_VALUE = 0;

//============================================================================
// IMPLEMENTATION:  ascertainment_sub_model
//============================================================================

/// Calculates the thresholds, if needed
bool ascertainment_sub_model::calculate_thresholds
    (const FPED::Multipedigree& rmp, size_t primary_trait)
{
  // Determine which, if any, of the thresholds need calculated:

  switch(v_option())
  {
    case not_specified :
    case actual        :
    case onset         :
      break;

    case gte_thresh    :
      if(!my_th_specified)
        return calculate_threshold_high(rmp, primary_trait);
      break;

    case lte_thresh    :
      if(!my_tl_specified)
        return calculate_threshold_low(rmp, primary_trait);
      break;

    case thresh_indic  :
      if(!my_tl_specified)
        if(!calculate_threshold_low(rmp, primary_trait))
          return false;
      if(!my_th_specified)
        if(!calculate_threshold_high(rmp, primary_trait))
          return false;
      
      // If either threshold isn't computed, we're ok (no one in the data
      // set uses that threshold).  Otherwise, we must verify that the 
      // Tlow < Thigh.

      return SAGE::isnan(my_thresh_low) || SAGE::isnan(my_thresh_high) ||
             my_thresh_low < my_thresh_high;
  }

  return true;
}

/// Calculates the threshold_high
bool ascertainment_sub_model::calculate_threshold_high
    (const FPED::Multipedigree& rmp, size_t primary_trait)
{
  my_thresh_high = POSITIVE_INF;

  FPED::PedigreeConstIterator ped  = rmp.pedigree_begin();

  for( ; ped != rmp.pedigree_end(); ++ped)  
  {
    FPED::MemberConstIterator mem = ped->member_begin();

    for( ; mem != ped->member_end(); ++mem)
    {
      // Determine if this individual uses the high threshold.
      if(ind_uses_thresh_high(&*mem))
      {
        double t = ped->info().trait(mem->index(), primary_trait);

        if(t < my_thresh_high)
          my_thresh_high = t;
      }
    }
  }

  if(finite(my_thresh_high))
    my_th_specified = true;
  else
    my_thresh_high = QNAN;

  return true;
}

/// Calculates the threshold_low
bool ascertainment_sub_model::calculate_threshold_low
    (const FPED::Multipedigree& rmp, size_t primary_trait)
{
  my_thresh_low = NEGATIVE_INF;

  FPED::PedigreeConstIterator ped  = rmp.pedigree_begin();

  for( ; ped != rmp.pedigree_end(); ++ped)  
  {
    FPED::MemberConstIterator mem = ped->member_begin();

    for( ; mem != ped->member_end(); ++mem)
    {
      // Determine if this individual uses the high threshold.
      if(ind_uses_thresh_low(&*mem))
      {
        double t = ped->info().trait(mem->index(), primary_trait);

        if(t > my_thresh_low)
          my_thresh_low = t;
      }
    }
  }

  if(finite(my_thresh_low))
    my_tl_specified = true;
  else
    my_thresh_low = QNAN;

  return true;
}

//
// ===============  Setting the sub-model
//
bool  
ascertainment_sub_model::set(s_sm_option so, v_sm_option vo,
                             const RPED::MultiPedigree* mp, const string& pi, const vector<double>& incl, 
                             const string& ti, double thrsh, double thrsh_high, double thrsh_low, double it)
{
  my_s_option = none;
  my_v_option = actual;
  my_psf_indicator = string();
  my_includes.resize(0);
  my_includes.push_back(ASCER_INCLUDE_DEFAULT_VALUE);
  my_thresh_indicator = string();
  my_thresh = QNAN;
  my_thresh_high = QNAN;
  my_thresh_low = QNAN;
  my_indic_thresh = ASCER_INDIC_THRESH_DEFAULT_VALUE;
  my_t_specified = false;
  my_th_specified = false;
  my_tl_specified = false;

  bool  return_value = false;
  switch(so)
  {
    case none:
      return_value = set_none(vo, mp, pi, incl, ti, thrsh, thrsh_high, thrsh_low, it);
      break;
    case founders:
      return_value = set_founders(vo, mp, pi, incl, ti, thrsh, thrsh_high, thrsh_low, it); 
      break;
    case psf:
      return_value = set_psf(vo, mp, pi, incl, ti, thrsh, thrsh_high, thrsh_low, it); 
      break;
    case founders_plus_psf:
      return_value = set_founders_plus_psf(vo, mp, pi, incl, ti, thrsh, thrsh_high, thrsh_low, it); 
      break;
    default:
      assert(false);
  }
  
  return return_value;
}

// - No parameters relevant for this option.
//
bool  
ascertainment_sub_model::set_none(v_sm_option vo, const RPED::MultiPedigree* , 
                                  const string& pi, const vector<double>& incl, const string& ti, 
                                  double thrsh, double thrsh_high, 
                                  double thrsh_low, double it)
{
  if(vo != not_specified || pi.size() || ti.size() || (! SAGE::isnan(thrsh)) || 
    (! SAGE::isnan(thrsh_high)) || (! SAGE::isnan(thrsh_low)) || incl.size() || (! SAGE::isnan(it)))
  {
    my_errors << priority(error) << "Unnecessary parameter(s) specified in ascertainment sub-block "
              << "with cond_set equal to 'none'.  Ignoring ..." << endl;
  }

  return true;
}
             
bool  
ascertainment_sub_model::set_founders(v_sm_option vo, const RPED::MultiPedigree* mp, 
                                      const string& pi, const vector<double>& incl, const string& ti, 
                                      double thrsh, double thrsh_high, 
                                      double thrsh_low, double it)
{
  my_s_option = founders;
  
  if(pi.size())
  {
    my_errors << priority(error) << "psf_indic specified in ascertainment sub-block "
              << "with cond_set equal to 'founders'.  Ignoring ..." << endl; 
  }
  
  if(incl.size())
  {
    my_errors << priority(error) << "psf_indic attribute, include, specified in ascertainment sub-block "
              << "with cond_set equal to 'founders'.  Ignoring ..." << endl; 
  }
  
  
  bool  return_value = false;
  switch(vo)
  {
    case actual:
      return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
      break;
    case gte_thresh:
      return_value = set_gte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
      break;
    case lte_thresh:
      return_value = set_lte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
      break;
    case thresh_indic:
      return_value = set_thresh_indic(mp, ti, thrsh, thrsh_high, thrsh_low, it);
      break;
    case onset:
      return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
      my_v_option = onset;
      break;

    case not_specified:
    default:
      assert(false);
  }
  
  return return_value;
}

bool  
ascertainment_sub_model::set_psf(v_sm_option vo, const RPED::MultiPedigree* mp, 
                                 const string& pi, const vector<double>& incl, const string& ti, 
                                 double thrsh, double thrsh_high, 
                                 double thrsh_low, double it)
{
  my_s_option = psf;
  
  if(trait_valid(pi, RPED::RefTraitInfo::binary_trait, RPED::RefTraitInfo::continuous_trait, mp))
  {
    my_psf_indicator = pi;
    set_includes(incl);

    bool  return_value = false;
    switch(vo)
    {
      case actual:
        return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case gte_thresh:
        return_value = set_gte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case lte_thresh:
        return_value = set_lte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case thresh_indic:
        return_value = set_thresh_indic(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case onset:
        return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        my_v_option = onset;
        break;

      case not_specified:
      default:
        assert(false);
    }
    
    return return_value;
  }
  else
  {
    return false;
  }
}

bool  
ascertainment_sub_model::set_founders_plus_psf(v_sm_option vo, const RPED::MultiPedigree* mp, 
                                               const string& pi, const vector<double>& incl, 
                                               const string& ti, double thrsh, double thrsh_high, 
                                               double thrsh_low, double it)
{
  my_s_option = founders_plus_psf;

  if(trait_valid(pi, RPED::RefTraitInfo::binary_trait, RPED::RefTraitInfo::continuous_trait, mp))
  {
    my_psf_indicator = pi;
    set_includes(incl); 
      
    bool  return_value = false;
    switch(vo)
    {
      case actual:
        return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case gte_thresh:
        return_value = set_gte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case lte_thresh:
        return_value = set_lte_thresh(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case thresh_indic:
        return_value = set_thresh_indic(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        break;
      case onset:
        return_value = set_actual(mp, ti, thrsh, thrsh_high, thrsh_low, it);
        my_v_option = onset;
        break;
      case not_specified:
      default:
        assert(false);
    }
    
    return return_value;
  }
  else
  {
    return false;
  }
}


//
// ===============  Ancillary functions
//
bool  
ascertainment_sub_model::set_actual(const RPED::MultiPedigree* , const string& ti, 
                                    double thrsh, double thrsh_high, 
                                    double thrsh_low, double it)
{
  my_v_option = actual;
  
  if(ti.size() || (! SAGE::isnan(thrsh)) || 
    (! SAGE::isnan(thrsh_high)) || (! SAGE::isnan(thrsh_low)) || (! SAGE::isnan(it)))
  {
    my_errors << priority(error) << "Unnecessary parameter(s) specified in ascertainment "
              << "sub-block with cond_val equal to 'actual'.  Ignoring ...";
  }
  
  return true;
}
                                 
bool  
ascertainment_sub_model::set_gte_thresh(const RPED::MultiPedigree* , const string& ti, 
                                        double thrsh, double thrsh_high, 
                                        double thrsh_low, double it)
{
  my_v_option = gte_thresh;
  
  if(! SAGE::isnan(thrsh))
  {
    my_thresh = thrsh;
    
    // - SCR 138.
    //
    my_thresh_low = thrsh;
    my_thresh_high = thrsh;
    
    my_t_specified = true;
  }

  if(ti.size() || (! SAGE::isnan(thrsh_high)) || (! SAGE::isnan(thrsh_low)) ||
                  (! SAGE::isnan(it)))
  {
    my_errors << priority(error) << "Unnecessary parameter(s) specified in ascertainment "
              << "sub-block with thresh equal to 'gte_thresh'.  Ignoring ...";
  }

  return true;
}

bool  
ascertainment_sub_model::set_lte_thresh(const RPED::MultiPedigree* , const string& ti, 
                                        double thrsh, double thrsh_high, 
                                        double thrsh_low, double it)
{
  my_v_option = lte_thresh;
  
  if(! SAGE::isnan(thrsh))
  {
    my_thresh = thrsh;
    
    // - SCR 138.
    //
    my_thresh_low = thrsh;
    my_thresh_high = thrsh;
    
    my_t_specified = true;
  }

  if(ti.size() || (! SAGE::isnan(thrsh_high)) || (! SAGE::isnan(thrsh_low)) ||
                  (! SAGE::isnan(it)))
  {
    my_errors << priority(error) << "Unnecessary parameter(s) specified in ascertainment "
              << "sub-block with thresh equal to 'gte_thresh'.  Ignoring ...";
  }

  return true;
}
                                 
bool  
ascertainment_sub_model::set_thresh_indic(const RPED::MultiPedigree* mp, const string& ti, 
                                          double thrsh, double thrsh_high, 
                                          double thrsh_low, double it)
{
  my_v_option = thresh_indic;
  
  if(trait_valid(ti, RPED::RefTraitInfo::continuous_trait, RPED::RefTraitInfo::continuous_trait, mp))
  {
    my_thresh_indicator = ti;
    
    if(! SAGE::isnan(thrsh_high) && ! SAGE::isnan(thrsh_low))
    {
      if(thrsh_low >= thrsh_high)
      {
        my_errors << priority(critical) << "thresh_indic_low is greater than or "
                  << "equal to thresh_indic_high in ascertainment sub-block.  "
                  << "Skipping analysis ..." << endl;
        return false;
      }
    }
  
    if(! SAGE::isnan(thrsh_high))
    {
      my_thresh_high = thrsh_high;
      my_th_specified = true;
    }
    
    if(! SAGE::isnan(thrsh_low))
    {
      my_thresh_low = thrsh_low;
      my_tl_specified = true;
    }
    
    if(! SAGE::isnan(it))
    {
      my_indic_thresh = it;
    }
     
    if(! SAGE::isnan(thrsh))
    {
      my_errors << priority(error) << "cond_val attribute, thresh, specified in ascertainment "
                << "sub-block with cond_val equal to 'thresh_indic'.  Ignoring attribute, thresh ...";
    }
    
    return true;
  }
  else
  {
    return false;
  }
}

void
ascertainment_sub_model::set_includes(const vector<double>& incl)
{
  my_includes.resize(0);

  vector<double>::const_iterator  iter;
  for(iter = incl.begin(); iter != incl.end(); iter++)
  {
    if(! SAGE::isnan(*iter))
    {
      vector<double>::const_iterator  result;
      result = find(my_includes.begin(), my_includes.end(), *iter);
      if(result == my_includes.end())
      {
        my_includes.push_back(*iter);
      }
    }
  }
  
  if(! my_includes.size())
  {
    my_includes.push_back(ASCER_INCLUDE_DEFAULT_VALUE);
  }
}

bool 
ascertainment_sub_model::trait_valid(const string& trait_name, RPED::RefTraitInfo::trait_t req_type1,
                                     RPED::RefTraitInfo::trait_t req_type2, const RPED::MultiPedigree* mp)
{
  bool  return_value = false;

  // - Per rce, there are no restrictions on using the same trait
  //   for more than one purpose in the ascertainment block or on
  //   using the primary trait in the ascertainment block.
  //
  if(trait_name.size())
  {
    if(mp->info().trait_exists(trait_name))
    {
      if(mp->info().get_trait_type(trait_name) == req_type1 ||
         mp->info().get_trait_type(trait_name) == req_type2   )
      {
        return_value = true;
      }
      else
      {
        my_errors << priority(critical) << "Trait '" << trait_name
                  << "' specified for ascertainment sub-block parameter, "
                  << "psf_indic or thresh_indic, is not the correct type.  "
                  << "Skipping analysis ..." << endl;
      }  
    }
    else
    {
      my_errors << priority(critical) << "Could not find trait '" << trait_name
                << "' specified for ascertainment sub-block parameter, psf_indic or thresh_indic."
                << "  Skipping analysis ..." << endl; 
    }
  }
  else
  {
    my_errors << priority(critical) << "Required trait not specified for ascertainment sub-block parameter, "
              << "psf_indic or thresh_indic.  Skipping analysis ..." << endl;
  }
  
  return return_value;
}

}
}
