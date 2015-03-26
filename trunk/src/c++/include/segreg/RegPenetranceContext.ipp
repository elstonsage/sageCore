#ifndef REG_PENETRANCE_CONTEXT_H
#include "segreg/RegPenetranceContext.h"
#endif

namespace SAGE {
namespace SEGREG {

//======================================================================
//
//  PenetranceContext(...)
//
//======================================================================

inline
PenetranceContextStatistics::PenetranceContextStatistics
    (const residual_correlation_sub_model& rsm,
     const continuous_member_calculator&   cmc,
     const RegPenetranceCommonAdjustments& corr_adj)
  : my_corr_sub_model  (rsm),
    my_member_calc     (cmc),
    my_corr_adjustments(corr_adj)
{ }

/// Base Constructor
///
/// \param rsm      The residual correlation sub model which supplies
///                 many of the necessary values
/// \param cmc      The continuous member calculator which supplies many
///                 of the other values needed.
/// \param corr_adj Supplies values which are common to all contexts. 
inline
PenetranceContext::PenetranceContext
    (const residual_correlation_sub_model& rsm,
     const continuous_member_calculator&   cmc,
     const RegPenetranceCommonAdjustments& corr_adj)
  : my_stats(rsm, cmc, corr_adj)
{
  clear_family_info();
}

/// Remove all the infromation about the current context.  Must be done
/// if any of the individuals within the context are to be changed. 
/// States of context individuals may be changed without clearing the
/// family info.
inline 
void
PenetranceContext::clear_family_info()
{
  my_current_family = NULL;

  my_link_ind       = NULL;
  my_link_spouse    = NULL;

  my_include_siblings = false;

  my_stats.my_mother_informative = false;
  my_stats.my_father_informative = false;

  my_stats.my_alpha_m = 0;
  my_stats.my_alpha_f = 0;

  my_stats.my_spousal_mean_adjustment = 0;
  my_stats.my_maternal_mean_adjustment = 0;
  my_stats.my_paternal_mean_adjustment = 0;

  my_stats.my_sibling_mean_adjustments.clear();
}

//======================================================================
//
//  Penetrance initialization functions start here.......
//
//======================================================================

/// Set the current nuclear family.  The nuclear family in question is
/// the family which includes the link individual as a child.
///
/// \param family       The nuclear family to add to the context.
///
/// \param include_sibs In class D models, sibling correlations are
///                     sometimes added.  This can be controlled by the
///                     include_sibs flag.  Has no effect with class A models.
inline
void
PenetranceContext::set_nuclear_family(const family_type& family, bool include_sibs)
{
  // Set basic variables

  my_current_family = &family;

  my_include_siblings = include_sibs;

  assert(my_current_family->get_mother());
  
  my_stats.calculate_family_statistics(*this);
}

/// Set the couple (used for founder penetrances and SL3) for which
/// penetrances are calculated
inline
void
PenetranceContext::set_link_couple
    (const member_type& ind, const member_type& spouse)
{
  my_link_ind    = &ind;
  my_link_spouse = &spouse;
}

/// Set the state of the link spouse.  Because of the way the anterior
/// function works, we only need to set this once for each call to SL3.
inline
void
PenetranceContext::set_link_spouse_state
    (const TypeDescription::State& spouse_state)
{
  my_stats.calculate_spousal_statistics(*this,spouse_state);
}


/// Set the mother's state
///
inline
void
PenetranceContext::set_mother_state
    (const TypeDescription::State& mother_state)
{
  my_mother_state = mother_state;

  my_stats.calculate_maternal_statistics(*this,mother_state);
}

/// Set the father's state
///
inline
void
PenetranceContext::set_father_state
    (const TypeDescription::State& father_state)
{
  my_father_state = father_state;

  my_stats.calculate_paternal_statistics(*this,father_state);
}

    
/// Returns true if the context includes a nuclear family component
///
inline
bool
PenetranceContext::includes_family() const
{
  return my_current_family != NULL;
}

/// Returns true if the context includes a spouse component
///
inline
bool
PenetranceContext::includes_spouse() const
{
  return my_link_ind != NULL;
}

/// Returns true if the context includes a sibling component
///
inline
bool
PenetranceContext::includes_siblings() const
{
  return my_include_siblings;
}

/// Returns \c true if spouse adjustments are needed for
/// individual \em i, \c false otherwise.
///
/// \param i individual who may require spousal adjustments to the penetrance
inline
bool
PenetranceContext::incorporate_spousal_adjustments(const member_type& i) const
{
  return (&i == my_link_ind);
}

inline
const PenetranceContext::member_type&
PenetranceContext::get_individual  () const
{
  return *my_link_ind;
}

inline
const PenetranceContext::member_type& 
PenetranceContext::get_mother      () const
{
  return *my_current_family->get_mother();
}

inline
const PenetranceContext::member_type& 
PenetranceContext::get_father      () const
{
  return *my_current_family->get_father();
}

inline
const PenetranceContext::member_type& 
PenetranceContext::get_link_spouse () const
{
  return *my_link_spouse;
}

inline
const TypeDescription::State&
PenetranceContext::get_mother_state() const
{
  return my_mother_state;
}

inline
const TypeDescription::State&
PenetranceContext::get_father_state() const
{
  return my_father_state;
}

inline
const PenetranceContext::family_type&
PenetranceContext::get_family() const
{
  return *my_current_family;
}

/// Returns the adjustment of the mean due to parents
///
inline
double
PenetranceContextStatistics::get_parental_mean_adjustment() const
{
  return my_maternal_mean_adjustment + my_paternal_mean_adjustment;
}

/// Returns the adjustment of the mean due to the spouse
///
inline
double
PenetranceContextStatistics::get_spousal_mean_adjustment() const
{
  return my_spousal_mean_adjustment;
}

/// Get the sibling mean adjustment for this individual
///
inline
double
PenetranceContextStatistics::get_sibling_mean_adjustment
    (size_t                   prior_sib_count,
     const PenetranceContext& context) const
{
  if(!context.includes_siblings())
    return 0.0;

  if(my_sibling_mean_adjustments.empty())
    calculate_sibling_mean_adjustments(context);

  return my_sibling_mean_adjustments[prior_sib_count];
}


/// Returns \c true if mother is informative, \c false otherwise.
///
inline
bool
PenetranceContextStatistics::get_mother_informativity() const
{
  return my_mother_informative;
}

/// Returns \c true if father is informative, \c false otherwise.
///
inline
bool
PenetranceContextStatistics::get_father_informativity() const
{
  return my_father_informative;
}

}}
