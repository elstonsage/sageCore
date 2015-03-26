#include "segreg/RegPenetranceContext.h"

namespace SAGE {
namespace SEGREG {

void
PenetranceContextStatistics::calculate_family_statistics
    (const PenetranceContext&      context)
{
  // Determine parental informativity

  my_mother_informative = my_member_calc.is_member_valid(context.get_mother());
  my_father_informative = my_member_calc.is_member_valid(context.get_father());

  // Lookup alpha_m and alpha_f (see eqs. 53 and 54)

  my_alpha_m = my_corr_sub_model.alpha_mother(my_father_informative,my_mother_informative);
  my_alpha_f = my_corr_sub_model.alpha_father(my_father_informative,my_mother_informative);
}

void
PenetranceContextStatistics::calculate_spousal_statistics
    (const PenetranceContext&      context,
     const TypeDescription::State& spouse_state)
{
  double Pfm        = my_corr_sub_model.father_mother_correlation();
  double std_spouse = my_member_calc.get_standardization(context.get_link_spouse(), spouse_state.get_index());

  my_spousal_mean_adjustment = - Pfm * std_spouse;
}

void
PenetranceContextStatistics::calculate_maternal_statistics
    (const PenetranceContext&      context,
     const TypeDescription::State& mstate)
{
  double std_mother = my_member_calc.get_standardization(context.get_mother(), mstate.get_index());

  my_maternal_mean_adjustment = - my_alpha_m * std_mother;

  my_sibling_mean_adjustments.clear();
}

void
PenetranceContextStatistics::calculate_paternal_statistics
    (const PenetranceContext&      context,
     const TypeDescription::State& fstate)
{
  double std_father = my_member_calc.get_standardization(context.get_father(), fstate.get_index());

  my_paternal_mean_adjustment = - my_alpha_f * std_father;

  my_sibling_mean_adjustments.clear();
}

/// Calculates the sibling mean adjustments.
///
/// This function is const because it is called only when sibling means
/// are requested by the get_sibling_mean_adjustment function.  It
/// modifies the mutable storage to calculate all sibling mean
/// adjustments (not just the specific one requested).
void
PenetranceContextStatistics::calculate_sibling_mean_adjustments
    (const PenetranceContext& context) const
{
  my_sibling_mean_adjustments.resize(context.my_current_family->offspring_count());

  // The estimated standardization, \hat{d}_i (eq. 65) is calculated by the
  // continuous_member_calculator.  We wish to calculate the negative sums
  // of the \hat{d}_i values for all prior sibs for each sibling.

  // The first sibling has no prior sibs, so the sum is always 0.0

  my_sibling_mean_adjustments[0] = 0.0;

  // Each later sib's negative sum is simply the previous sib's sum minus
  // the previous sib's \hat{d}_i

  FPED::OffspringConstIterator offspring =
      context.my_current_family->offspring_begin();

  size_t i = 1;

  for( ; i < context.my_current_family->offspring_count(); ++offspring, ++i)
  {
    double est_stan = 
        my_member_calc.get_estimated_standardization(*offspring,
                                                     context.my_mother_state.get_index(),
                                                     context.my_father_state.get_index());

    my_sibling_mean_adjustments[i] = my_sibling_mean_adjustments[i-1] - est_stan;
  }

  // Adjust sums by betas (see eqs. 63 and 64)

  i = 1;

  for( ; i < context.my_current_family->offspring_count(); ++i)
  {
    my_sibling_mean_adjustments[i] *= 
        my_corr_adjustments.get_prior_sib_coefficient(my_mother_informative,
                                                      my_father_informative, i);
  }
}

}}
