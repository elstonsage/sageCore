#ifndef POLYGENIC_PENETRANCE_CALCULATOR_H
#include "segreg/polygenic_penetrance_calculator.h"
#endif

namespace SAGE {
namespace SEGREG {

//======================================================================
//
//  polygenic_penetrance_calculator(...)
//
//======================================================================
inline
polygenic_penetrance_calculator::polygenic_penetrance_calculator(
      const FPED::Multipedigree & ped_data,
      const model               & modp,
      bool                        use_ascertainmentp) :
      mod(modp),
      my_cont_member_calc  (ped_data,modp,use_ascertainmentp),
      my_bin_member_calc   (ped_data,modp,use_ascertainmentp),
      my_onset_member_calc (ped_data,modp,use_ascertainmentp)
{ }

//======================================================================
//
//  update()
//
//======================================================================
inline int
polygenic_penetrance_calculator::update()
{
  int err_code = pt_NONE;

  switch(mod.get_primary_trait_type())
  {
    case pt_CONTINUOUS :
      err_code = my_cont_member_calc.update();
      break;

    case pt_BINARY     :
      err_code = my_bin_member_calc.update();
      break;

    case pt_ONSET      :
      err_code = my_onset_member_calc.update();
      break;

    case pt_NONE       :
    default            :
      SAGE_internal_error();  // Should never happen!
  }

  return err_code;
}

// P(Ti|Ui,Vi) new version
inline double
polygenic_penetrance_calculator::get_polygenic_penetrance
    (penetrance_info        indiv) const
{
  double penetrance = 0.0;

  // Determine type

  switch(mod.get_primary_trait_type())
  {
    case pt_CONTINUOUS : penetrance = continuous_penetrance (indiv); break;
    case pt_BINARY     : penetrance = binary_penetrance     (indiv); break;
    case pt_ONSET      : penetrance = onset_penetrance      (indiv); break;
    case pt_NONE       :
    default            : SAGE_internal_error();    // Should never happen once onset done.
  }

  return penetrance;
}

inline double
polygenic_penetrance_calculator::continuous_penetrance
    (penetrance_info        indiv) const
{
  if(!my_cont_member_calc.is_member_valid(*indiv.member))
    return 1.0;

  const member_type& i = *indiv.member;
  genotype_index     g = indiv.genotype;
  size_t             p = indiv.polygenotype;

  double analysis_t = my_cont_member_calc.get_composite_trait  (i    );
  double exp_mean   = my_cont_member_calc.get_expected_mean    (i,g,p);
  double exp_var    = my_cont_member_calc.get_expected_variance(i,g  );
  double z          = analysis_t - exp_mean;
  double w          = exp_var;

  return calculate_continuous_penetrance(z,w, my_cont_member_calc.get_member_class(i));
}

inline double
polygenic_penetrance_calculator::binary_penetrance
    (penetrance_info        indiv)    const
{
  const member_type& i = *indiv.member;
  genotype_index     g = indiv.genotype;
  size_t             p = indiv.polygenotype;

  return my_bin_member_calc.get_penetrance(i,g,p);
}

inline double
polygenic_penetrance_calculator::onset_penetrance
    (penetrance_info        indiv)    const
{
  const member_type& i = *indiv.member;
  genotype_index     g = indiv.genotype;
  size_t             p = indiv.polygenotype;

  onset_member_calculator::member_class mem_type =
      my_onset_member_calc.get_member_class(i);

  if(mem_type == onset_member_calculator::missing)
    return 1.0;

  double onset_penetrance = 0.0;

  // There are five basic variables in the age of onset calculation.  Three
  // are calculated: expected susceptibility, expected age of onset, and
  // alpha (a term based on variance).  The other two (affection and age)
  // are not.  Note that age is either age of onset or age at exam depending
  // on the user choice.

  double suscept   = my_onset_member_calc.get_expected_susc      (i, g, p);
  double exp_onset = my_onset_member_calc.get_expected_age_onset (i, g, p);
  double alpha     = my_onset_member_calc.get_alpha_i            (i, g);

  double age = 0.0;

  if(mem_type == onset_member_calculator::age_at_exam_aff ||
     mem_type == onset_member_calculator::age_at_exam_unaff)
    age = my_onset_member_calc.get_age_exam  (i);
  else
    age = my_onset_member_calc.get_age_onset (i);

  // The age of onset penetrance can be divided into two parts, the
  // susceptibility term and the age of onset term, each of which is a
  // logistic.

  // The susceptibility is the same regardless of whether we have age of
  // exam or age of onset.

  double susc_term = exp(suscept) / (1.0 + exp(suscept));

  // The age of onset term varies depending on affection and age of onset
  // being available, but is always based on a consistent exponent to the
  // logistic.
  
  double exponent = alpha * (age - exp_onset);

  double onset_term = 0.0;

  // There are two variations to the age of onset term.  The first is for
  // all 'at' types (age_at_exam_* and age_at_onset) and the other is only
  // for age_OF_onset members.
  
  if(mem_type != onset_member_calculator::age_of_onset)
  {
    onset_term = exp(exponent) / (1.0 + exp(exponent));
  }
  else
  {
    double denom = 1.0 + exp(exponent);

    denom = denom * denom;

    onset_term = alpha * exp(exponent) / denom;
  }

  // Now we calculate the onset penetrance.  For all but age_at_exam_unaff, this
  // is just the susceptibility term by the onset term.  For unaffected, it is
  // 1 - the product of the two terms

  onset_penetrance = susc_term * onset_term;

  if(mem_type == onset_member_calculator::age_at_exam_unaff)
    onset_penetrance = 1.0 - onset_penetrance;

  return onset_penetrance;
}

}}
