#ifndef REG_PENETRANCE_CALCULATOR_H
#include "segreg/RegPenetranceCalculator.h"
#endif

namespace SAGE {
namespace SEGREG {

//======================================================================
//
//  RegPenetranceCalculator(...)
//
//======================================================================
inline
RegPenetranceCalculator::RegPenetranceCalculator
    (const FPED::Multipedigree & ped_data,
     const model               & modp,
     bool                        use_asc)
  : my_model_class   (modp.get_model_class()),
    my_residuals     (modp.resid_sub_model),
    my_member_calc   (ped_data,
                      modp,
                      use_asc),
    my_corr_adjs     (modp.get_model_class(),
                      modp.resid_sub_model,
                      MPED::mp_utilities::getMaxSibshipSize(ped_data))
{ }

inline int
RegPenetranceCalculator::update()
{
  int err;

  err = my_member_calc.update();  if(err) return err;
  err = my_corr_adjs.update();    if(err) return err;

  return 0;
}

inline PenetranceContext RegPenetranceCalculator::get_context() const
{
  return PenetranceContext(my_residuals, my_member_calc, my_corr_adjs);
}

inline size_t
RegPenetranceCalculator::get_prior_sib_count(const member_type& m) const
{
  FPED::OffspringConstIterator offspring = m.family()->offspring_begin();

  size_t i = 0;

  for( ; offspring != m.family()->offspring_end(); ++offspring, ++i)
  {
    if(&*offspring == &m) return i;
  }

  return i;
}

}}
