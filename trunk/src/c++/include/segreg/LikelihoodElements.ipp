#ifndef LIKELIHOOD_ELEMENTS_H
#include "segreg/LikelihoodElements.h"
#endif

namespace SAGE {
namespace SEGREG {

inline
  LikelihoodElements::LikelihoodElements
    (const FPED::Multipedigree& mped,
     const model&               mdl,
     bool                       use_ascertainment)
  : my_freq_function(boost::bind(&genotype_frequency_sub_model::prob,
                                 boost::ref(mdl.freq_sub_model), _1)),
    my_transm_function(boost::bind(&transmission_sub_model::prob,
                                 boost::ref(mdl.transm_sub_model), _1, _2, _3)),
    my_pen_calc(mped,mdl,use_ascertainment)
{
  build_type_description();
}

inline int
  LikelihoodElements::update()
{
  int err;

  err = my_pen_calc.update();  if(err) return err;

  return 0;
}

inline const TypeDescription&
  LikelihoodElements::get_type_description() const
{
  return my_type_description;
}

inline double
  LikelihoodElements::get_frequency
    (const TypeDescription::State& state) const
{
  return my_freq_function(state.get_index());
}

inline double
  LikelihoodElements::get_frequency
    (const TypeDescription::StateIterator& state) const
{
  return get_frequency(*state);
}

inline double 
  LikelihoodElements::get_transmission
    (const TypeDescription::State& istate,
     const TypeDescription::State& mstate,
     const TypeDescription::State& fstate)    const
{
  return my_transm_function(istate.get_index(), mstate.get_index(), fstate.get_index());
}

inline double 
  LikelihoodElements::get_transmission
    (const TypeDescription::StateIterator& istate,
     const TypeDescription::StateIterator& mstate,
     const TypeDescription::StateIterator& fstate)    const
{
  return get_transmission(*istate, *mstate, *fstate);
}

inline PenetranceContext
  LikelihoodElements::get_penetrance_context() const
{
  return my_pen_calc.get_context();
}

inline double
  LikelihoodElements::get_penetrance
    (const FPED::Member&                ind,
     const TypeDescription::State&      state) const
{
  return my_pen_calc.get_penetrance(
            RegPenetranceCalculator::penetrance_info(ind, state.get_index()));
}

inline double
  LikelihoodElements::get_penetrance
    (const FPED::Member&                        ind,
     const TypeDescription::StateIterator&      state) const
{
  return get_penetrance(ind, *state);
}

inline double
  LikelihoodElements::get_penetrance
    (const FPED::Member&                        ind,
     const TypeDescription::State&              state,
     const PenetranceContext&                   context)  const
{
  return my_pen_calc.get_penetrance(
            RegPenetranceCalculator::penetrance_info(ind, state.get_index()),
            context);
}

inline double
  LikelihoodElements::get_penetrance
    (const FPED::Member&                     ind,
     const TypeDescription::StateIterator&   state,
     const PenetranceContext&                context)  const
{
  return get_penetrance(ind, *state, context);
}

inline void
  LikelihoodElements::build_type_description()
{
  my_type_description.set_name("Main Type");
  
  my_type_description.add_state("AA");
  my_type_description.add_state("AB");
  my_type_description.add_state("BB");
}
  
}
}

