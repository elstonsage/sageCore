#ifndef FREQ_ESTIMATES_H
#include "segreg/freq_estimates.h"
#endif

namespace SAGE
{
namespace SEGREG
{

inline size_t genotype_frequency_initial_estimates::get_model_count() const
{
  return my_model_count;
}

inline const genotype_frequency_sub_model&
    genotype_frequency_initial_estimates::get_model(size_t s) const
{
  return my_estimates[s];
}

inline genotype_frequency_initial_estimates::genotype_frequency_initial_estimates()
  : my_model_count(0)
{
}

}}
