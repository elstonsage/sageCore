#ifndef FREQ_ESTIMATES_H
#define FREQ_ESTIMATES_H

#include <string>
#include "rped/rped.h"
#include "segreg/freq_sub_model.h"
#include "segreg/mean_split.h"

namespace SAGE
{
namespace SEGREG
{

class genotype_frequency_initial_estimates
{
public:

  typedef genotype_frequency_sub_model freq_sub_model;

  genotype_frequency_initial_estimates(const mean_split& splitter, const freq_sub_model& fsm);

  size_t get_model_count() const;

  const freq_sub_model& get_model(size_t) const;

 static bool cont;
private:

  size_t         my_model_count;

  freq_sub_model my_estimates[6]; // ++++ needed for new binary initial values

  //lint -e(1704) This is private to prevent default construction
  genotype_frequency_initial_estimates();
};


} // end SEGREG namespace
} // end SAGE namespace

#include "segreg/freq_estimates.ipp"

#endif
