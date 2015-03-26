#ifndef MLM_FRA_UTIL_H
#define MLM_FRA_UTIL_H

#include "segreg/member_calculator.h"
#include "pedcalc/fam_resid_adj.h"

namespace SAGE {
namespace SEGREG {

struct ResidualGetter
{
  public:
  
    ResidualGetter(const residual_correlation_sub_model& resid)
      : my_resids(resid) { }
      
    double operator() (PED_CALC::FraResidType tp)
      { return my_resids.correlation((residual_correlation_sub_model::corr) tp); }
      
  private:
  
    const residual_correlation_sub_model& my_resids;
};

}
}

#endif
