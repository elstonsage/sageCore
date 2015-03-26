#include "maxfunapi/DivideByParamTransformer.h"

namespace SAGE   {
namespace MAXFUN {

double DivideByParamTransformer::transformToMaximizationScale(const ParameterMgr * mgr, double val) const
{
  return val * mgr->operator()(my_param_idx);
}

double DivideByParamTransformer::transformToReportedScale(const ParameterMgr * mgr, double val) const
{
  return val / mgr->operator()(my_param_idx);
}

} // End namespace MAXFUN
} // End namespace SAGE
