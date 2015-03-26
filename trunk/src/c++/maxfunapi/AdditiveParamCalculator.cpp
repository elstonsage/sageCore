#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

AdditiveParamCalculator::AdditiveParamCalculator(const vector<int> & ids)
{
  my_ids = ids;
}

double 
AdditiveParamCalculator::calculateParam (const ParameterMgr * mgr) const
{
  double x = 0.0;

  for(size_t i = 0; i < my_ids.size(); ++i)
    x += mgr->operator()(my_ids[i]);

  return x;
}

} // End namespace MAXFUN
} // End namespace SAGE
