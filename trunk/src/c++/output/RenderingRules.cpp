#include "output/RenderingRules.h"

namespace SAGE   {
namespace OUTPUT {

const RenderingRules default_rules;

//===================================
//
//  exceedsThreshold(...)
//
//===================================
bool 
RenderingRules::exceedsThreshold(double v) const
{
  for(HasLowerThreshold::ThresholdSet::const_iterator i = getLowerThresholds().begin(); i != getLowerThresholds().end(); ++i)
    if(v < *i)
      return true;

  if(getUpperThresholds().size())
  {
    HasUpperThreshold::ThresholdSet::const_iterator i = getUpperThresholds().end();

    while(1)
    {
      i--;

      if(v > *i)
        return true;
      
      if(i == getUpperThresholds().begin())
        break;
    }
  }

  return false;
}

} // End namespace OUTPUT
} // End namespace SAGE
