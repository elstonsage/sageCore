#ifndef CONTINUOUS_PENETRANCE_FUNCTION_H
#include "segreg/ContinuousPenetranceFunction.h"
#endif

namespace SAGE {
namespace SEGREG {

inline
double calculate_continuous_penetrance
    (double                               z,
     double                               w,
     member_calculator_base::member_class mt)
{
  // Determine the individual type
  if(SAGE::isnan(w) || SAGE::isnan(z))
  {
    cout << "NaN" << w << ' ' << z << endl;
    exit(1);
    return 0.0;
  }

  // Small restriction on w, to hopefully make things work better.

  if(w < 0.01) w = 0.01;

  // Determine the type and perform the appropriate calculation

  switch(mt)
  {
    case member_calculator_base::actual :
      return NUMERICS::normal_pdf(z,w);

    case member_calculator_base::missing :
      return 1.0;

    case member_calculator_base::gte_thresh :
      return normal_cdf(-z / sqrt(w));

    case member_calculator_base::lte_thresh :
      return normal_cdf(z / sqrt(w));
      
    case member_calculator_base::age_at_onset :
      // Never happens.
      break;
  }

  return 1.0;
}

}
}
