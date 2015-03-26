#include "math.h"

#include "numerics/mt.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

int geometric(float  q, SAGE::MersenneTwister& r)
{
  if(q == 1) return 0;

  double x = r.uniform_real();

  while(x == 0 || x == 1)
    x=r.uniform_real();

  return  (int) (log(1-x) / log(q) - 1);
}
