#ifndef TDTEX_NUMERICS_H
#define TDTEX_NUMERICS_H

#include "numerics/cephes.h"

namespace SAGE  {
namespace TDTEX {

/// Compute factorial n! by conventional means if n < 32 or by gamma function if larger.
inline double fact(int n)
{
  static int top = 1;
  static double a[500] = { 1.0, 1.0 };
  if(n>=500 || (!a[n] && n>32) )
  {
    double res = exp( lgam(n+1.0) );
    if(n < 500) a[n] = res;
    return res;
  }
  if(n<top) return a[n];
  for( ; top <= n ; ++top)
    a[top+1]=a[top-1]*top;
  return a[top];
}

/// Compute factorial log(n!) by gamma function, with cache
inline double log_fact(int n)
{
  static double la[500];
  double res;

  if( n >= 500 || !la[n])
  {
    res = lgam(n+1.0);
    if(n < 500) la[n] = res;
    return res;
  }
  return la[n];
}

} // End namespace TDTEX
} // End namespace SAGE

#endif
