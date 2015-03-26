#ifndef CONSTANTS_H
#define CONSTANTS_H

/* DESCRIPTION:
 *
 * This file contains a number of mathematical constants and
 * also some needed size parameters of the computer arithmetic.
 *
 * The default size parameters are as follows.
 *
 * For IEEE arithmetic:
 * MACHEP =  1.11022302462515654042E-16       2**-53
 * MAXLOG =  7.09782712893383996843E2         log(2**1024)
 * MINLOG = -7.08396418532264106224E2         log(2**-1022)
 * MAXNUM =  1.7976931348623158E308           2**1024

 * The global symbols for mathematical constants are
 * PI     =  3.14159265358979323846           pi
 * PIO2   =  1.57079632679489661923           pi/2
 * PIO4   =  7.85398163397448309616E-1        pi/4
 * SQRT2  =  1.41421356237309504880           sqrt(2)
 * SQRTH  =  7.07106781186547524401E-1        sqrt(2)/2
 * LOG2E  =  1.4426950408889634073599         1/log(2)
 * SQ2OPI =  7.9788456080286535587989E-1      sqrt( 2/pi )
 * LOGE2  =  6.93147180559945309417E-1        log(2)
 * LOGSQ2 =  3.46573590279972654709E-1        log(2)/2
 * THPIO4 =  2.35619449019234492885           3*pi/4
 * TWOOPI =  6.36619772367581343075535E-1     2/pi
 *
 * These lists are subject to change.
 */

#include <limits>
#include "globals/config.h"

namespace SAGE
{

template <class T>
class numeric_constants : public std::numeric_limits<T>
{
};

template<>
class numeric_constants<double> : public std::numeric_limits<double>
{
public:
  static double maxlog()              { return (max_exponent-1.0)*log(2.0);    }
  static double minlog()              { return (min_exponent-digits)*log(2.0); }
  static double pi()                  { return 3.14159265358979323846;      }
};

}

#endif
