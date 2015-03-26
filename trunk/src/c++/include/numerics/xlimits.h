#ifndef EXTENDED_LIMITS_H
#define EXTENDED_LIMITS_H

#include <limits>
#include <cmath>
#include "globals/config.h" 

#if defined(__OSF1__) || defined(__osf__)
//#include <fp.h>
#endif

namespace xlimits
{

template <class T> T    min_value()           { return std::numeric_limits<T>::min();             }
template <class T> T    max_value()           { return std::numeric_limits<T>::max();             }
template <class T> int  digits()              { return std::numeric_limits<T>::digits;            }
template <class T> int  digits10()            { return std::numeric_limits<T>::digits10;          }
template <class T> bool is_integer()          { return std::numeric_limits<T>::is_integer;        }
template <class T> bool is_signed()           { return std::numeric_limits<T>::is_signed;         }
template <class T> bool is_exact()            { return std::numeric_limits<T>::is_exact;          }
template <class T> int  radix()               { return std::numeric_limits<T>::radix;             }
template <class T> T    epsilon()             { return std::numeric_limits<T>::epsilon();         }
template <class T> int  round_error()         { return std::numeric_limits<T>::round_error();     }
template <class T> int  min_exponent()        { return std::numeric_limits<T>::min_exponent;      }
template <class T> int  min_exponent10()      { return std::numeric_limits<T>::min_exponent10;    }
template <class T> int  max_exponent()        { return std::numeric_limits<T>::max_exponent;      }
template <class T> int  max_exponent10()      { return std::numeric_limits<T>::max_exponent10;    }
template <class T> bool has_infinity()        { return std::numeric_limits<T>::has_infinity;      }
template <class T> bool has_quiet_NaN()       { return std::numeric_limits<T>::has_quiet_NaN;     }
template <class T> bool has_signaling_NaN()   { return std::numeric_limits<T>::has_signaling_NaN; }
template <class T> T    signaling_NaN()       { return std::numeric_limits<T>::signaling_NaN();   }
template <class T> T    infinity()            { return std::numeric_limits<T>::infinity();        }
template <class T> T    quiet_NaN()           { return std::numeric_limits<T>::quiet_NaN();       }

/*

#ifndef isnan
#  if defined(ISNAN)
#    define isnan ISNAN
#  elif defined(IS_NAN)
#    define isnan IS_NAN
#  else
#    define SAGE::isnan(x) ::SAGE::isnan(x)
#  endif
#endif

#ifndef finite
#  if defined(FINITE)
#    define finite FINITE
#  elif defined(IS_FINITE)
#    define finite IS_FINITE
#  else
#    define finite(x) ::finite(x)
#  endif
#endif

*/

}

#endif
