#ifndef SAGE_ISNAN
#define SAGE_ISNAN

#include <cmath>

namespace SAGE
{

#ifdef __LINUX__
inline bool isnan(double d) { return ::isnan(d); }
#endif

#ifdef __SOLARIS__
inline bool isnan(double d) { return ::isnan(d); }
#endif

#ifdef __WIN32__
inline bool isnan(double d) { return std::isnan(d); }
#endif

#ifdef __MACOS__
inline bool isnan(double d) { return std::isnan(d); }
#endif

}

#endif
