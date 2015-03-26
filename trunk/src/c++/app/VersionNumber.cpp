#include "app/VersionNumber.h"

#if !defined(SAGE_MAIN_VERSION) || !defined(SAGE_SUB_VERSION) || !defined(SAGE_MICRO_VERSION) || !defined(SAGE_BETA_VERSION)
#error You must declare all four components of this version for app to compile correctly (ie: gcc -DSAGE_MAIN_VERSION=1 -DSAGE_SUB_VERSION=2 -DSAGE_MICRO_VERSION=3 -DSAGE_BETA_VERSION=0)
#endif

namespace SAGE {
namespace APP  {

VersionNumber version(SAGE_MAIN_VERSION, SAGE_SUB_VERSION, SAGE_MICRO_VERSION, SAGE_BETA_VERSION);

} // End namespace APP
} // End namespace SAGE

