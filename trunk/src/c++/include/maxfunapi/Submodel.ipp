#ifndef SUBMODEL_H
#include "maxfun/Submodel.h"
#endif

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//
//  getParameterMgr()
//
//=======================================================================
inline       ParameterMgr * Submodel::getParameterMgr()       { return my_info; }
inline const ParameterMgr * Submodel::getParameterMgr() const { return my_info; }

} // End MAXFUN namespace
} // End SAGE namespace
