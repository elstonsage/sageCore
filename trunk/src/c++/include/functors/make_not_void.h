#ifndef FUNCTORS_MAKE_NOT_VOID_H
#define FUNCTORS_MAKE_NOT_VOID_H

#include "boost/type_traits/is_same.hpp"

namespace SAGE     {
namespace FUNCTORS {

struct NotVoidType { };

/// If T is void, makes the typedef 'type' point to NotVoidType. Otherwise, makes the typedef
/// 'type' point to T.
///
template<typename T> struct make_not_void { typedef T type; };

template<> struct make_not_void<void> { typedef NotVoidType type; };

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
