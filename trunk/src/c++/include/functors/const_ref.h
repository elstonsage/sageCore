#ifndef FUNCTORS_CONST_REF_H
#define FUNCTORS_CONST_REF_H

#include "boost/type_traits/add_const.hpp"
#include "boost/type_traits/add_reference.hpp"

namespace SAGE     {
namespace FUNCTORS {

///
/// Turns any type T into a const T &. If T is void, returns void (since a const void & won't work).
template<typename T> struct const_ref { typedef typename boost::add_const<typename boost::add_reference<T>::type>::type type; };

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
