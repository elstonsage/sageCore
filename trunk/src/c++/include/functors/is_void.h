#ifndef FUNCTORS_IS_VOID_H
#define FUNCTORS_IS_VOID_H

#include "boost/type_traits/is_same.hpp"
#include "boost/mpl/bool.hpp"

namespace SAGE     {
namespace FUNCTORS {

/// Indicates whether or not R is of type 'void' (via the boolean 'value').
template<typename R> struct is_void
{
  static const bool value = boost::is_same<R, void>::type::value;

  typedef typename boost::mpl::bool_<value>::type type;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
