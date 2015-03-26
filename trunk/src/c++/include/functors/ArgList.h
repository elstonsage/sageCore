#ifndef FUNCTORS_ARGLIST_H
#define FUNCTORS_ARGLIST_H

#include "boost/type_traits/is_same.hpp"
#include "boost/mpl/equal_to.hpp"
#include "boost/mpl/greater_equal.hpp"
#include "boost/mpl/size.hpp"
#include "boost/static_assert.hpp"
#include "functors/at_idx.h"

namespace SAGE     {
namespace FUNCTORS {

struct null_type { };

template<typename ARGS, int IDX, bool AVAILABLE> struct Arg;

template<typename ARGS, int IDX> struct Arg<ARGS, IDX, true>  { typedef typename at_idx<ARGS, IDX>::type type; };
template<typename ARGS, int IDX> struct Arg<ARGS, IDX, false> { typedef null_type type; };

/// \brief Provides helper functionality for accessing the arguments in an argument list
template<typename ARGS>
struct ArgList
{
  static const int arity = boost::mpl::size<ARGS>::value;

  typedef typename Arg<ARGS, 0, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<1> >::value>::type arg1_type;
  typedef typename Arg<ARGS, 1, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<2> >::value>::type arg2_type;
  typedef typename Arg<ARGS, 2, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<3> >::value>::type arg3_type;
  typedef typename Arg<ARGS, 3, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<4> >::value>::type arg4_type;
  typedef typename Arg<ARGS, 4, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<5> >::value>::type arg5_type;
  typedef typename Arg<ARGS, 5, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<6> >::value>::type arg6_type;
  typedef typename Arg<ARGS, 6, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<7> >::value>::type arg7_type;
  typedef typename Arg<ARGS, 7, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<8> >::value>::type arg8_type;
  typedef typename Arg<ARGS, 8, boost::mpl::greater_equal<boost::mpl::size<ARGS>, boost::mpl::int_<9> >::value>::type arg9_type;
};

/// \brief Performs a compile-time assertion that the given sequence is of the given size
/// 
/// SEQ A sequence
///
/// SIZE The required size of the sequence
///
template<typename SEQ, int SIZE>
void assertSeqSize()
{
  BOOST_STATIC_ASSERT(( boost::mpl::equal_to<typename boost::mpl::size<SEQ>, boost::mpl::int_<SIZE> >::value )); // "You have attempted to invoke a functor with the incorrect number of Arguments.";
}

#define DECLARE_ARGLIST                        \
                                               \
typedef ARGS                        args;      \
typedef ArgList<args>               arglist;   \
typedef typename arglist::arg1_type arg1_type; \
typedef typename arglist::arg2_type arg2_type; \
typedef typename arglist::arg3_type arg3_type; \
typedef typename arglist::arg4_type arg4_type; \
typedef typename arglist::arg5_type arg5_type; \
typedef typename arglist::arg6_type arg6_type; \
typedef typename arglist::arg7_type arg7_type; \
typedef typename arglist::arg8_type arg8_type; \
typedef typename arglist::arg9_type arg9_type; \
                                               \
static const int arity = arglist::arity;

  

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
