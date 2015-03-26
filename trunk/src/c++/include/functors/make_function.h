#ifndef FUNCTORS_MAKE_FUNCTION_H
#define FUNCTORS_MAKE_FUNCTION_H

#include "boost/mpl/vector.hpp"
#include "boost/mpl/at.hpp"
#include "boost/function.hpp"
#include "functors/at_idx.h"

namespace SAGE     {
namespace FUNCTORS {

template<int ARGCOUNT, typename R, typename ARGS> struct make_function_impl { };

template<typename R, 
         typename ARGS> struct make_function_impl<0, R, ARGS> { typedef typename boost::function0<R>                                 type; };
template<typename R, 
         typename ARGS> struct make_function_impl<1, R, ARGS> { typedef typename boost::function1<R, typename at_idx<ARGS, 0>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<2, R, ARGS> { typedef typename boost::function2<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<3, R, ARGS> { typedef typename boost::function3<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<4, R, ARGS> { typedef typename boost::function4<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<5, R, ARGS> { typedef typename boost::function5<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type,
                                                                                                     typename at_idx<ARGS, 4>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<6, R, ARGS> { typedef typename boost::function6<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type,
                                                                                                     typename at_idx<ARGS, 4>::type,
                                                                                                     typename at_idx<ARGS, 5>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<7, R, ARGS> { typedef typename boost::function7<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type,
                                                                                                     typename at_idx<ARGS, 4>::type,
                                                                                                     typename at_idx<ARGS, 5>::type,
                                                                                                     typename at_idx<ARGS, 6>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<8, R, ARGS> { typedef typename boost::function8<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type,
                                                                                                     typename at_idx<ARGS, 4>::type,
                                                                                                     typename at_idx<ARGS, 5>::type,
                                                                                                     typename at_idx<ARGS, 6>::type,
                                                                                                     typename at_idx<ARGS, 7>::type> type; };
template<typename R, 
         typename ARGS> struct make_function_impl<9, R, ARGS> { typedef typename boost::function9<R, typename at_idx<ARGS, 0>::type,
                                                                                                     typename at_idx<ARGS, 1>::type,
                                                                                                     typename at_idx<ARGS, 2>::type,
                                                                                                     typename at_idx<ARGS, 3>::type,
                                                                                                     typename at_idx<ARGS, 4>::type,
                                                                                                     typename at_idx<ARGS, 5>::type,
                                                                                                     typename at_idx<ARGS, 6>::type,
                                                                                                     typename at_idx<ARGS, 7>::type,
                                                                                                     typename at_idx<ARGS, 8>::type> type; };

///
/// Construct a boost::function with the given result_type and argument list.
/// \param R The result type
/// \param ARGS A boost::mpl::vector storing the the argument list
template<typename R, typename ARGS = boost::mpl::vector<> > struct make_function
{
  typedef typename make_function_impl<boost::mpl::size<ARGS>::value, R, ARGS>::type type;
};

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
