#ifndef FUNCTORS_AT_IDX_H
#define FUNCTORS_AT_IDX_H

#include "boost/mpl/at.hpp"
#include "boost/function.hpp"

namespace SAGE     {
namespace FUNCTORS {

/// Extracts the type at the indexed entry in the given sequence.
template<typename SEQ, int IDX> struct at_idx { typedef typename boost::mpl::at<SEQ, boost::mpl::int_<IDX> >::type type; };

} // End namespace FUNCTORS
} // End namespace SAGE

#endif
