#ifndef SIBPAL_DEFS_H
#define SIBPAL_DEFS_H

#include <utility>
#include <string>
#include <sstream>
#include "LSF/parse_ops.h"
#include "LSF/LSFsymbol.h"
#include "numerics/corinfo.h"
#include "numerics/histogram.h"
#include "numerics/functions.h"
#include "numerics/fmatrix.h"
#include "numerics/mt.h"
#include "numerics/print_util.h"
#include "containers/ldatamap.h"
#include "palbase/pair_filter.h"
#include "palbase/pal_ibd.h"
#include "palbase/pair_info_file.h"


namespace SAGE   {
namespace SIBPAL {

typedef PALBASE::rel_pair_data                  sib_pair_data;
typedef PALBASE::sibcluster_info                sibcluster_info;

typedef PALBASE::sibship_cluster                sibship_cluster;
typedef PALBASE::sibship_cluster_iterator       sibship_cluster_iterator;
typedef PALBASE::sibship_cluster_const_iterator sibship_cluster_const_iterator;
typedef PALBASE::sib_set                        sib_set;
typedef PALBASE::sib_mean_map                   sib_mean_map;

typedef PALBASE::relative_pairs                 relative_pairs;
typedef PALBASE::rel_pair                       sib_pair;
typedef PALBASE::pair_filter                    pair_filter;

typedef double                                  stable_double;
typedef FortranMatrix<double>                   matrix;
typedef FortranMatrix<stable_double>            stable_matrix;
typedef FortranMatrix<double>                   double_vector;

typedef map<pair<size_t, size_t>, matrix>       correlation_matrix_map;

typedef map<size_t, pair<double, double> >      w_map;

enum regression_type    { SINGLE_MARKER=0, MULTIPLE_MARKER, ZERO_MARKER };
enum weight_status_type { NORMALW = 0, INVERSEW, BESTW };
enum pair_type          { MM = 0, MF, FF, MIXED };
enum use_pair_type      { FSIB = 0, HSIB, SIB };

inline string
get_regression_type_to_string(regression_type reg_ty)
{
  switch(reg_ty)
  {
    case SINGLE_MARKER:   return "single marker regression"; break;

    case MULTIPLE_MARKER: return "multiple marker regression"; break;

    case ZERO_MARKER:     return "zero marker regression"; break;

    default:              return "invalid"; break;
  }

  return "invalid";
}

inline string
get_bool_to_string(bool opt)
{
  if( opt )
    return "yes";

  return "no";
}

struct pair_index
{
  pair_index()
   : sibship_size( (size_t)-1 ), sibship_number( (size_t)-1 ),
     sib_number( (size_t)-1 ) {}

  pair_index(const size_t s, const size_t n, const size_t i)
   : sibship_size(s), sibship_number(n), sib_number(i) {}

  size_t sibship_size;
  size_t sibship_number;
  size_t sib_number;
};

typedef vector<pair_index>                      simulation_map_type;

typedef vector<size_t>                          permutation_vector_type;
typedef vector<permutation_vector_type>         sibship_permutation_vector_type;
typedef vector<sibship_permutation_vector_type> group_permutation_vector_type;

struct MTRandomizer
{
  unsigned int operator()(size_t N) { return mt.uniform_integer(N); }
  MersenneTwister mt;
};

} // end of namespace SIBPAL
} // end of namespace SAGE

#endif

