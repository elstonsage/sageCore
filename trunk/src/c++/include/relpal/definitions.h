#ifndef RELPAL_DEFS_H
#define RELPAL_DEFS_H

#include "numerics/fmatrix.h"
#include "numerics/trimatrix.h"
#include "numerics/sinfo.h"
#include "numerics/mt.h"
#include "numerics/print_util.h"
#include "palbase/pair_filter.h"
#include "palbase/pal_ibd.h"
#include "palbase/pair_info_file.h"


namespace SAGE   {
namespace RELPAL {

#define DEFAULT_PARAM_NAME_MAX  25

typedef FPED::Subpedigree               subpedigree;
typedef FPED::Member                    member;

typedef PALBASE::rel_pair_data          rel_pair_data;
typedef PALBASE::relative_pairs         relative_pairs;
typedef PALBASE::rel_pair               rel_pair;
typedef PALBASE::pair_filter            pair_filter;

typedef FortranMatrix<double>           matrix;
typedef TriangleMatrix<double>          trimatrix;

// Pair type option supported
//
enum use_pair_type { FSIB=0, HSIB, SIB, ALL };

// Member data option
//
enum use_member_type { EVERY=0, INFORMATIVE_LOCAL, INFORMATIVE_REGION, DEFAULT };

// Type of regression analysis supported
//
// STSM = single trait single marker regression
// STMM = single trait multiple marker regression
// STZM = single trait zero marker regression
// MTSM = multiple trait single marker regression
// MTMM = multiple trait multiple marker regression
// MTZM = multiple trait zero marker regression
//
enum regression_type { STSM=0, STMM, STZM, MTSM, MTMM, MTZM };

//
// Input Data
//
struct member_less : public std::binary_function<mem_pointer, mem_pointer, bool>
{
  member_less() {}
  bool operator()(mem_pointer i1, mem_pointer i2) const
  {
    //return i1->name() < i2->name();
    return i1->index() < i2->index();
  }
};

typedef set< mem_pointer, member_less >  member_set;
typedef member_set::iterator             member_set_iterator;
typedef member_set::const_iterator       member_set_const_iterator;

struct subped_map_info
{
  member_set               members;
  map< id_pair, size_t >   member_to_pair;
};

typedef map< FPED::SubpedigreeConstPointer, subped_map_info >      pair_map_type;
typedef pair_map_type::iterator                                    pair_map_iterator;
typedef pair_map_type::const_iterator                              pair_map_const_iterator;

struct subpedigree_data
{
  vector< mem_pointer >     members;
  TriangleMatrix< size_t >  member_to_pair;
};

typedef vector< subpedigree_data >     analysis_data_type;

//
// Result
//
struct regression_result
{
  matrix   ind_beta;
  matrix   ind_variance;

  matrix   ped_beta;
  matrix   ped_variance;
};

struct score_test_result
{
  matrix          U;
  vector<matrix>  sigma;

  vector<double>  T;
  vector<double>  correction;
  vector<double>  emp_pvalue;

  vector<size_t>  rep_count;
};

struct analysis_result
{
  regression_result H0_result;
  regression_result H1_result;

  score_test_result score_result;
};

//
//
//
struct resid_info
{
  size_t sp_id;
  size_t ind_id;
  double res_val;
  double rank_val;
  double inr_val;

  resid_info(size_t sp, size_t i, double r)
   : sp_id(sp), ind_id(i), res_val(r), rank_val(QNAN), inr_val(QNAN) {}
};

struct resid_info_less : public std::binary_function<resid_info, resid_info, bool>
{
  resid_info_less() {}

  bool operator()(resid_info r1, resid_info r2) const
  {
    return r1.res_val < r2.res_val;
  }
};

//
// Inline functions
//
inline string
get_regression_type_to_string(regression_type reg_ty)
{
  switch(reg_ty)
  {
    case STSM: return "single trait, single marker";  break;
    
    case STMM: return "single trait, multiple markers";  break;

    case STZM: return "single trait, zero marker";  break;

    case MTSM: return "multiple traits, single marker";  break;

    case MTMM: return "multiple traits, multiple markers";  break;

    case MTZM: return "multiple traits, zero marker";  break;

    default:   return "invalid"; break;
  }

  return "invalid";
}

inline void
print_trimatrix(const trimatrix& w, ostream& o, string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';
  
  for( size_t i = 0; i < w.size(); ++i )
  {
    if( i )
      o << std::endl << std::setw(offset) << "";

    for( size_t j = 0; j <= i; ++j )
    {
      if( fabs(w(i,j)) < 10e-15 )
        o << std::setw(8) << 0.0 << " ";
      else
        o << std::setw(8) << w(i,j) << " ";
    }
  }

  o << " ]" << std::endl;

  return;
}

inline void
print_trimatrix_first10(const trimatrix& w, ostream& o, string title = "")
{
  int offset = 1;

  if(title.size())
  {
    o << title << " = [";
    offset += title.size() + 3;
  }
  else
    o << '[';

  size_t size = 10;
  if( w.size() < 10 )
    size = w.size();

  for( size_t i = 0; i < size; ++i )
  {
    if( i )
      o << std::endl << std::setw(offset) << "";

    for( size_t j = 0; j <= i; ++j )
    {
      if( fabs(w(i,j)) < 10e-15 )
        o << std::setw(8) << 0.0 << " ";
      else
        o << std::setw(8) << w(i,j) << " ";
    }
  }

  o << " ]" << std::endl;

  return;
}


} // end of namespace RELPAL
} // end of namespace SAGE

#endif

