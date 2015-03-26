#ifndef FCORDEFS_H
#define FCORDEFS_H

//****************************************************************************
//* File:      definitions.h                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                    1.0  var-cov structure added               yjs Jan 01 *
//*                    2.0  structure for maintype added          yjs May 01 *
//*                                                                          *
//* Notes:     This header file contains various type definition statements  *
//*            used through out the fcor files.                              *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "pairs/relmatrix.h"
#include "pairs/reltypename.h"
#include "containers/ldatamap.h"
#include "numerics/xlimits.h"
#include "numerics/corinfo.h"
#include "numerics/functions.h"
#include "numerics/print_util.h"
#include "numerics/matfunctions.h"

namespace SAGE {
namespace FCOR {

enum pair_type    {SUBTYPES = 0, MAINTYPES, BOTH, SELF, UNRELATED, MATE, INTRA, OTHER};
enum weight_type  {PAIR_WISE = 0, UNIFORM, MEAN, WEIGHT_COUNT};
enum output_type  {SUMMARY = 0, DETAILED, PAIR_COUNT, XLS_FORMAT};

typedef KahanAdder<double>     internal_real_type;
typedef pair<string, size_t>   name_index_type;
typedef pair<string, string>   name_name_type;

//
// definitions for variance_covariance matrix.
//
struct var_cov_param
{
  enum matrix_type {SINGLE = 0, JOINT};
  
  var_cov_param();
  var_cov_param(matrix_type m, name_name_type c, const vector<name_index_type>& t);
  
  string name()             const;
  string trait_name()       const;
  string correlation_name() const;
  string matrix_type_name() const;

  matrix_type              m_type;
  name_name_type           correlations;
  vector<name_index_type>  traits;
};

//
// options for analysis
//
struct analysis_option_type
{
  analysis_option_type();

  pair_type             pairset;

  weight_type           class_weight;

  int                   generation_limit;

  bool                  standard_error;
  bool                  conservative;
  bool                  gender_name;
  bool                  homogeneity_test;
  bool                  individual_homog;
  bool                  valid_only_homog;
  bool                  detailed_output;
  bool                  pair_output;
  bool                  xls_output;
  bool                  KE_with_optimal;

  vector<var_cov_param> var_covs;
};

//
// weighted memeber pair.
//
struct WeightedMemberPair
{
  WeightedMemberPair(size_t w_ct);
  WeightedMemberPair(const_pedigree_member_pair, size_t w_ct);
  
  const_pedigree_member_pair  member_pair;
  vector<double>              pair_weight;


  bool operator()(const WeightedMemberPair& m1,
                  const WeightedMemberPair& m2 ) const;
                  
  bool operator< (const WeightedMemberPair& m)   const;
  bool operator==(const WeightedMemberPair& m)   const;
};

//
// structures for analysis result.
//
struct correlation_result
{
  correlation_result();
  correlation_result(size_t w_count);

  size_t          count;
  vector<double>  correlation;
};

struct standard_error_result
{
  standard_error_result();
  standard_error_result(size_t w_count);

  size_t          least_pair_count;
  bool            replaced;
  double          w1;
  vector<double>  standard_error;
};

struct pooled_correlation_result
{
  pooled_correlation_result();

  double  correlation;
  double  variance;
  double  chi_square;
};

struct pairset_result
{
  pairset_result();
  pairset_result(size_t m, size_t w, bool se);

  Matrix2D< correlation_result >         corr;
  Matrix2D< standard_error_result >      std_err;
  Matrix2D< pooled_correlation_result >  pooled_cross_corr;

  double  pooled_cross_corr_chi_square;
  double  pooled_cross_corr_p_value;
};

struct Htest_result
{
  Htest_result();

  vector<double>  chi_square;
  vector<double>  p_value;

  bool            replaced;
};

struct var_cov_result
{
  var_cov_result();

  double  var_cov;
  bool    replaced;
  size_t  least_pair_count;
};

typedef Matrix2D<var_cov_result> var_cov_matrix;

struct pairset_info
{
  pair_type  type;
  string     name;
  string     gname;
  string     rcode;
  size_t     distance1;
  size_t     distance2;

  size_t     total_pair_count;
  size_t     distinctive_pair_count;
};

//
// type definitions.
//
typedef vector< WeightedMemberPair >              pairset_type;
typedef vector< pairset_type >                    pairset_by_pedigree_type;

typedef pairset_type::iterator                    pairset_iterator;
typedef pairset_type::const_iterator              pairset_const_iterator;

typedef pairset_by_pedigree_type::iterator        pairset_by_pedigree_iterator;
typedef pairset_by_pedigree_type::const_iterator  pairset_by_pedigree_const_iterator;

typedef vector< pairset_by_pedigree_type >        pairset_vector;
typedef vector< pairset_info >                    pairset_info_vector;
typedef vector< pairset_result >                  pairset_result_vector;

typedef vector< Htest_result >                    Htest_result_vector;
typedef vector< var_cov_matrix >                  var_cov_result_vector;

typedef vector< CorrelationInfo >                 corinfo_by_weight_type;
typedef vector< corinfo_by_weight_type >          corinfo_vector;
typedef vector< TriangleMatrix< size_t > >        pair_to_corinfo_map;
//
// definitions for optimal_weight.
//
typedef pair< Matrix2D<double>, Matrix2D<double> > weight_matrix;
typedef vector< weight_matrix >                   weight_matrix_by_pedigree;
typedef vector< weight_matrix_by_pedigree >       weight_matrix_vector;

//
//
typedef vector< CompleteRelationType >            subtype_vector;

typedef map< CompleteRelationType,
             pairset_by_pedigree_type,
             CompleteRelationMateLess<> >         sub_pairset_map;

typedef sub_pairset_map::const_iterator           sub_pairset_const_iterator;

typedef map< CompleteRelationType,
             size_t,
             CompleteRelationMateLess<> >         sub_count_map;

typedef LinearDatamap< RelationType,
                       size_t,
                       RelationMateLess >         main_to_sub_map;

typedef main_to_sub_map::iterator                 main_to_sub_iterator;
typedef main_to_sub_map::const_iterator           main_to_sub_const_iterator;
typedef main_to_sub_map::type_iterator            main_to_sub_type_iterator; 
typedef main_to_sub_map::type_const_iterator      main_to_sub_type_const_iterator;

// ---------------------------------------------------------------------------
// Inline implementation
// ---------------------------------------------------------------------------

inline
var_cov_param::var_cov_param()
{
  m_type = SINGLE;
  traits.resize(0);
}

inline
var_cov_param::var_cov_param(matrix_type m, name_name_type c, const vector<name_index_type>& t)
{
  m_type       = m;
  correlations = c;
  traits       = t;
}

inline string
var_cov_param::name() const
{
  string vc = "Variance-Covariance matrix for Correlations of\n";

  vc.append(correlation_name());
  vc.append(trait_name());
  vc.append(matrix_type_name());
  
  return vc;
}

inline string
var_cov_param::trait_name() const
{
  string vc = "     trait(s) : ";
  for( size_t t = 0; t < traits.size(); ++t )
  {
    vc.append(traits[t].first);
    vc.append(" ");
  }
  
  return vc;
}

inline string
var_cov_param::correlation_name() const
{
  string vc = "     ";
  vc.append(correlations.first);
  vc.append("\n       with\n");
  vc.append("     ");
  vc.append(correlations.second);
  vc.append("\n\n");

  return vc;
}

inline string
var_cov_param::matrix_type_name() const
{
  string vc;

  if( m_type == SINGLE )
    vc = "- SINGLY";
  else  
    vc = "- JOINTLY";

  vc.append("\n");
  
  return vc;
}

// ---------------------------------------------------------------------------

inline
analysis_option_type::analysis_option_type()
{
  pairset           = SUBTYPES;
  class_weight      = WEIGHT_COUNT;

  generation_limit  = 2;

  standard_error    = true;
  conservative      = false;
  gender_name       = true;
  homogeneity_test  = false;
  individual_homog  = false;
  valid_only_homog  = false;
  detailed_output   = false;
  pair_output       = false;
  xls_output        = false;
  KE_with_optimal   = false;

  var_covs.resize(0);
}

// ---------------------------------------------------------------------------

inline
WeightedMemberPair::WeightedMemberPair(size_t w_ct)
{ 
  pair_weight.resize(w_ct, 1.0);
}

inline
WeightedMemberPair::WeightedMemberPair(const_pedigree_member_pair mp, size_t w_ct)
{ 
  member_pair = mp;
  pair_weight.resize(w_ct, 1.0);
}

inline bool
WeightedMemberPair::operator()(const WeightedMemberPair& m1,
                               const WeightedMemberPair& m2 ) const
{
  if( m1.member_pair.first->pedigree() == m2.member_pair.first->pedigree() )
    if( m1.member_pair.first == m2.member_pair.first )
      return m1.member_pair.second < m2.member_pair.second;
    else
      return m1.member_pair.first  < m2.member_pair.first;

  return m1.member_pair.first->pedigree() < m2.member_pair.first->pedigree();
}                                   

inline bool
WeightedMemberPair::operator<(const WeightedMemberPair& m) const
{
  for( size_t i = 0; i < pair_weight.size(); ++i )
    if( pair_weight[i] >= m.pair_weight[i] )
      return false;
    
  return true;
}

inline bool
WeightedMemberPair::operator==(const WeightedMemberPair& m) const
{
  if(    member_pair.first  != m.member_pair.first
      || member_pair.second != m.member_pair.second )
    return false;

  for( size_t i = 0; i < pair_weight.size(); ++i )
    if( pair_weight[i] != m.pair_weight[i] )
      return false;

  return true;
}

// ---------------------------------------------------------------------------

inline
correlation_result::correlation_result()
{
  count = 0;
  correlation.resize(0);
}

inline
correlation_result::correlation_result(size_t w_ct)
{
  count = 0;
  correlation.resize(w_ct, std::numeric_limits<double>::quiet_NaN());
}

// ---------------------------------------------------------------------------

inline
standard_error_result::standard_error_result()
{
  least_pair_count = 0;
  replaced         = false;
  w1               = 0.;
  standard_error.resize(0);
}

inline
standard_error_result::standard_error_result(size_t w_ct)
{
  least_pair_count = 0;
  replaced         = false;
  w1               = 0.;
  standard_error.resize(w_ct, std::numeric_limits<double>::quiet_NaN());
}

// ---------------------------------------------------------------------------

inline
pooled_correlation_result::pooled_correlation_result()
{
  correlation = std::numeric_limits<double>::quiet_NaN();
  variance    = std::numeric_limits<double>::quiet_NaN();
  chi_square  = std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------

inline
pairset_result::pairset_result()
{
  corr.resize(0, 0);
  std_err.resize(0, 0);
  pooled_cross_corr.resize(0, 0);

  pooled_cross_corr_chi_square = std::numeric_limits<double>::quiet_NaN();
  pooled_cross_corr_p_value    = std::numeric_limits<double>::quiet_NaN();
}

inline
pairset_result::pairset_result(size_t m, size_t w, bool se)
{
  corr.resize(m, m, correlation_result(w));

  if( se )
  {
    std_err.resize(m, m, standard_error_result(w));
    pooled_cross_corr.resize(m, m, pooled_correlation_result());
  }

  pooled_cross_corr_chi_square = std::numeric_limits<double>::quiet_NaN();
  pooled_cross_corr_p_value    = std::numeric_limits<double>::quiet_NaN();
}

// ---------------------------------------------------------------------------

inline
Htest_result::Htest_result()
{
  chi_square.resize(0);
  p_value.resize(0);

  replaced = false;
}

// ---------------------------------------------------------------------------

inline
var_cov_result::var_cov_result()
{
  var_cov = std::numeric_limits<double>::quiet_NaN();
  replaced = false;
  least_pair_count = 0;
}

// ---------------------------------------------------------------------------

inline bool is_selfclass(const CompleteRelationType& sub_rel_type)
{
  if(    sub_rel_type.relationship.distance1() + sub_rel_type.relationship.distance2() == 0
      && sub_rel_type.relationship.bridge() == SELF_BRIDGE )
    return true;

  return false;
}

inline bool is_selfclass(const RelationType& main_rel_type)
{
  if(    main_rel_type.distance1() + main_rel_type.distance2() == 0
      && main_rel_type.bridge() == SELF_BRIDGE )
    return true;

  return false;
}

inline bool is_mate(const CompleteRelationType& sub_rel_type)
{
  if(    sub_rel_type.relationship.distance1() + sub_rel_type.relationship.distance2() == 0
      && sub_rel_type.relationship.bridge() == NO_BRIDGE
      && sub_rel_type.relationship.has_offspring() )
    return true;

  return false;
}

inline bool is_mate(const RelationType& main_rel_type)
{
  if(    main_rel_type.distance1() + main_rel_type.distance2() == 0
      && main_rel_type.bridge() == NO_BRIDGE
      && main_rel_type.has_offspring() )
    return true;

  return false;
}

inline bool is_unrelated(const CompleteRelationType& sub_rel_type)
{
  if(     sub_rel_type.relationship.distance1() + sub_rel_type.relationship.distance2() == 0
      &&  sub_rel_type.relationship.bridge() == NO_BRIDGE
      && !sub_rel_type.relationship.has_offspring() )
    return true;

  return false;
}

inline bool is_unrelated(const RelationType& main_rel_type)
{
  if(     main_rel_type.distance1() + main_rel_type.distance2() == 0
      &&  main_rel_type.bridge() == NO_BRIDGE
      && !main_rel_type.has_offspring() )
    return true;

  return false;
}

inline bool is_intraclass(const CompleteRelationType& sub_rel_type)
{
  if(    is_unrelated(sub_rel_type)
      || is_mate(sub_rel_type)
      || is_selfclass(sub_rel_type) )
    return false;

  if( sub_rel_type.relationship.distance1() == sub_rel_type.relationship.distance2() )
  {
    std::pair<string, string> gs = SAGE::gender_strings(sub_rel_type);
    
    if( gs.first == gs.second )
      return true;
  }
  
  return false;
}

inline bool is_intraclass(const RelationType& main_rel_type)
{
  if(    is_unrelated(main_rel_type)
      || is_mate(main_rel_type)
      || is_selfclass(main_rel_type)  )
    return false;

  if( main_rel_type.distance1() == main_rel_type.distance2() )
      return true;
  
  return false;
}

inline bool is_invalid_pair_type(const pairset_info& pinfo, size_t gen_limit, size_t t_count)
{
  if( pinfo.total_pair_count < 3 )
    return true;

  if(    pinfo.distance1 > gen_limit
      || pinfo.distance2 > gen_limit )
    return true;

  if( pinfo.type == UNRELATED )
    return true;

  if( t_count < 2 && pinfo.type == SELF )
    return true;

  return false;
}

inline void print_vec_matrix(const Matrix2D< vector<size_t> >& m, ostream &out, string title = "")
{
  out << title << " : "
      << m.rows() << " x " << m.cols() << " matrix" << endl;

  for( size_t i = 0; i < m.rows(); ++i )
  {
    for( size_t j = 0; j < m.cols(); ++j )
    {
      out << "(";
      for( size_t w = 0; w < m(i, j).size(); ++w )
      {
        if( m(i, j)[w] != (size_t)-1 )
          out << setw(10) << m(i, j)[w] << " ";
        else
          out << "********** ";
      }
      out << ") ";
    }
    out << endl;
  }
  out << endl;
}

inline void print_vec_matrix(const Matrix2D< vector<double> >& m, ostream &out, string title = "")
{
  out << title << " : "
      << m.rows() << " x " << m.cols() << " matrix" << endl;

  for( size_t i = 0; i < m.rows(); ++i )
  {
    for( size_t j = 0; j < m.cols(); ++j )
    {
      out << "(";
      for( size_t w = 0; w < m(i, j).size(); ++w )
      {
        if( !SAGE::isnan(m(i, j)[w]) )
          out << setw(10) << m(i, j)[w] << " ";
        else
          out << "********** ";
      }
      out << ") ";
    }
    out << endl;
  }
  out << endl;
}

inline void print_matrix(const Matrix2D<double>& m, ostream &out, string title = "")
{
  out << title << " : " << endl;

  for( size_t i = 0; i < m.rows(); ++i )
  {
    for( size_t j = 0; j < m.cols(); ++j )
    {
      if( !SAGE::isnan(m(i, j)) )
        out << setw(10) << m(i, j) << " ";
      else
        out << "********** ";
    }
    out << endl;
  }
  out << endl;
}

inline void print_matrix(const Matrix2D<internal_real_type>& m, ostream &out, string title = "")
{
  out << title << " : " << endl;

  for( size_t i = 0; i < m.rows(); ++i )
  {
    for( size_t j = 0; j < m.cols(); ++j )
    {
      if( !SAGE::isnan(m(i, j)) )
        out << setw(10) << m(i, j) << " ";
      else
        out << "********** ";
    }
    out << endl;
  }
  out << endl;
}

} // end of namespace FCOR
} // end of namespace SAGE

#endif

