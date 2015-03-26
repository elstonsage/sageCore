#ifndef LODLINK_LINKAGE_RESULTS_H
#define LODLINK_LINKAGE_RESULTS_H
//============================================================================
// File:      linkage_results.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/25/2 created         -djb
//                                                                          
// Notes:     Defines classes for storing and writing results of linkage tests.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "maxfun/sub_model.h"
#include "maxfun/maxfun.h"
#include "lodlink/instructions.h"
#include "lodlink/results.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"
#include "lodlink/output.h"

using std::vector;
using std::string;
using std::ostream;
using std::endl;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    lod_ratio_result
//                                                                          
//  Purpose:  represents a lod score likelihood ratio test result.
//                                                                          
//----------------------------------------------------------------------------
//
struct lod_ratio_result : public hypothesis_result
{
  lod_ratio_result();
  virtual ~lod_ratio_result() = 0;
  
  double  p_value_u_bound() const;
  double  chi_sq_stat() const;
  void  set_var_cov_relaxed(const MAXFUN::Results& data);
  void  write_var_cov(ostream& out) const;
  
  // - Linkage ratio test done only if certain conditions are met.
  //
  bool  perform_test;             
  Matrix2D<double>  var_cov_relaxed;  
};


//----------------------------------------------------------------------------
//  Class:    non_ss_lod_ratio_result
//                                                                          
//  Purpose:  represents a lod score likelihood ratio test result.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_lod_ratio_result : public lod_ratio_result
{
  non_ss_lod_ratio_result();

  double  p_value() const;
  
  static void  build_headers();
  
  double  restricted_alt_theta;
  double  restricted_alt_theta_ub;
  double  relaxed_alt_theta;
  double  relaxed_alt_theta_ub;
  
  void  set_perform_test(const MAXFUN::Results& data);
  static bool  bound_at_point_five(const MAXFUN::Results& data);  
  
  void  dump(ostream& out) const;
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_vc_matrix(ostream& out) const;
  
  static header meta;
  static header columns;
  static header vc_meta;
};


//----------------------------------------------------------------------------
//  Class:    ss_lod_ratio_result
//                                                                          
//  Purpose:  represents a lod score likelihood ratio test result.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_lod_ratio_result : public lod_ratio_result
{
  ss_lod_ratio_result();

  double  p_value() const;
  
  static void  build_headers();
  
  theta_pair  restricted_alt_thetas;
  theta_pair  restricted_alt_theta_ubs;
  theta_pair  relaxed_alt_thetas;
  theta_pair  relaxed_alt_theta_ubs;
  
  void  set_perform_test(const MAXFUN::Results& data);
  static bool  either_bound_at_point_five(const MAXFUN::Results& data);    
  
  void  dump(ostream& out) const;
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_vc_matrix(ostream& out) const;  
  
  static header meta;
  static header columns;
  static header vc_meta;
};


//----------------------------------------------------------------------------
//  Class:    cleves_elston_result
//                                                                          
//  Purpose:  represents Cleves-Elston likelihood ratio test result.
//                                                                          
//----------------------------------------------------------------------------
//
struct cleves_elston_result : public hypothesis_result
{
  cleves_elston_result();

  double  p_value() const;

  static void build_headers();
  
  theta_pair  alt_thetas;
  theta_pair  alt_theta_ubs;
  theta_pair  null_thetas;
  theta_pair  null_theta_ubs;
  
  void  dump(ostream& out) const;
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_vc_matrix(ostream& out) const;  
  
  static header meta;
  static header columns;
  static header vc_meta;
};

#include "lodlink/linkage_results.ipp"

}
}

#endif
