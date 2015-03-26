#ifndef LODLINK_HOMOGENEITY_RESULTS_H
#define LODLINK_HOMOGENEITY_RESULTS_H
//============================================================================
// File:      homogeneity_results.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/13/3 created         -djb
//                                                                          
// Notes:     Defines classes for storing and writing results of homogeneity tests.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <utility>
#include <iostream>
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "maxfun/sub_model.h"
#include "lodlink/instructions.h"
#include "lodlink/results.h"
#include "lodlink/linkage_results.h"
#include "lodlink/likelihood.h"
#include "lodlink/max_opt.h"

using std::vector;
using std::string;
using std::ostream;
using std::endl;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Class:    non_ss_mortons_result
//                                                                         
//  Purpose:  represent result of Morton's test for linkage homogeneity using
//            an average recombination fraction.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_mortons_result : hypothesis_result
{
  non_ss_mortons_result();
  
  double  p_value() const;
  
  static void  build_headers();
  
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_group_detail(ostream& out, const string& group) const;
  void  write_vc_matrix(ostream& out) const;
  
  // - <group name, pair<theta, ln maximum likelihood> >
  //
  std::map<string, pair<double, double> >  group_results;
  double  group_theta_ub;
  double  null_theta;
  double  null_theta_ub;
  
  static header  meta;
  static header  columns;
  static header  detail_columns;
};

//----------------------------------------------------------------------------
//  Class:    ss_mortons_result
//                                                                         
//  Purpose:  represent result of Morton's test for linkage homogeneity using
//            separate male and female recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_mortons_result : hypothesis_result
{
  ss_mortons_result();
  
  double  p_value() const;
  
  static void  build_headers();
  
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_group_detail(ostream& out, const string& group) const;
  void  write_vc_matrix(ostream& out) const;
  
  // - <group name, pair<thetas, ln maximum likelihood> >
  //
  std::map<string, pair<theta_pair, double> >  group_results;
  theta_pair  group_theta_ubs;
  theta_pair  null_thetas;
  theta_pair  null_theta_ubs;
  
  static header meta;
  static header columns;
  static header detail_columns;
};

#include "lodlink/homogeneity_results.ipp"

}
}

#endif
