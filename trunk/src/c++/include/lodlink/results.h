#ifndef LODLINK_RESULTS_H
#define LODLINK_RESULTS_H
//============================================================================
// File:      results.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/18/2 - created.                                   djb
//                                                                          
// Notes:     base classes for program task results.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <ostream>
#include <string>
#include <vector>
#include <utility>
#include <limits>
#include "numerics/matrix.h"
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "rped/rped.h"
#include "lodlink/mle_sub_model.h"
#include "lodlink/output.h"
#include "lodlink/definitions.h"

using std::ostream;
using std::string;
using std::vector;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

bool  converged(const MAXFUN::Results& data);
double  max_checked_value(const MAXFUN::Results& data);
bool  one_non_zero_std_error(const MAXFUN::Results& data);
double  lod_score(double alt_ln_like, double null_ln_like);

//----------------------------------------------------------------------------
//  Class:    task_result
//                                                                          
//  Purpose:  base class for results of tasks.
//                                                                          
//----------------------------------------------------------------------------
//
struct task_result
{
  task_result();
  virtual ~task_result();
  
  virtual void  write_summary(ostream& out) const = 0;
  virtual void  write_detail(ostream& out) const = 0;
  virtual void  write_vc_matrix(ostream& out) const = 0;
};

//----------------------------------------------------------------------------
//  Class:    hypothesis_result
//                                                                         
//  Purpose:  base class for linkage and homogeneity test results.
//                                                                          
//----------------------------------------------------------------------------
//
struct hypothesis_result : public task_result
{
  hypothesis_result();
  virtual ~hypothesis_result() = 0;

  double  chi_sq_stat() const;
  double  lod_score() const;
  void  dump(ostream& out) const;
  void  set_var_cov(const MAXFUN::Results& data);
  void  write_var_cov(ostream& out) const;
    static void  write_vc_row(ostream& out, const Matrix2D<double>& matrix, size_t row);
  
  double  null_ln_like;      // Natural log of likelihood.
  double  alt_ln_like;       // Natural log of likelihood.
  
  Matrix2D<double>  var_cov;
  
  string  trait;
  string  marker;
};

/* 
    Posterior data structures -
   
    test    per analysis
      |
      vector<result>    per locus
              |
            locus
            vector<pedigree_posterior>    per pedigree
                           |
                      pedigree name
                      posterior
                      vector<subpedigree_posterior>   per subpedigree
                                      |
                                member name
                                posterior
                                                  
*/

//----------------------------------------------------------------------------
//  Class:    subpedigree_posterior
//                                                                          
//  Purpose:  container for posterior likelihood of a subpedigree calculated 
//            using estimates for theta(s) and alpha obtained under Smith's 
//            heterogeneity model.
//                                                                          
//----------------------------------------------------------------------------
//
struct subpedigree_posterior
{
  subpedigree_posterior(const string& mn = "", double l = QNAN);
  
  bool  operator==(const subpedigree_posterior& right) const;

  string  member_name;
  double  posterior;
};


//----------------------------------------------------------------------------
//  Class:    pedigree_posterior
//                                                                          
//  Purpose:  container for posterior likelihoods of a pedigree and its sub-
//            pedigrees calculated using estimates for theta(s) and alpha
//            obtained under Smith's heterogeneity model.
//                                                                          
//----------------------------------------------------------------------------
//
struct pedigree_posterior
{
  pedigree_posterior(const string& pn = "", double l = QNAN);
  
  bool  operator==(const pedigree_posterior& right) const;

  string  pedigree_name;
  double  posterior;
  vector<subpedigree_posterior>  sub_posteriors;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_alt_result
//                                                                         
//  Purpose:  store data common to Smith's test for homogeneity and Faraway's
//            test for linkage under Smith's model so that this information
//            will not be calculated twice.
//
//  Note:     NOT PART OF THE RESULT HIERARCHY. 
//            
//----------------------------------------------------------------------------
//
struct non_ss_alt_result 
{
  non_ss_alt_result();

  double  alt_ln_like;
  double  alt_theta;
  double  alt_theta_ub;
  double  alpha;                // Proportion of families w. linkage.
  Matrix2D<double>  var_cov;      
  string  marker;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_smiths_faraways_result
//                                                                         
//  Purpose:  base class for result of Smith's test for homogeneity and 
//            Faraway's test for linkage assuming heterogeneity w. a single
//            recombinations fraction.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_smiths_faraways_result : hypothesis_result
{
  non_ss_smiths_faraways_result();
  
  virtual double  p_value() const = 0;
  
  double  alt_theta;
  double  alt_theta_ub;
  double  null_theta;
  double  null_theta_ub;
  double  alpha;         // Proportion of families w. linkage.  
  
  vector<pedigree_posterior>  posteriors;
};


//----------------------------------------------------------------------------
//  Class:    ss_alt_result
//                                                                         
//  Purpose:  store data common to Smith's test for homogeneity and Faraway's
//            test for linkage under Smith's model so that this information
//            will not be calculated twice.
//
//  Note:     NOT PART OF THE RESULT HIERARCHY. 
//            
//----------------------------------------------------------------------------
//
struct ss_alt_result 
{
  ss_alt_result();

  double      alt_ln_like;  
  theta_pair  alt_thetas;
  theta_pair  alt_theta_ubs;
  double      alpha;           // Proportion of families w. linkage.  
  Matrix2D<double>  var_cov;  
  string      marker;
};


//----------------------------------------------------------------------------
//  Class:    ss_smiths_faraways_result
//                                                                         
//  Purpose:  base class for result of Smith's test for homogeneity and 
//            Faraway's test for linkage assuming heterogeneity w. two
//            recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_smiths_faraways_result : hypothesis_result
{
  ss_smiths_faraways_result();
  
  virtual double  p_value() const = 0;
  
  theta_pair  alt_thetas;
  theta_pair  alt_theta_ubs;
  theta_pair  null_thetas;
  theta_pair  null_theta_ubs;
  double  alpha;               // Proportion of families w. linkage.  
  
  vector<pedigree_posterior>  posteriors;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_faraways_result
//                                                                          
//  Purpose:  represents results of Faraway's test for linkage under Smith's
//            heterogeneity model w average recombination fraction.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_faraways_result : public non_ss_smiths_faraways_result
{
  double  p_value() const;
  double  get_subpedigree_posterior(const string& pedigree, const string& member) const;

  static void  build_headers();

  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_family_detail(ostream& out, const string& pedigree, const string& member) const;
  void  write_vc_matrix(ostream& out) const;
  
  static header  meta;
  static header  columns;  
  static header  detail_columns;
  static header  vc_meta;
};


//----------------------------------------------------------------------------
//  Class:    ss_faraways_result
//                                                                          
//  Purpose:  represents results of Faraway's test for linkage under Smith's
//            heterogeneity model with sex specific recombination fractions.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_faraways_result : public ss_smiths_faraways_result
{
  double  p_value() const;
  double  get_subpedigree_posterior(const string& pedigree, const string& member) const;  
  
  static void  build_headers();

  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_family_detail(ostream& out, const string& pedigree, const string& member) const;  
  void  write_vc_matrix(ostream& out) const;
  
  static header  meta;
  static header  columns;
  static header detail_columns;
  static header  vc_meta;  
};


//----------------------------------------------------------------------------
//  Class:    non_ss_smiths_result
//                                                                         
//  Purpose:  represents results of non-sex-specific Smith's test for linkage
//            homogeneity.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_smiths_result : public non_ss_smiths_faraways_result
{
  double  p_value() const;

  static void  build_headers();
  
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_vc_matrix(ostream& out) const;
  
  static header  meta;
  static header  columns;
  static header  vc_meta;  
};

//----------------------------------------------------------------------------
//  Class:    ss_smiths_result
//                                                                         
//  Purpose:  represents results of sex-specific Smith's test for linkage
//            homogeneity.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_smiths_result : public ss_smiths_faraways_result
{
  double  p_value() const;
  
  static void  build_headers();
  
  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_vc_matrix(ostream& out) const;
  
  static header  meta;
  static header  columns;
  static header  vc_meta;  
};

#include "lodlink/results.ipp"
}
}

#endif
