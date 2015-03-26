#ifndef LODLINK_LODS_RESULTS_H
#define LODLINK_LODS_RESULTS_H
//============================================================================
// File:      lods_results.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/10/03 created         -djb
//                                                                          
// Notes:     Defines classes for storing and writing results of lod score
//            calculations.
//                                                                          
// Copyright (c) 2003 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <utility>
#include <algorithm>
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "maxfun/sub_model.h"
#include "lodlink/instructions.h"
#include "lodlink/results.h"
#include "lodlink/likelihood.h"

using std::vector;
using std::string;
using std::ostream;
using std::endl;
using std::pair;

namespace SAGE
{

namespace LODLINK
{

/* 
    Lod score data structures -
   
    test    per analysis
      |
      vector<result>    per locus
              |
            locus
            vector<theta, multipedigree_lod_score>     per theta
                                    |
                              lod score
                              vector<pedigree_lod_score>    per pedigree
                                             |
                                        pedigree name
                                        lod score
                                        vector<subpedigree_lod_score>   per subpedigree
                                                        |
                                                  member name
                                                  lod score
                                                  
*/

//----------------------------------------------------------------------------
//  Class:    subpedigree_lod_score
//                                                                          
//  Purpose:  represents a lod score of for a subpedigree.
//                                                                          
//----------------------------------------------------------------------------
//
struct subpedigree_lod_score
{
  subpedigree_lod_score(const string& mn = "", const double& lod = QNAN);

  bool  operator==(const subpedigree_lod_score& right) const;


  string      member_name;
  double  lod_score;
};


//----------------------------------------------------------------------------
//  Class:    pedigree_lod_score
//                                                                          
//  Purpose:  container for lod scores of a pedigree and its sub-
//            pedigrees.
//                                                                          
//----------------------------------------------------------------------------
//
struct pedigree_lod_score
{
  pedigree_lod_score(const string& pn = "", const double& lod = QNAN);
  
  bool  operator==(const pedigree_lod_score& right) const;

  void  write(ostream& out);
  
  string      pedigree_name;
  double  lod_score;
  vector<subpedigree_lod_score>  sub_lod_scores;
};

//----------------------------------------------------------------------------
//  Class:    multipedigree_lod_score
//                                                                          
//  Purpose:  container for lod scores of a multipedigree and its
//            pedigrees.
//                                                                          
//----------------------------------------------------------------------------
//
struct multipedigree_lod_score
{
  multipedigree_lod_score(const double& lod = QNAN);

  void  write(ostream& out);
  
  double  lod_score;
  vector<pedigree_lod_score>  ped_lod_scores;
};


//----------------------------------------------------------------------------
//  Class:    lods_result
//                                                                          
//  Purpose:  base class for lod score tasks.
//                                                                          
//----------------------------------------------------------------------------
//
struct lods_result : public task_result
{
  virtual ~lods_result() = 0;
  
  static double  get_subpedigree_lod_score(const multipedigree_lod_score& mpls,
                                           const string& pedigree, const string& member);    

  string trait;
  string marker;
};


//----------------------------------------------------------------------------
//  Class:    non_ss_lods_result
//                                                                          
//  Purpose:  represents a lod score result for a sex averaged recombination
//            fraction.
//                                                                          
//----------------------------------------------------------------------------
//
struct non_ss_lods_result : public lods_result
{
  non_ss_lods_result();

  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;
  void  write_family_detail(ostream& out, const string& pedigree, const string& member) const;
  void  write_vc_matrix(ostream& out) const;
  
  vector<pair<double, multipedigree_lod_score> >   lod_scores;
};

//----------------------------------------------------------------------------
//  Class:    ss_lods_result
//                                                                          
//  Purpose:  represents a lod score result for sex specific recombination
//            fractions.
//                                                                          
//----------------------------------------------------------------------------
//
struct ss_lods_result : public lods_result
{
  ss_lods_result();

  void  write_summary(ostream& out) const;
  void  write_detail(ostream& out) const;  
  void  write_family_detail(ostream& out, const string& pedigree, const string& member) const;  
  void  write_vc_matrix(ostream& out) const;
  
  vector<pair<theta_pair, multipedigree_lod_score> >   lod_scores;
};

#include "lodlink/lods_results.ipp"

}
}

#endif
