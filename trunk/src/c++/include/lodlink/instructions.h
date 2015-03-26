#ifndef LODLINK_INSTRUCTIONS_H
#define LODLINK_INSTRUCTIONS_H
//============================================================================
// File:      instructions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   9/5/2 created         -djb
//                                                                          
// Notes:     Defines struct, instructions, for holding instructions
//            supplied by the user in a parameter file LODLINK analysis block.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <set>
#include <map>
#include <string>
#include <limits>
#include <ostream>
#include "rped/rped.h"
#include "fped/fped.h"
#include "maxfun/sub_model.h"
#include "lodlink/output.h"

using std::ostream;
using std::endl;

namespace SAGE
{

namespace LODLINK
{

//----------------------------------------------------------------------------
//  Struct:    theta_pair
//                                                                          
//  Purpose:  Pair of recombination fractions for sex specific case.
//                                                                          
//----------------------------------------------------------------------------
//
struct theta_pair
{
  theta_pair(double ml_theta = QNAN, double fml_theta = QNAN);
  
  double  male_theta;
  double  female_theta;
};

bool      operator<(const theta_pair& left, const theta_pair& right);
bool      operator==(const theta_pair& left, const theta_pair& right);
ostream&  operator<<(ostream& out, const theta_pair& p);

extern const double THETAS[];
extern const theta_pair  THETA_PAIRS[];

extern const double* FIRST_THETA;
extern const double* LAST_THETA;
extern const theta_pair* FIRST_THETA_PAIR;
extern const theta_pair* LAST_THETA_PAIR;

typedef set<string>  group;

//----------------------------------------------------------------------------
//  Struct:    instructions
//                                                                          
//  Purpose:  repository for LODLINK analysis information.
//                                                                          
//----------------------------------------------------------------------------
//
struct instructions
{
  typedef RPED::RefMultiPedigree::member_const_pointer member_const_pointer;

  instructions(cerrorstream& errors = sage_cerr);
  void  reset();
  void  reset_thetas();
  
  void  write(ostream& out) const;
    void  write_groups(ostream& out) const;
    void  write_thetas(ostream& out) const;
  
  enum linkage_option  { MARKER, TRAIT };
  
  static string  linkage_option_2_string(linkage_option l);
  
  string  link_option() const;
  
  bool  sex_specific() const;
  
  // Data members.
  string  file_name_root;
  string  title;
  
  linkage_option  linkage;
  string  trait;  // Modeled trait or marker against wh. linkage is tested.
  
  bool  linkage_test;
  bool  linkage_sex_specific;
  bool  linkage_homog;                        // Assume homogeneity.
  bool  smiths_test;
  bool  smiths_sex_specific;
  bool  mortons_test;
  bool  mortons_sex_specific;
  std::map<string, group>  groups;                 // For Morton's homogeneity test.
  
  std::vector<theta_pair>  male_female_thetas;
  std::vector<double>      average_thetas;

  bool  genotypes;
  bool  genotypes_sex_specific;
  
  bool valid;
};

ostream& operator<<(ostream& out, const instructions& m);

#include "lodlink/instructions.ipp"

}
}

#endif
