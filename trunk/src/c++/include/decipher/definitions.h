#ifndef DECIPHER_DEFINITIONS_H
#define DECIPHER_DEFINITIONS_H

//============================================================================
// File:      definitions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/28/5 - created.                                   djb
//                                                                          
// Notes:     Definitions used in multiple files.
//
//
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <ostream>
#include <sstream>
#include <utility>

#include "mlocus/imodel.h"
#include "fped/fped.h"
#include "globals/SAGEConstants.h"
#include "numerics/mt.h"

#define NULL_STRING string("")

namespace SAGE
{

namespace DECIPHER
{

using namespace std;

const size_t  MAX_SIZE = (size_t)(-1);
extern MersenneTwister  rand_src;

typedef FPED::FilteredMultipedigree::member_const_pointer        member;
typedef FPED::FilteredMultipedigree::member_const_iterator       member_iterator;
typedef FPED::FilteredMultipedigree::subpedigree_const_iterator  subpedigree_iterator;
typedef FPED::FilteredMultipedigree::pedigree_const_iterator     pedigree_iterator;

typedef vector<pair<size_t, const MLOCUS::inheritance_model*> >  locus_group;

// - Contains begin and end indices into a locus_group.
//
typedef vector<pair<size_t, size_t> >  block_vector;

struct output_state
{
  output_state();
  
  bool  msg_interrupted;
  bool  allow_progress_msg;
  bool  first_block;
};

struct p_values
{
  p_values(const string& name = "", double asymp = QNAN, double empir = QNAN);

  string  outer_sub_pop_name;
  double  asymptotic;
  double  composite_ln_like;
  double  whole_ln_like;
  double  test_statistic;
  size_t  degrees_of_freedom;
  double  empirical;
  size_t  permutations;
};

struct ld_record
{
  ld_record(const string& m1, const string& m2, double linkage_disequilibrium);

  string  marker1;
  string  marker2;
  double  ld;
};

string  double_to_string(double value);    

#include "definitions.ipp"

}
}

#endif

