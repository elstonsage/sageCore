#ifndef DECIPHER_INSTRUCTIONS_H
#define DECIPHER_INSTRUCTIONS_H

//============================================================================
// File:      instructions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/2/4 created         -djb
//                                                                          
// Notes:     Defines struct, instructions, for holding instructions
//            supplied by the user in a parameter analysis block.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "rped/rped.h"
#include "decipher/definitions.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace RPED;

const size_t  PARTITION_COUNT = 2;

const size_t  STARTING_POINTS_DEFAULT = 10;
const double  EPSILON_DEFAULT = .00001;
const double  FREQ_CUTOFF_DEFAULT = .001;
const double  DUMP_CUTOFF_DEFAULT = .001;
const double  LIKELY_CUTOFF_DEFAULT = .05;
const size_t  PERMUTATIONS_DEFAULT = 0;
const int     SEED_DEFAULT = 0;
const size_t  MIN_PERMUTATIONS_DEFAULT = 50;
const size_t  MAX_PERMUTATIONS_DEFAULT = 10000;
const double  WIDTH_DEFAULT = .2;
const double  CONFIDENCE_DEFAULT = .95;
const size_t  POOL_SIZE_DEFAULT = 2;
const size_t  WINDOW_WIDTH_DEFAULT = 3;
const double  FG_THRESHOLD_DEFAULT = .01;
const double  LD_THRESHOLD_DEFAULT = .8;
const double  MAF_THRESHOLD_DEFAULT = .1;

//----------------------------------------------------------------------------
//  Struct:    instructions
//                                                                          
//  Purpose:  repository for DECIPHER analysis information.
//                                                                          
//----------------------------------------------------------------------------
//
struct instructions
{
  enum  value_type { TRAIT, STRING };
  enum  unit       { EACH_INDIVIDUAL, FAMILY_REP, FAMILY_FOUNDERS, POOL };
  enum  field_type { CONTINUOUS_FIELD, BINARY_FIELD, STRING_FIELD, INVALID_FIELD };  
  
  static string  field_type_2_string(instructions::field_type type);
  static string  unit_2_string(instructions::unit u);  
  
  struct value
  {
    value();
    void  reset();  // Set dbl and str to QNAN and "" respectively.
    bool  empty();
    bool  operator ==(const value& other) const;
  
    double  dbl;
    string  str;
  };
  
  struct partition_data
  {
    partition_data();
    void  reset();
    bool  valid() const;
  
    field_type          type;
    size_t              field_index;
    std::map<string, value>  sub_pops;  
  };
  
  // - If valid ( != to MAX_SIZE) first and last locus indices are given,
  //   they define the region, otherwise name refers to a region in the
  //   genome_description file.
  //
  struct region_data     // Defines 'haplotyping domains.'
  {
    region_data(const string& nm = "", size_t f = MAX_SIZE, 
                size_t l = MAX_SIZE, bool v = true, bool xl = false);
  
    string  name;  
    size_t  first;
    size_t  last;
    bool  valid;
    bool  x_linked;
    locus_group  loci;
  };
  
  static bool  invalid_region(const region_data& region);

  instructions(cerrorstream& errors = sage_cerr);
  void  reset();
  void  reset_family_rep();

  
  void  write(ostream& out) const;
  
  // Data members.
  string  file_name_root;
  string  title;
  
  list<region_data>  regions;
  
  double  epsilon;                // EM algorithm convergence criterion.
  size_t  starting_points;        // Number of sets of starting conditions for
                                  // which em algorithm is run.
  // Block sub-block
  bool  sliding_window;
    size_t  window_width;
  bool  four_gamete_rule;
    double  fg_threshold;    
  bool  ld_blocks;
    double  ld_threshold;        
  
  // Filter sub-block
  bool  maf_filter;
    double  maf_threshold;

  // Data sub-block
  unit  analysis_unit;  
  
  // - partitions[0] is 'outer' partition (specified second by user)
  //   partitions[1] is 'inner' partition (specified first by user)
  //
  vector<partition_data>  partitions;
  
  size_t      family_rep;
    field_type  rep_field_type;
    value       family_rep_value;

  
  // Pools sub-block
  size_t  pool_size;
  size_t  pool_size_trait;      // Allows specification of pool size on a per record basis.
  
  struct pool_locus
  {
    bool  operator ==(const pool_locus& other) const;
    
    string  name;
    
    struct allele
    {
      allele(const string& nm, size_t idx = (size_t)(-1));
      bool  operator <(const allele& other) const;
      
      string  name;
      size_t  index;
    };
    
    set<allele>  alleles;  
  };
  
  vector<pool_locus>  pool_loci;
  locus_group  loci;                // Used only for pools.
 
  // Tasks sub-block
  bool  pop_freq;
    double  freq_cutoff;
  bool  dump;
    double  dump_cutoff;
  bool  all_possible_diplotypes;
  bool  all_possible_diplotypes_table;
  bool  most_likely_diplotypes;
    double  likely_cutoff;
    bool  all_possible_haplotypes;    // Hidden feature.

  bool  likelihood_ratio_test;
  bool  compute_empirical_pvalue;
    size_t  permutations;
    int  seed;                      // Hidden feature.
    size_t  min_permutations;       // Hidden feature.
    size_t  max_permutations;
    double  width;
    double  confidence;
  
 
  bool  valid;
};

size_t  valid_partition_count(const vector<instructions::partition_data>& partitions);

ostream& operator <<(ostream& out, const instructions& instr);
ostream& operator <<(ostream& out, const instructions::value& v);
ostream& operator <<(ostream& out, const vector<instructions::partition_data>& partitions);
ostream& operator <<(ostream& out, const std::map<string, instructions::value>& sub_pops);
ostream& operator <<(ostream& out, const locus_group& loci);
ostream& operator <<(ostream& out, const instructions::pool_locus::allele& a);
ostream& operator <<(ostream& out, const instructions::pool_locus& pl);
ostream& operator <<(ostream& out, const instructions::region_data& rd);
ostream& operator <<(ostream& out, const list<instructions::region_data>& rds);

#include "decipher/instructions.ipp"

}
}

#endif
