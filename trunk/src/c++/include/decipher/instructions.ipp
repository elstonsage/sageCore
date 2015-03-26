//============================================================================
// File:      instructions.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/2/4 created        -djb
//                                                                          
// Notes:     Inline implementation of struct, instructions.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



inline size_t
valid_partition_count(const vector<instructions::partition_data>& partitions)
{
  size_t  valid_count = 0;
  
  size_t  total_count = partitions.size();
  assert(total_count == 2);
  
  for(size_t p = 0; p < total_count; p++)
  {
    if(partitions[p].valid())
    {
      valid_count++;
    }
  }
  
  // - 'Inner' partition should be the valid one.
  //
  if(valid_count == 1)
  {
    assert(! partitions[0].valid());
  }
  
  return  valid_count;
}

//============================================================================
// IMPLEMENTATION:  instructions::value
//============================================================================
//
inline
instructions::value::value()
{
  reset();
}

inline void
instructions::value::reset()
{
  dbl = QNAN;
  str = "";
}

inline bool
instructions::value::empty()
{
  return  SAGE::isnan(dbl) && str.empty();
}

inline bool
instructions::value::operator ==(const value& other) const
{
  return  (SAGE::isnan(dbl) && SAGE::isnan(other.dbl) ||
           dbl == other.dbl                 ) &&
          str == other.str;
}

//============================================================================
// IMPLEMENTATION:  instructions::partition_data
//============================================================================
//
inline
instructions::partition_data::partition_data()
{
  reset();
}

inline void
instructions::partition_data::reset()
{
  type = instructions::INVALID_FIELD;
  field_index = (size_t)(-1);
  sub_pops.clear();
}

inline bool
instructions::partition_data::valid() const
{
  bool  field_index_valid = field_index != (size_t)(-1);
  bool  has_sub_pops = ! sub_pops.empty();
  bool  type_valid = type == instructions::CONTINUOUS_FIELD ||
                     type == instructions::BINARY_FIELD     ||
                     type == instructions::STRING_FIELD       ;
                     
  return  field_index_valid && has_sub_pops && type_valid;
}


//============================================================================
// IMPLEMENTATION:  instructions::region_data
//============================================================================
//
inline
instructions::region_data::region_data(const string& nm, size_t f, size_t l, bool v, bool xl)
      : name(nm), first(f), last(l), valid(v), x_linked(xl)
{}


//============================================================================
// IMPLEMENTATION:  instructions::pool_locus
//============================================================================
//
inline bool
instructions::pool_locus::operator ==(const pool_locus& other) const
{
  return  name == other.name;
}

//============================================================================
// IMPLEMENTATION:  instructions::pool_locus::allele
//============================================================================
//
inline
instructions::pool_locus::allele::allele(const string& nm, size_t idx)
      : name(nm), index(idx)
{}

inline bool
instructions::pool_locus::allele::operator <(const allele& other) const
{
  return  name < other.name;
}

//============================================================================
// IMPLEMENTATION:  instructions
//============================================================================
//
inline
instructions::instructions(cerrorstream& errors)
      : partitions(PARTITION_COUNT)
{
  reset();
}

inline void
instructions::reset()
{
  file_name_root = "";
  title          = "";
  
  regions.clear();
  
  epsilon = EPSILON_DEFAULT;
  starting_points = STARTING_POINTS_DEFAULT;

  sliding_window = false;
  window_width = WINDOW_WIDTH_DEFAULT;
  four_gamete_rule = false;
  fg_threshold = FG_THRESHOLD_DEFAULT;
  ld_blocks = false;
  ld_threshold = LD_THRESHOLD_DEFAULT;  

  maf_filter = false;
  maf_threshold = MAF_THRESHOLD_DEFAULT;  
  
  analysis_unit = FAMILY_FOUNDERS;
  
  size_t  partition_count = partitions.size();
  for(size_t p = 0; p < partition_count; ++p)
  {
    partitions[p].reset();
  }
  
  reset_family_rep();  
  
  pool_size = POOL_SIZE_DEFAULT;
  pool_size_trait = (size_t)(-1);
  pool_loci.clear();
  loci.clear();
  
  pop_freq                      = true;
  freq_cutoff                   = FREQ_CUTOFF_DEFAULT;
  dump                          = false;
  dump_cutoff                   = DUMP_CUTOFF_DEFAULT;
  all_possible_diplotypes       = false;
  all_possible_diplotypes_table = false;
  most_likely_diplotypes        = false;
  likely_cutoff                 = LIKELY_CUTOFF_DEFAULT;
  all_possible_haplotypes       = false;
  likelihood_ratio_test         = false;
  compute_empirical_pvalue      = false;
  permutations                  = PERMUTATIONS_DEFAULT;
  seed                          = SEED_DEFAULT;
  min_permutations              = MIN_PERMUTATIONS_DEFAULT;
  max_permutations              = MAX_PERMUTATIONS_DEFAULT;
  width                         = WIDTH_DEFAULT;
  confidence                    = CONFIDENCE_DEFAULT;
  
  // - Note:  not *really* valid until parser::init_parse() is called and user supplied
  //   regions are set.
  //
  valid = true;
}

inline void
instructions::reset_family_rep()
{
  family_rep  = (size_t)(-1);
  rep_field_type = INVALID_FIELD;
  family_rep_value.reset();
}


inline string
instructions::field_type_2_string(instructions::field_type type)
{
  string  result = "";
  switch(type)
  {
    case CONTINUOUS_FIELD:
      result = "continuous";
      break;
      
    case BINARY_FIELD:
      result = "binary";
      break;
      
    case STRING_FIELD:
      result = "string";
      break;
      
    case INVALID_FIELD: 
      result = "invalid";
      break;
      
    default:
      assert(false);
  }
  
  return  result;
}

inline string
instructions::unit_2_string(instructions::unit u)
{
  string  result = "";
  switch(u)
  {
    case EACH_INDIVIDUAL:
      result = "each individual";
      break;
      
    case FAMILY_REP:
      result = "family rep";
      break;
      
    case FAMILY_FOUNDERS:
      result = "family founders";
      break;
      
    case POOL:
      result = "pool";
      break;      
      
    default:
      assert(false);
  }
  
  return  result;
}

inline bool
instructions::invalid_region(const region_data& region)
{
  return  ! region.valid;
}


