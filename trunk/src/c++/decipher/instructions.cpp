//============================================================================
// File:      instructions.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   4/2/4 - created.                         djb
//                                                                          
// Notes:     Implementaion of instructions class.   
//               
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/instructions.h"

namespace SAGE
{

namespace DECIPHER
{

ostream& 
operator <<(ostream& out, const instructions::value& v)
{
  out << v.dbl << ", " << v.str;

  return  out;
}


ostream&
operator <<(ostream& out, const instructions& instr)
{
  out << std::boolalpha;
  out << "valid:                            " << instr.valid               << "\n"
      << "Output root name:                 " << instr.file_name_root      << "\n"      
      << "Title:                            " << instr.title               << "\n"
      << "Regions:                          " << instr.regions             << "\n"
      << "Epsilon:                          " << instr.epsilon             << "\n"                                                  
      << "Starting points:                  " << instr.starting_points     << "\n"
      << "Dump:                             " << instr.dump                << "\n"
      << "Dump cutoff:                      " << instr.dump_cutoff         << "\n"
      << "Sliding window:                   " << instr.sliding_window      << "\n"
      << "Window width:                     " << instr.window_width        << "\n"  
      << "Four gamete blocks:               " << instr.four_gamete_rule    << "\n"
      << "Four gamete threshold:            " << instr.fg_threshold        << "\n"
      << "LD blocks:                        " << instr.ld_blocks           << "\n"
      << "LD threshold:                     " << instr.ld_threshold        << "\n"      
      << "MAF filtration:                   " << instr.maf_filter          << "\n"
      << "MAF filter threshold              " << instr.maf_threshold       << "\n"        
      << "Partitions:                       " ; 
   
  assert(instr.partitions.size() == PARTITION_COUNT);
  if(! (instr.partitions[0].valid() || instr.partitions[1].valid()))
  {
    out << "none\n";
  }
  else
  {
    out << instr.partitions;
  }   

  out << "Family rep:                       ";
  if(instr.family_rep == (size_t)(-1))
  {
    out << "none\n";
  }
  else
  {
    out << instr.family_rep << "\n";
  }
  
  string  family_rep_type = instructions::field_type_2_string(instr.rep_field_type);
  out << "Family rep type:                  " << family_rep_type << endl;  
                                                
  out << "Family rep value:                 " << instr.family_rep_value              << "\n"
      << "Analysis unit:                    " << instructions::unit_2_string(instr.analysis_unit)  << "\n"
      << "Pool size:                        " << instr.pool_size                     << "\n"
      << "Pool size trait:                  " << instr.pool_size_trait               << "\n"
      << "Pop freq:                         " << instr.pop_freq                      << "\n"
      << "Freq cutoff:                      " << instr.freq_cutoff                   << "\n"
      << "All possible combinations:        " << instr.all_possible_diplotypes       << "\n"
      << "All possible combinations table:  " << instr.all_possible_diplotypes_table << "\n"      
      << "Most likely combinations:         " << instr.most_likely_diplotypes        << "\n"
      << "Likely cutoff:                    " << instr.likely_cutoff                 << "\n"
      << "All possible haplotypes:          " << instr.all_possible_haplotypes       << "\n"      
      << "Likelihood ratio test:            " << instr.likelihood_ratio_test         << "\n"
      << "Compute empirical p-value:        " << instr.compute_empirical_pvalue      << "\n"
      << "Permutations:                     " << instr.permutations                  << "\n"
      << "Seed:                             " << instr.seed                          << "\n"
      << "Minimum permutations:             " << instr.min_permutations              << "\n"      
      << "Maximum permutations:             " << instr.max_permutations              << "\n"
      << "Width:                            " << instr.width                         << "\n"
      << "Confidence:                       " << instr.confidence                    << endl;
      
  vector<instructions::pool_locus>::const_iterator  pl_iter     = instr.pool_loci.begin();
  vector<instructions::pool_locus>::const_iterator  pl_end_iter = instr.pool_loci.end();
  for(; pl_iter != pl_end_iter; ++pl_iter)
  {
    out << *pl_iter << endl;
  }
  
  out << "loci (for pools only):\n " << instr.loci << endl;

  return  out;
}

ostream&
operator <<(ostream& out, const vector<instructions::partition_data>& partitions)
{
  size_t  partition_count = partitions.size();
  for(size_t p = 0; p < partition_count; ++p)
  {
    out << "\npartition\n"
        << "field type:             " << instructions::field_type_2_string(partitions[p].type) << "\n"
        << "index:                  " << partitions[p].field_index << "\n"
        << "subpopulations:\n"        << partitions[p].sub_pops
        << endl;
  }

  return  out;
}

ostream&
operator <<(ostream& out, const map<string, instructions::value>& sub_pops)
{
  map<string, instructions::value>::const_iterator  sp_iter     = sub_pops.begin();
  map<string, instructions::value>::const_iterator  sp_end_iter = sub_pops.end();
  for(; sp_iter != sp_end_iter; ++sp_iter)
  {
    out << "sub pop " << sp_iter->first 
        << "  value " << sp_iter->second << endl;
  }
  
  return  out;
}

ostream&
operator <<(ostream& out, const locus_group& loci)
{
  size_t  locus_count = loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    out << loci[l].first << " : " << (loci[l].second)->name() << " " << endl;
    out << "alleles  ";
    
    MLOCUS::allele_iterator  a_iter     = (loci[l].second)->allele_begin();
    MLOCUS::allele_iterator  a_end_iter = (loci[l].second)->allele_end();
    for(; a_iter != a_end_iter; ++a_iter)
    {
      out << a_iter->name() << " ";
    }
    
    out << endl;
  }
  
  out << endl;

  return  out;
}

ostream& 
operator <<(ostream& out, const instructions::pool_locus::allele& a)
{
  out << a.name << "  " << a.index << endl;
  
  return  out;
}

ostream& 
operator <<(ostream& out, const instructions::pool_locus& pl)
{
  out << pl.name << endl;
  
  set<instructions::pool_locus::allele>::const_iterator      a_iter = pl.alleles.begin();
  set<instructions::pool_locus::allele>::const_iterator  a_end_iter = pl.alleles.end();  
  for(; a_iter != a_end_iter; ++a_iter)
  {
    out << "  " << *a_iter;
  }
  
  out << endl;
  
  return  out;
}

ostream& 
operator <<(ostream& out, const instructions::region_data& rd)
{
  out << "\n" << rd.name << endl;
  out << "first marker index: " << rd.first << endl;
  out << "last marker index: " << rd.last << endl;
  out << "valid: " << (rd.valid ? "yes" : "no") << endl;
  out << "x_linked: " << (rd.x_linked ? "yes" : "no") << endl;
  out << rd.loci << endl;  
  
  return  out;
}

ostream& 
operator <<(ostream& out, const list<instructions::region_data>& rds)
{
  list<instructions::region_data>::const_iterator  r_iter     = rds.begin();
  list<instructions::region_data>::const_iterator  r_end_iter = rds.end();
  for(; r_iter != r_end_iter; ++r_iter)
  {
    out << *r_iter;
  }
  
  out << endl;
  
  return out;
}


}
}

