#ifndef DECIPHER_ANALYSIS_H
#define DECIPHER_ANALYSIS_H

//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/8/4 - created.                                   djb
//                                                                          
// Notes:     declaration of a analysis class. 
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "numerics/functions.h"
#include "decipher/partitioning.h"
#include "decipher/pool.h"
#include "decipher/rebuilt.h"
#include "decipher/sub_pop_shuffler.h"
#include "decipher/output.h"

namespace SAGE
{

namespace DECIPHER
{

const string  double_line = "\n  ===========================================================================\n";

//----------------------------------------------------------------------------
//  Class:    analysis
//                                                                          
//  Purpose:  execute and write results of a set of user instructions 
//            corresponding to a parameter file decipher analysis block.
//                                                                          
//----------------------------------------------------------------------------
//
class analysis
{
  public:
  
    analysis(APP::Output_Streams& streams, ostream& detail_file, ostream& summary_file,
             ostream& dump_file, 
             const RefMultiPedigree& mped, const genome_description* gen, 
             const instructions& instr);
    
    void  analyze();
    void  report_inconsistencies(const family_generator::iterator& iter); 
    const string&  title() const;
    
  private:
    friend class sliding_window_blocks;
    friend class four_gamete_rule_blocks;
    friend class linkage_disequilibrium_blocks;
  
    // ----- Construction -----

    void      build_filtered_mped();
      
    // ----- Preliminaries -----

    void  set_test_seed();
    bool  maximization_needed() const;
    bool  nothing_to_do() const;
    void  write_user_options();    
   
    // ----- Non-pool Analysis -----
    
      void  do_non_pool_analysis();
 
      // ----- Locus Filtration -----
   
      void  build_current_region(const instructions::region_data& region);
        bool  low_maf_freq(const inheritance_model* locus);

      // ----- Blocks -----

      bool  no_block_options_specified() const;

      void  analyze_by_sliding_window();
      void  analyze_by_four_gamete_rule_blocks();
      void  analyze_by_ld_blocks();
      
      const locus_group  build_block(const block_vector& bv, size_t index) const;
      void  analyze_block(const locus_group& loci, 
                          const string& region = "", size_t block_number = 0);
                        
      void  write_done();
    
    // ----- Tasks -----

    //   M = phenotype map
    //   V = value of innermost map
    //
    template<typename M, typename V>    
    void   do_tasks(const std::map<string, std::map<string, V> >& members, 
                          const locus_group& loci, 
                          const string& region, size_t block_number = 0);        
      p_values  do_tests(const vector<const member_em_phenotype_map*>& phenotype_maps, const locus_group& loci);
        double    do_permutation_test(double test_statistic,
                                      double whole_ln_likelihood,
                                      sub_pop_shuffler& suf,
                                      const vector<const member_em_phenotype_map*>& phenotype_maps,
                                      const locus_group& loci,
                                      size_t& total_permutations);    

    // Data members.
    APP::Output_Streams&  my_streams;
    cerrorstream&  my_errors;
    ostream&  my_messages;
    ostream&  my_detail_file;
    ostream&  my_summary_file;
    ostream&  my_dump_file;
    const RefMultiPedigree&  my_mped;
    const genome_description*  my_genome_description;
    const instructions&  my_instructions;
    
    locus_group  my_current_region;
    string  my_current_region_name;
    vector<ld_record>  my_current_region_lds;
    
    FilteredMultipedigree  my_filtered_mped;
    partitioner  my_partitioner;
    output_state  my_output_state;
};


//----------------------------------------------------------------------------
//  Class:    blocks
//                                                                          
//  Purpose:  base class for block determining classes. 
//                                                                          
//----------------------------------------------------------------------------
//
class blocks
{
  public:
    blocks(analysis& a);
    virtual ~blocks() = 0;    
    
    const block_vector&  operator()() const;
  
  protected:
  
    // Data members
    analysis&  my_analysis;
    block_vector  my_blocks;
};


//----------------------------------------------------------------------------
//  Class:    sliding_window_blocks
//                                                                          
//  Purpose:  determine sliding window blocks. 
//                                                                          
//----------------------------------------------------------------------------
//
class sliding_window_blocks : public blocks
{
  public:
    sliding_window_blocks(analysis& a);
};


//----------------------------------------------------------------------------
//  Class:    four_gamete_rule_blocks
//                                                                          
//  Purpose:  determine blocks using the four gamete rule. 
//                                                                          
//----------------------------------------------------------------------------
//
class four_gamete_rule_blocks : public blocks
{
  public:
    four_gamete_rule_blocks(analysis& a);
    
  private:
    bool  recomb(size_t begin, size_t end);
    size_t  common_freqs(const base_em_phenotype_map* pheno_map, double threshold) const;
};


//----------------------------------------------------------------------------
//  Class:    linkage_disequilibrium_blocks
//                                                                          
//  Purpose:  determine blocks using linkage disequilibrium. 
//                                                                          
//----------------------------------------------------------------------------
//
class linkage_disequilibrium_blocks : public blocks
{
  public:
    linkage_disequilibrium_blocks(analysis& a);
    const vector<ld_record>&  ld_data() const;
    
  private:
    bool  ld(size_t locus1, size_t locus2);
    double  lewontins_ld(const base_em_phenotype_map* em_map);
    void  check_allele_frequencies(const base_em_phenotype_map* em_map);
    void  test_allele_frequency(const MLOCUS::inheritance_model* locus, double allele_frequency);
    
    // Data members
    vector<ld_record>  my_lds;
    string  my_last_freq_warning;        // Holds name of last locus for which a frequency warning was issued.
};

void  write_blocks(ostream& out, const block_vector& bv);

#include "decipher/analysis.ipp"
}
}

#endif

