#ifndef DECIPHER_PARSER_H
#define DECIPHER_PARSER_H

//============================================================================
// File:      parser.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/2/4 created                     - djb
//                                                                          
// Notes      Defines class, parser, for parsing DECIPHER block of
//            a SAGE parameter file.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "rped/genome_description.h"
#include "app/aparser.h"
#include "LSF/LSF.h"
#include "LSF/Attr.h"
#include "LSF/parse_ops.h"
#include "mped/sp.h"
#include "mped/mp.h"
#include "decipher/instructions.h"
#include "decipher/input.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace RPED;
using namespace APP;
using namespace MPED;

const int  MAX_INT = 32767;
const string  DUMMY_SUB_POP_NAME = "family_rep_value";
const double  MAX_FG_THRESHOLD = 0.25;
const double  MAX_LD_THRESHOLD = 1.0;
const double  MAX_MAF_THRESHOLD = .5;

//----------------------------------------------------------------------------
//  Class:    parser
//                                                                          
//  Purpose:  extract information from decipher block of SAGE parameter
//            file.
//                                                                          
//----------------------------------------------------------------------------
//
class parser : public BasicParser
{
  public:
    
    typedef  instructions::partition_data  partition_data;
    typedef  instructions::value  value;
    typedef  instructions::pool_locus  pool_locus;
    
    // Constructor/destructor.
    parser(const RefMultiPedigree* mp, const decipher_data& data,
           ostream& messages = std::cout,
           SAGE::cerrorstream& errors = SAGE::sage_cerr,
           genome_description* gd = 0);
  
    virtual void  parse_symbols(const SymbolTable* syms);
    virtual void  parse_parameter(const LSFBase* param);
    virtual void  parse_test_parameter_section(const LSFBase* params);
    virtual void  parse_test_parameter(const LSFBase* param);
    void  parse(const LSFBase* params);
    
    const instructions&  user_instructions() const;  
    size_t  analysis_id() const;
    
  private:
  
    void  init_parse();
    void  print_header();
    void  print_footer();

    // General ancillary functions.
    void  reset();
    
    // Constraint checking.
    void  check_for_required_sub_pops();
    void  check_dump_conditions();
    
    void  check_for_regions();
    void  validate_regions();
      void  build_loci(instructions::region_data& region); 
        bool  check_first_last_order(const instructions::region_data& region);   
      void  check_locus_count(instructions::region_data& region);  
      void  check_for_codominance(instructions::region_data& region);  
      void  check_x_linkage(instructions::region_data& region);     
      void  check_related_x_linkage(instructions::region_data& region); 
      void  check_snp_requirement(instructions::region_data& region);

    void  build_loci_for_pools();
    void  check_pool_locus_count();
    
    void  check_sliding_window();
    void  check_four_gamete_rule();
    void  check_ld();
    void  check_maf_filter();
    void  check_four_gamete_threshold();
    
    void  check_for_mld_file();   

    // Parsing functions.
    void  parse_out_attr(const LSFBase* param);
    void  parse_title(const LSFBase* param);
    
    void  parse_region(const LSFBase* param);
      size_t  marker_index(const string& marker_name) const;
      
    void  parse_epsilon(const LSFBase* param);
    void  parse_starting_points(const LSFBase* param);
    void  parse_dump(const LSFBase* param);
    void  parse_seed(const LSFBase* param);      // For testing.  Undocumented.
    
    class duplication
    {
      public:
        duplication(const pair<string, value>& sub_pop);
        
        bool  operator ()(const pair<string, value>& elem);
        
      private:
        const pair<string, value>  my_sub_pop;
    };
    
    // Filters.
    void  parse_filters(const LSFBase* param);
      void  parse_maf_filter(const LSFBase* param);

    // Blocks.
    void  parse_blocks(const LSFBase* param);
      void  parse_sliding_window(const LSFBase* param);
      void  parse_four_gamete_rule(const LSFBase* param);
      void  parse_ld(const LSFBase* param); 
    
    void  parse_data(const LSFBase* param);
      void  parse_analysis_unit(const LSFBase* param);          
    
      // Pools.
      void  parse_pools(const LSFBase* param);
        void  parse_pool_size(const LSFBase* param);  
        void  parse_pool_size_trait(const LSFBase* param);    
        void  parse_locus(const LSFBase* param);      
          void  parse_allele(const LSFBase* param, pool_locus& pl);
          void  parse_last_allele(const LSFBase* param, const pool_locus& pl, pool_locus::allele& a);
            bool  add_allele(pool_locus& pl, const pool_locus::allele& a, const string& trait_name = "");
          
      // Partitions.
      void  parse_partition(const LSFBase* param);
        void  get_field_type_and_index(const string& partition_name, partition_data& parse_results);
        bool  partition_repeated(const partition_data& parse_results);
        void  parse_sub_pops(const LSFBase* param, partition_data& parse_results);
          void  parse_continuous_sub_pop(const LSFBase* param, partition_data& parse_results);                
          void  parse_binary_sub_pop(const LSFBase* param, partition_data& parse_results);        
          void  parse_string_sub_pop(const LSFBase* param, partition_data& parse_results);
            static bool  new_sub_pop(const pair<string, instructions::value>&, const partition_data& parse_results);                  
        void  fill_continuous_sub_pops(partition_data& parse_results);
        void  fill_binary_sub_pops(partition_data& parse_results);
        void  fill_string_sub_pops(partition_data& parse_results);        
        
      // Family representatives.
      void  parse_family_rep(const LSFBase* param);
        void  parse_continuous_family_rep_value(const LSFBase* param, partition_data& parse_results);
        void  parse_binary_family_rep_value(const LSFBase* param, partition_data& parse_results);
        void  parse_string_family_rep_value(const LSFBase* param, partition_data& parse_results);

    
    void  parse_tasks(const LSFBase* param);
      void  parse_pop_freq(const LSFBase* param);
      void  parse_all_possible_diplotypes(const LSFBase* param);      
      void  parse_all_possible_diplotypes_table(const LSFBase* param);
      void  parse_most_likely_diplotypes(const LSFBase* param);      
      void  parse_all_possible_haplotypes(const LSFBase* param);            
      void  parse_likelihood_ratio_test(const LSFBase* param);
      void  parse_compute_empirical_pvalue(const LSFBase* param);
      
    enum  trait_usage { FOR_FAMILY_REP, FOR_PARTITIONING, 
                        FOR_POOL_SIZE, FOR_ALLELE_PCT, NONE };
    static string  trait_usage_2_string(trait_usage);

    struct trait_info
    {
      trait_info(const string& n = "", parser::trait_usage u = NONE);
    
      string  name;
      trait_usage  usage;
    };

    bool  trait_registered(size_t t) const;
    bool  string_registered(size_t s) const;
    
    static void  write_trait_info(ostream& out, const trait_info& ti);
    void  write_trait_registry(ostream& out) const;  
    
    // Data members.
    instructions  my_instructions;
    
    size_t  my_analysis_id;       
    const RefMultiPedigree*  my_ref_mped;
    const decipher_data&  my_data;
    ostream&  my_messages;
    genome_description*  my_genome;
    
    vector<MLOCUS::inheritance_model>  my_pool_loci;
    
    map<size_t, trait_info>  my_trait_registry;
    map<size_t, trait_info>  my_string_registry;
};

#include "decipher/parser.ipp"

}
}

#endif

