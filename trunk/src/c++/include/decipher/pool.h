#ifndef DECIPHER_POOL_H
#define DECIPHER_POOL_H

//============================================================================
// File:      pool.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   1/6/5 - created.                                   djb
//                                                                          
// Notes:     Classes for generating and storing phenotypes of pools
//            as used in the EM algorithm.
//
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "rped/genome_description.h"
#include "decipher/em.h"
#include "decipher/instructions.h"


namespace SAGE
{

namespace DECIPHER
{

// - 3 missing triallelic markers or 5 missing snp's w. pool size
//   size of 4.
//
const size_t MAX_EXPANDED_PHENO_SEQ_COUNT = 3375;

typedef map<size_t, double>  allele_weights;      // <allele_id, weight>
typedef map<size_t, size_t>  allele_counts;       // <allele_id, count>
typedef pair<const MLOCUS::inheritance_model*, allele_counts>  pool_pheno;     
typedef vector<pool_pheno>  pool_pheno_seq;

typedef vector<hap_seq>  pool_comb;

ostream&  operator <<(ostream& out, const pool_pheno pheno);
ostream&  operator <<(ostream& out, const pool_pheno_seq& seq);
void  write_allele_weights(ostream& out, const allele_weights& weights, const MLOCUS::inheritance_model* locus);
void  write_allele_counts(ostream& out, const allele_counts& counts, const MLOCUS::inheritance_model* locus);


//----------------------------------------------------------------------------
//  Class:    pool_pheno_seq_generator
//                                                                         
//  Purpose:  Provide iteration through phenotype sequences in pool
//            data.
//                                                                          
//----------------------------------------------------------------------------
//
class pool_pheno_seq_generator
{
  public:
  
    class iterator
    {
      public:
        enum reconcile { GREATEST, LEAST, RANDOM, OVER_UNDER };
      
        friend class pool_pheno_seq_generator;
      
        bool  operator ==(const iterator& other) const;
        bool  operator !=(const iterator& other) const;
        
        member  get_member() const;
        
        iterator  operator ++();
        pool_pheno_seq&  operator *();
        pool_pheno_seq*  operator ->();
        
        string  pool_name() const;
        size_t  phenotype_count() const;
        size_t  pool_size() const;        
      
      private:
        iterator(pool_pheno_seq_generator* generator, bool end = false);
        
        bool  build_seq();
          static  size_t  genotype_count(size_t n, size_t k);
          void  build_pool_phenotype(size_t l, allele_counts& counts, size_t pool_size);
            void  get_allele_weights(size_t l, allele_weights& weights);
            void  get_allele_counts(size_t l, const allele_weights& weights, 
                                      allele_counts& counts, size_t pool_size);
              size_t  get_pool_size() const;
              static size_t  round_2_integer(double num);
              void  reconcile_counts(size_t l, allele_counts& counts, 
                                            size_t pool_size, reconcile  method);
                static size_t  count_alleles(const allele_counts& counts);                                                                     
                static allele_counts::iterator  greatest(allele_counts& counts);
                static allele_counts::iterator  least(allele_counts& counts);                
                static allele_counts::iterator  random(allele_counts& counts);
                static allele_counts::iterator  over_under(allele_counts& counts, bool over);                
                  static const allele_counts::iterator&  random_pick(const vector<allele_counts::iterator>& candidates);
      
        // Data members.
        pool_pheno_seq  my_pool_pheno_seq;
        
        // - If all possibilities for missing data are considered, how many
        //   pheno_seqs are there?
        //
        size_t  my_expanded_pheno_seq_count;
        vector<member>::const_iterator  my_member_iter;
        vector<member>::const_iterator  my_member_end_iter;
        
        pool_pheno_seq_generator*  my_generator;        
    };
    
    // Constructor/destructor.
    pool_pheno_seq_generator(const vector<member>& members, const locus_group& loci,
                             cerrorstream& errors, ostream& messages, 
                             const instructions& instr, output_state& messages_ostate);
    
    cerrorstream&  errors();
    ostream&  messages();
    const instructions&  instr() const;
    const vector<member>&  members() const;
    bool& messages_interrupted();
    const locus_group&  loci() const;
    
    // Iteration.
    iterator  begin();
    iterator  end();
    
  private:
  
    // Data members.
    const vector<member>&  my_members;
    const locus_group&  my_loci;
    cerrorstream&  my_errors;
    ostream&  my_messages;
    const instructions&  my_instructions;
    output_state&  my_output_state;
};


//----------------------------------------------------------------------------
//  Class:    pool_comb_generator
//                                                                         
//  Purpose:  Generate combinations of haplotypes consistent w. phenotype
//            sequences.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class pool_comb_generator
{
  public:
    typedef vector<size_t>  allele_group;
    
    // For missing data.
    typedef pair<const MLOCUS::inheritance_model*, vector<allele_counts> >  multi_pool_pheno;
    typedef vector<multi_pool_pheno>  multi_pool_pheno_seq;
    typedef pair<const MLOCUS::inheritance_model*, const allele_counts*>  pool_pheno_ptrs;    
    
    static allele_group   counts_2_group(const allele_counts& counts);
    static allele_counts  group_2_counts(const MLOCUS::inheritance_model* locus, const allele_group& group);
    
    static ostream&  write_allele_group(ostream& out, const allele_group& group, 
                                             const MLOCUS::inheritance_model* locus);
  
    static void  generate(pool_pheno_seq_generator::iterator& seq_iter, vector<pool_comb>& combs);
    static void  write_comb(ostream& out, const pool_comb& comb, const pool_pheno_seq& seq);
    
  private:
    static void  expand_pheno_seq(pool_pheno_seq_generator::iterator& seq_iter, 
                                               multi_pool_pheno_seq& expanded_seq);
                                               
    static void  enumerate_genotypes(const MLOCUS::inheritance_model* locus, size_t pool_size, 
                                                                vector<allele_counts>& genotypes);
                                                                
    static void  enumerate_genotypes(const MLOCUS::inheritance_model* locus, size_t pool_size,
                                                                 const allele_group& genotype, 
                                                                 vector<allele_counts>& genotypes);
                                                                 
    static void  enumerate_sequences(multi_pool_pheno_seq::const_iterator mp_iter,
                                     multi_pool_pheno_seq::const_iterator& mp_end_iter,
                                     const vector<pool_pheno_ptrs>& sequence,
                                     vector<vector<pool_pheno_ptrs> >& sequences);
                                                         
    static pool_pheno_seq  make_dereferenced_sequence(const vector<pool_pheno_ptrs>& sequence);
  
    static void  init_comb(pool_comb& comb, const allele_counts& counts);
    static void  generate(pool_comb& comb, 
                          const pool_pheno_seq& seq, vector<pool_comb>& combs);
      static size_t  find_start(const pool_comb& comb);
      static size_t  get_extension_size(const pool_comb& comb, size_t start, size_t locus_count);
      static allele_counts  get_remaining_alleles(const pool_comb& comb, const pool_pheno_seq& seq,
                                                                  size_t start, size_t current_locus);
      static void  generate_allele_groups(const allele_counts& counts, size_t group_size, 
                                                          set<allele_group>& results);
        static void  generate_allele_groups(allele_group& group, size_t group_size,
                                            allele_group& partial_group, set<allele_group>& results);
};


//----------------------------------------------------------------------------
//  Class:    pool_em_phenotype_map
//                                                                         
//  Purpose:  Container for phenotypes of pools in the EM 
//            algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class pool_em_phenotype_map : public member_em_phenotype_map
{
  public:
  
    // Constructor/destructor.
    pool_em_phenotype_map(output_state& ostate,
                     APP::Output_Streams& streams,
                     const FilteredMultipedigree& mped,
                     const vector<member>& members, 
                     const instructions& instr, 
                     const locus_group& loci, 
                     const string& inner_sub_pop_name = "",
                     const string& outer_sub_pop_name = "" );

    const em_phenotype&  operator [](member m) const;
    set<member, member_order<member> >  members() const;
    const em_phenotype&  final_phenotype(member m) const;    
    
    const FilteredMultipedigree&  multipedigree() const;    

    void  dump(ostream& out) const;
    
  private:
  
    void  build();
    
    size_t  pheno_index(member m) const;        
    
    // Data members.
    const FilteredMultipedigree&  my_mped;
    const instructions&  my_instructions;
    
    std::map<member, pool_pheno_seq>  my_member_directory;
    std::map<pool_pheno_seq, size_t>  my_phenotype_directory;
};


#include "decipher/pool.ipp"

}
} 

#endif

