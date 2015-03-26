#ifndef DECIPHER_FOUNDERS_H
#define DECIPHER_FOUNDERS_H

//============================================================================
// File:      founders.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/11/5 - created.                                   djb
//                                                                          
// Notes:     Classes for generating and storing "phenotype" a set of founders
//            for each subpedigree as used in the EM algorithm.  See phenotype
//            class in em.h.
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "boost/smart_ptr.hpp"
#include "mped/mp_utilities.h"
#include "fped/fped_func.h"
#include "rped/genome_description.h"
#include "mcmc/founder_allele_graph.h"
#include "gelim/ped_imodel_gen.h"
#include "decipher/family.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace MCMC;
using namespace MLOCUS;
using namespace FPED;

typedef std::vector<size_t> allele_group;   // Ordered group of alleles.  Analogous to phased genotype.
typedef pair<MLOCUS::inheritance_model, set<allele_group> >  ags;  // (allele_group_set) Allele groups 
                                                                   // at a given locus. 
typedef vector<ags>  ags_seq;

ostream& operator <<(ostream& out, const ags_seq& seq);
void  write_allele_group(ostream& out, MLOCUS::inheritance_model* locus, allele_group ag);

//----------------------------------------------------------------------------
//  Class:    founders_generator
//                                                                         
//  Purpose:  Provide iteration through sets of haplotype combinations from
//            families in the pedigree data file.
//                                                                          
//----------------------------------------------------------------------------
//
class founders_generator
{
  public:
    typedef FilteredMultipedigree::pedigree_const_pointer  pedigree;
    typedef FilteredMultipedigree::member_const_iterator  member_iterator;

    class iterator
    {
      public:
        friend class founders_generator;
      
        bool  operator ==(const iterator& other) const;
        bool  operator !=(const iterator& other) const;
        
        member  get_member() const;
        pedigree  get_pedigree() const;
        const vector<size_t>&  get_inconsistent_loci() const;
        
        iterator  operator  ++();
        set<hap_seq_comb>&  operator  *();
        set<hap_seq_comb>*  operator ->();
      
      private:
        iterator(founders_generator* generator, bool end = false);

        void  dump_phenotypes(ostream& out) const;
        void  memory_error();
        void  recombination_error();
        bool  build_combs();
          void  set_hap_count();
          void  init_ags_seq();
          void  reset_ags_seq(const vector<boost::shared_ptr<FounderAlleleGraph> >& graphs);
          void  build_combinations(bool& messages_interrupted); // throw(bad_alloc);
      
        // Data members.
        ags_seq  my_ags_seq;
        set<hap_seq_comb>  my_hap_seq_combs;
        vector<member>::const_iterator  my_member_iter;
        vector<member>::const_iterator  my_member_end_iter;
        size_t  my_hap_count;
        founders_generator*  my_generator;
        
        FPED::FilteredMultipedigree::subpedigree_const_pointer  my_subped_ptr;
        vector<size_t>  my_inconsistent_loci;           // indices into marker info
        vector<size_t>  my_inconsistent_group_loci;     // indices into locus group
    };
    
    // Constructor/destructor.
    founders_generator(const vector<member>& members, const locus_group& loci,
                      cerrorstream& errors, ostream& messages, output_state& ostate);
    
    cerrorstream&  errors();
    ostream&  messages();
    bool&  messages_interrupted();
    const locus_group&  loci() const;
    const vector<member>&  members() const;
    
    // Iteration.
    iterator  begin();
    iterator  end();
    
  private:
  
    // Data members.
    const vector<member>&  my_members;
    const locus_group&  my_loci;
    cerrorstream&  my_errors;
    ostream&  my_messages;
    output_state&  my_output_state;
};


//----------------------------------------------------------------------------
//  Class:    founders_comb_generator
//                                                                         
//  Purpose:  Generate combinations of haplotypes consistent w. sequences of
//            sets of ordered allele id groups.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class founders_comb_generator
{
  public:
  
    // - Set of ordered allele groups.  Analogous to comb_generator::geno_set
    //   (see unrelated.h).
    // 
    typedef vector<vector<size_t> >  group_set;
    
    // - Sequence of group sets. Analogous to comb_generator::geno_seq
    //   (see unrelated.h).
    //
    typedef vector<group_set>  group_seq;
    
    typedef MLOCUS::allele_iterator  allele_iterator;
    
    static void  generate(const ags_seq& a_seq, set<hap_seq_comb>& combs, size_t hap_count);
    
    static void  write_group_set(ostream& out, const group_set& grp);
    static void  write_group_seq(ostream& out, const group_seq& seq);
    
  private:
    static void  init_group_seq(const ags_seq& a_seq, group_seq& grp_seq);
      static void  expand_unknown_alleles(const ags& grp_set, set<allele_group>& working_set);
        static void  expand_group(const allele_group& ag,
                                  const MLOCUS::inheritance_model& locus,
                                  set<allele_group>& working_set   );
          static void  populate_new_group(allele_group& new_group,
                                          allele_group::const_iterator old_begin_iter,
                                          allele_group::const_iterator old_end_iter,
                                          allele_group::const_iterator npos_position,
                                          allele_iterator allele_iter                 );
    
    
    static void  generate(const ags_seq& a_seq, group_seq& grp_seq, set<hap_seq_comb>& combs, size_t hap_count);
    static hap_seq_comb  comb(const ags_seq& seq, const group_seq& grp_seq, size_t hap_count);
    
    /* For forcing bad_alloc() exception to test exception handling.
    static size_t  chunk_size;
    */
};

ostream& operator <<(ostream& out, founders_comb_generator::group_set grp);


//----------------------------------------------------------------------------
//  Class:    founders_em_phenotype_map
//                                                                         
//  Purpose:  Container for phenotypes of founders from a families 
//            in the EM algorithm.
//
//----------------------------------------------------------------------------
//
class founders_em_phenotype_map : public family_em_phenotype_map
{
  public:
  
    // Constructor/destructor.
    founders_em_phenotype_map(output_state& ostate,
                              APP::Output_Streams& streams,
                              const FilteredMultipedigree& filtered_mped,
                              const pair<vector<member>, vector<member> >& members, 
                              const instructions& instr, 
                              const locus_group& loci,
                              const string& inner_sub_ped_name = "",
                              const string& outer_sub_ped_name = "");

  private:
  
    void  build();
      void  build_filtered_member_pools();
      void  dump_filtered_member_pools(ostream& out) const;
      void  founders_sub_build(founders_generator::iterator iter, 
                             founders_generator::iterator end_iter);                                   

    // Data members.
    const vector<member>&  my_alt_member_pool;       // founder pool families    
    vector<member>  my_alt_filtered_member_pool;
};

#include "decipher/founders.ipp"

}
} 

#endif

