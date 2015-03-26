#ifndef DECIPHER_FAMILY_H
#define DECIPHER_FAMILY_H

//============================================================================
// File:      family.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   6/3/4 - created.                                   djb
//                                                                          
// Notes:     Classes for generating and storing "phenotype" of one individual
//            per subpedigree as used in the EM algorithm.  See phenotype
//            class in em.h.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <algorithm>
#include "boost/smart_ptr.hpp"
#include "mped/mp_utilities.h"
#include "fped/fped_func.h"
#include "mcmc/founder_allele_graph.h"
#include "gelim/ped_imodel_gen.h"
#include "decipher/unrelated.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace MCMC;
using namespace MLOCUS;
using namespace FPED;

const boost::shared_ptr<FounderAlleleGraph>  SKIP_MARKER(static_cast<FounderAlleleGraph*>(0));

typedef std::pair<size_t, size_t> allele_pair;
typedef pair<MLOCUS::inheritance_model, set<allele_pair> >  aps;  // (allele_pair_set) Allele pairs at a given locus 
typedef vector<aps>  aps_seq;

ostream& operator <<(ostream& out, const aps_seq& seq);
void  write_allele_pair(ostream& out, MLOCUS::inheritance_model* locus, allele_pair ap);

//----------------------------------------------------------------------------
//  Class:    family_generator
//                                                                         
//  Purpose:  Provide iteration through sets of haplotype combinations from
//            families in the pedigree data file.
//                                                                          
//----------------------------------------------------------------------------
//
class family_generator
{
  public:
    typedef FilteredMultipedigree::pedigree_const_pointer  pedigree;
    typedef FilteredMultipedigree::member_const_iterator  member_iterator;

    class iterator
    {
      public:
        friend class family_generator;
      
        bool  operator ==(const iterator& other) const;
        bool  operator !=(const iterator& other) const;
        
        member  get_member() const;
        pedigree  get_pedigree() const;
        const vector<size_t>&  get_inconsistent_loci() const;
        
        iterator  operator ++();
        set<hap_seq_comb>&  operator  *();
        set<hap_seq_comb>*  operator ->();
      
      private:
        iterator(family_generator* generator, bool end = false);

        void  dump_phenotypes(ostream& out) const;
        void  memory_error();
        void  recombination_error();
        bool  build_combs();
          void  init_aps_seq();
          void  reset_aps_seq();
          void  build_combinations(bool& messages_interrupted) throw(bad_alloc);
      
        // Data members.
        aps_seq  my_aps_seq;
        aps_seq  my_initialized_aps_seq;
        set<hap_seq_comb>  my_hap_seq_combs;
        vector<member>::const_iterator  my_member_iter;
        vector<member>::const_iterator  my_member_end_iter;
        family_generator*  my_generator;
        
        FPED::FilteredMultipedigree::subpedigree_const_pointer  my_subped_ptr;
        vector<size_t>  my_inconsistent_loci;        // indices into marker info
        vector<size_t>  my_inconsistent_group_loci;  // indices into locus_group
    };
    
    // Constructor/destructor.
    family_generator(const vector<member>& members, const locus_group& loci,
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
//  Class:    family_comb_generator
//                                                                         
//  Purpose:  Generate combinations of haplotypes consistent w. sequences of
//            groups of ordered allele id pairs.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class family_comb_generator
{
  public:
  
    // - Group of ordered allele pairs.  Analogous to comb_generator::geno_set
    //   (see unrelated.h).
    // 
    typedef vector<pair<size_t, size_t > >  pair_set;
    
    // - Sequence of pair sets. Analogous to comb_generator::geno_seq
    //   (see unrelated.h).
    //
    typedef vector<pair_set>  pair_seq;
    
    typedef MLOCUS::allele_iterator  allele_iterator;
    
    // - Characterizes allele pair w. regard to unknown alleles.
    //
    enum  UNKNOWNS  { none, first_allele, second_allele, both };
  
    static void  generate(const aps_seq& a_seq, set<hap_seq_comb>& combs);
    
    static void  write_pair_set(ostream& out, const pair_set& ps);
    static void  write_pair_seq(ostream& out, const pair_seq& ps);
    
  private:
    static void  init_pair_seq(const aps_seq& a_seq, pair_seq& pr_seq);
    static void  expand_unknown_alleles(pair_set& pr_set, const penetrance_model& locus);
    static UNKNOWNS  unknowns(const allele_pair& pair);
    static bool  is_new_pair(const pair_set& pr_set, const allele_pair& pr);
    static void  generate(const aps_seq& a_seq, pair_seq& pr_seq, set<hap_seq_comb>& combs);
    static hap_seq_comb  comb(const aps_seq& seq, const pair_seq& pr_seq);
};

ostream& operator <<(ostream& out, family_comb_generator::pair_set ps);

class family_em_phenotype_map;


//----------------------------------------------------------------------------
//  Class:    not_uninformative_leaf
//                                                                         
//  Purpose:  functor to be used with MPFilter::filter_multipedigree 
//            (see fped/fped_func.h).  Will be used to eliminate persons
//            who have no informative markers in the haplotype range and who
//            are not needed in the descent graphs (i. e. uniformative
//            subpedigree "leaves").  7-24-6:  modified to leave one child
//            regardless of informativity.
//                                                                          
//----------------------------------------------------------------------------
//
class not_uninformative_leaf : public unary_function<FilteredMultipedigree::member_type, bool>
{
  public:
  
    not_uninformative_leaf(FPED::has_informative_loci<FilteredMultipedigree::member_type>& base_filter);
    
    bool  operator()(const FilteredMultipedigree::member_type& member) const;
    
  private:
    
    FPED::has_informative_loci<FilteredMultipedigree::member_type>  genotyped;
      static bool  keeper(const FilteredMultipedigree::member_type& member);
};


//----------------------------------------------------------------------------
//  Class:    family_em_phenotype_map
//                                                                         
//  Purpose:  Container for phenotypes of selected individuals from a family 
//            in the EM algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class family_em_phenotype_map : public member_em_phenotype_map
{
  public:
  
    // Constructor/destructor.
    family_em_phenotype_map(output_state& ostate,
                            APP::Output_Streams& streams,
                            const FilteredMultipedigree& filtered_mped,
                            const vector<member>& members, 
                            const instructions& instr, 
                            const locus_group& loci,
                            const string& inner_sub_ped_name = "",
                            const string& outer_sub_ped_name = "" ,
                            bool build_phenotypes = true);

    const em_phenotype&  operator [](member m) const;
    set<member, member_order<member> > members() const;
    const em_phenotype&  final_phenotype(member m) const;    
    
    const FilteredMultipedigree&  multipedigree() const;    
    
    void  dump(ostream& out) const;
    void  dump_filtered_mped(ostream& out) const;    
    
  protected:
  
    void init_filter(FPED::has_informative_loci<FilteredMultipedigree::member_type>& base_filter) const;
    void  build();
      void  build_filtered_member_pool();
      void  dump_filtered_member_pool(ostream& out) const;
      void  divide_member_pool(vector<member>& family_members, vector<member>& singletons) const;
      
      template<typename T> 
      void  related_sub_build(T iter, T end_iter);            
      
      void  singleton_sub_build(pheno_seq_generator::iterator iter, 
                                pheno_seq_generator::iterator end_iter, 
                                comb_generator& comb_gen);

    size_t  pheno_index(member m) const;

    // Data members.
    FilteredMultipedigree  my_filtered_mped;
    const instructions&  my_instructions;
    
    vector<member>  my_filtered_member_pool;
    std::map<member, size_t>  my_phenotype_directory;
};

#include "decipher/family.ipp"

}
} 

#endif

