#ifndef DECIPHER_UNRELATED_H
#define DECIPHER_UNRELATED_H

//============================================================================
// File:      unrelated.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/26/4 - created.                                   djb
//                                                                          
// Notes:     Classes for generating and storing phenotypes of unrelated
//            individuals as used in the EM algorithm.
//
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "rped/genome_description.h"
#include "decipher/em.h"
#include "decipher/instructions.h"


namespace SAGE
{

namespace DECIPHER
{

typedef pair<const penetrance_model*, size_t>  pheno;     // locus, phenotype id
typedef vector<pheno>  pheno_seq;

ostream&  operator <<(ostream& out, const pheno_seq& seq);

//----------------------------------------------------------------------------
//  Class:    pheno_seq_generator
//                                                                         
//  Purpose:  Provide iteration through phenotype sequences in population
//            data.
//                                                                          
//----------------------------------------------------------------------------
//
class pheno_seq_generator
{
  public:
  
    class iterator
    {
      public:
        friend class pheno_seq_generator;
      
        bool  operator ==(const iterator& other) const;
        bool  operator !=(const iterator& other) const;
        
        member      get_member() const;
        
        iterator    operator ++();
        pheno_seq&  operator *();
        pheno_seq*  operator ->();
      
      private:
        iterator(const pheno_seq_generator* generator, bool end = false);
        
        void  build_seq();
      
        // Data members.
        pheno_seq  my_pheno_seq;
        vector<member>::const_iterator  my_member_iter;
        vector<member>::const_iterator  my_member_end_iter;
        
        const pheno_seq_generator*  my_generator;        
    };
    
    // Constructor/destructor.
    pheno_seq_generator(const vector<member>& members, const locus_group& loci);
    
    const vector<member>&  members() const;
    const locus_group&  loci() const;
    
    // Iteration.
    iterator  begin() const;
    iterator  end() const;
    
  private:
  
    // Data members.
    const vector<member>& my_members;
    const locus_group& my_loci;
};


//----------------------------------------------------------------------------
//  Class:    comb_generator
//                                                                         
//  Purpose:  Generate combinations of haplotypes consistent w. phenotype
//            sequences.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class comb_generator
{
  public:
    typedef penetrance_model::phased_penetrance_iterator  pen_iter;
    
    typedef vector<int>  geno_set;      // - Genotypes at a given locus 
                                        //   consistent w. a given phenotype.
    typedef vector<geno_set>  geno_seq;
    
    typedef vector<size_t>        x_allele_set;
    typedef vector<x_allele_set>  x_allele_seq;  
  
    static void  generate(pheno_seq_generator::iterator& p_seq_iter, set<hap_seq_comb>& combs);
    
  private:
    static bool  x_linked(pheno_seq_generator::iterator& p_seq_iter);  
    static bool  male(pheno_seq_generator::iterator& p_seq_iter);
    static bool  female(pheno_seq_generator::iterator& p_seq_iter);    
    static void  generate_x_linked_male(pheno_seq_generator::iterator& p_seq_iter, set<hap_seq_comb>& combs);
    static void  generate_x_linked_male(x_allele_seq& a_seq, set<hap_seq_comb>& combs);
    static void  init_geno_seq(const pheno_seq& p_seq, geno_seq& g_seq);
    static void  generate(const pheno_seq& p_seq,
                          geno_seq& g_seq, set<hap_seq_comb>& combs);
    static hap_seq_comb  comb(const pheno_seq& p_seq, const geno_seq& g_seq);
};


//----------------------------------------------------------------------------
//  Class:    unrelated_em_phenotype_map
//                                                                         
//  Purpose:  Container for phenotypes of unrelated individuals in the EM 
//            algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class unrelated_em_phenotype_map : public member_em_phenotype_map
{
  public:
  
    // Constructor/destructor.
    unrelated_em_phenotype_map(output_state& ostate,
                     APP::Output_Streams& streams,
                     const FilteredMultipedigree& mped,
                     const vector<member>& members, 
                     const instructions& instr, 
                     const locus_group& loci, 
                     const string& inner_sub_pop_name = "",
                     const string& outer_sub_pop_name = "");

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
    
    std::map<member, pheno_seq>  my_member_directory;
    std::map<pheno_seq, size_t>  my_phenotype_directory;
};

#include "decipher/unrelated.ipp"

}
} 

#endif

