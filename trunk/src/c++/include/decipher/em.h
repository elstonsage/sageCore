#ifndef DECIPHER_EM_H
#define DECIPHER_EM_H

//============================================================================
// File:      em.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   2/26/4 - created.                                   djb
//                                                                          
// Notes:     Classes representing haplotypes and phenotypes as used  
//            in the EM algorithm.
//
//
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "app/output_streams.h"
#include "rped/rped.h"
#include "numerics/log_double.h"
#include "numerics/cephes.h"
#include "decipher/definitions.h"


namespace SAGE
{

namespace DECIPHER
{

typedef FPED::FilteredMultipedigree  FilteredMultipedigree;

using namespace RPED;
using namespace MLOCUS;

typedef vector<size_t>  hap_seq;    // allele id's for a sequence of markers
typedef multiset<hap_seq>  hap_seq_comb;

std::string  hap_seq_string(const hap_seq& seq, const locus_group& loci);
void  write_hap_seq(ostream& out, const hap_seq& seq, const locus_group& loci);
void  write_markers(ostream& out, const locus_group& loci);

void  dump_locus(ostream& out, const inheritance_model* locus);

class em_phenotype;
class base_em_phenotype_map;

template <typename T>
class member_order
{
  public:
    bool  operator ()(T m1, T m2) const;
};

//----------------------------------------------------------------------------
//  Class:    em_haplotype
//                                                                         
//  Purpose:  Represent a haplotype in the EM algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class em_haplotype
{
  public:
    
    // Constructor/destructor.
    // Note: using default copy constructor with hap_seq* data member.  This should not
    //       pose a problem as haplotypes and haplotype_maps are used only in the context of a single
    //       em_phenotype_map and haplotype sequences don't change once an em_phenotype map
    //       is constructed.
    // 
    em_haplotype();
    em_haplotype(const hap_seq* seq);
    
    double  old_freq(size_t total_hap_count) const;
    double  new_freq(size_t total_hap_count) const;
    const hap_seq&  sequence() const;
    
    void  incr(double amount);
    void  reset();
    void  zero();
    bool  converged(size_t total_hap_count, double epsilon) const;
    bool  similar(const em_haplotype& other, size_t total_hap_count, double epsilon) const;    
    
  private:
    
    // Data members.
    const hap_seq*  my_sequence;
    
    double  my_old_count;
    double  my_new_count;
    double  my_static_count;
};

//----------------------------------------------------------------------------
//  Class:    hap_freq
//                                                                         
//  Purpose:  Represent a haplotype frequency for purpose of displaying
//            results.
//            
//----------------------------------------------------------------------------
//
struct hap_freq
{
  hap_freq(const string& haplotype, double frequency);
  bool  operator >(const hap_freq& other) const;

  const string  hap;
  const double  freq;
};

//----------------------------------------------------------------------------
//  Class:    em_haplotype_map
//                                                                         
//  Purpose:  Container for haplotypes in the EM algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class em_haplotype_map
{
  public:
    
    // Constructor/destructor.
    em_haplotype_map();
    
    bool  contains(const hap_seq& key) const;
    bool  converged(size_t total_hap_count, double epsilon) const;
    bool  similar(em_haplotype_map& other, size_t total_hap_count, double epsilon) const;    
    
    double  total_freq(size_t total_hap_count) const;
    
    em_haplotype&   operator [](size_t index);
    const  em_haplotype&  get_haplotype(size_t index) const;
    const hap_seq&  index_to_hap_seq(size_t index) const;
    size_t          hap_seq_to_index(const hap_seq& key) const;
    
    const std::map<size_t, em_haplotype>&  haplotypes() const;
    
    size_t  add_haplotype(const hap_seq& key);    
    void  reset();       // Called at ea. iteration of algorithm.
    void  zero();        // Called at beginning of algorithm.
    void  update_counts(const base_em_phenotype_map& pm);

    void  dump(ostream& out, size_t total_hap_count, const locus_group& loci) const;
  
  private:
    
    // Data members.
    mutable std::map<hap_seq, size_t>  my_hap_seqs;    
    mutable std::map<size_t, em_haplotype>  my_haplotypes;
    
    bool  counts_initialized;
};

//----------------------------------------------------------------------------
//  Class:    comb_prob
//                                                                         
//  Purpose:  Represent a haplotype combination probability for purpose of 
//            displaying results.
//            
//----------------------------------------------------------------------------
//
struct comb_prob
{
  comb_prob(string c, double p);
  bool  operator >(const comb_prob& other) const;

  string  comb;
  double  prob;
};


//----------------------------------------------------------------------------
//  Class:    em_phenotype
//                                                                         
//  Purpose:  Represent a phenotype in the EM algorithm.
//            
//----------------------------------------------------------------------------
//
class em_phenotype
{
  public:
    typedef vector<size_t>  combination;
    
    // Constructor/destructor.
    em_phenotype();
    em_phenotype(em_haplotype_map* haplotypes, const set<hap_seq_comb>& combinations);
    em_phenotype(em_haplotype_map* haplotypes, const vector<vector<hap_seq> >& combinations);

    const em_haplotype_map&  haplotypes() const;    
    size_t  count() const;
    bool  is_ambiguous() const;
    const vector<pair<combination, double> >&  combinations() const;
    size_t  comb_size() const;
    const set<comb_prob, greater<comb_prob> >&  probabilities(const locus_group& loci,
                                                                    const em_haplotype_map& haplotypes) const;
    const vector<string>&  comb_strs(const locus_group& loci, const em_haplotype_map& haplotypes) const;
    
    void  incr();
    void  init_weights();
    void  update_weights(size_t chromosome_count, em_haplotype_map& haplotypes);
    
    double  calc_prior(const combination& comb, size_t chromosome_count,
                       em_haplotype_map& haplotypes) const;
    
    void  dump(ostream& out, const locus_group& loci, const em_haplotype_map& haplotypes) const;
    
  private:
    
    string  combination_string(const vector<pair<combination, double> >::const_iterator& c_iter,
                               const locus_group& loci, const em_haplotype_map& haps) const;
    
    // - For iterating through distinct haplotypes in a combination
    //   and determining the number of copies.
    //
    class combination_iterator
    {
      public:
        combination_iterator(const combination& comb);
        
        pair<size_t, size_t>  operator *();
        std::map<size_t, size_t>::iterator  operator ++(int);
        bool  at_end() const;  
      
      private:
        bool  previously_found(vector<size_t>::const_iterator iter) const;
      
        std::map<size_t, size_t>  my_distinct_haplotypes;
        std::map<size_t, size_t>::iterator  my_internal_iterator;
    };
    
    // Data members.
    size_t  my_count;
    bool  ambiguous;      // Has more than one combination.
    vector<pair<combination, double> >  my_combinations;     // <group of haplotypes, weight>
    mutable set<comb_prob, greater<comb_prob> >  my_probabilities;
    mutable vector<string>  my_comb_strs;
};

class pop_freq_writer;

//----------------------------------------------------------------------------
//  Class:    base_em_phenotype_map
//                                                                         
//  Purpose:  base container for phenotypes in the EM algorithm.
//            
//                                                                          
//----------------------------------------------------------------------------
//
class base_em_phenotype_map
{
  public:

    typedef vector<em_phenotype>::const_iterator  const_iterator;

    base_em_phenotype_map(output_state& ostate, APP::Output_Streams&  streams, const locus_group& loci, 
                          const string& inner_sub_pop_name = "", 
                          const string& outer_sub_pop_name = "");
    virtual ~base_em_phenotype_map();
    
    void  maximize(double epsilon, size_t total_runs, ostream& dump_file,
                   bool dump = false, double cutoff = 1.0, bool silent = true);    
    bool  is_maximized() const;    


    // These don't change after constructions.
    bool  empty() const;
    const string&  inner_sub_pop_name() const;
    const string&  outer_sub_pop_name() const;
    const string&  sub_pop_name() const;    
    const locus_group&  loci() const;
    size_t  total_hap_count() const;  
    size_t  independent_param_count() const;    


    // These return current state.
    const em_haplotype_map&  haplotypes() const;
    const_iterator  begin() const;
    const_iterator  end() const; 
    const set<hap_freq, greater<hap_freq> >&  current_frequencies() const;    
    virtual void  dump(ostream& out) const;
    

    // These return final results.    
    const em_haplotype_map&  final_haplotypes() const;
    const_iterator  final_begin() const;
    const_iterator  final_end() const;     
    const set<hap_freq, greater<hap_freq> >&  final_frequencies() const;
    double  max_ln_likelihood() const;

  
  protected:
    void  init_weights();      // Randomly assign weights.    
    void  update_weights();    // Called at ea. iteration of algorithm.  
    double  ln_likelihood();    
    void  build_sub_pop_name();
    void  set_frequencies(set<hap_freq, greater<hap_freq> >& frequencies,
                                  const std::map<size_t, em_haplotype>& haps) const;
    void  init_dump(ostream& dump_file) const;
    
    // Data members.
    output_state&  my_output_state;
    APP::Output_Streams&  my_streams;
    cerrorstream&  my_errors;
    ostream&  my_messages;
   
      // These don't change after construction. 
      const locus_group&  my_loci;
      const string  my_inner_sub_pop_name;
      const string  my_outer_sub_pop_name;
      string  my_sub_pop_name;
      size_t  my_total_hap_count;     // Total number of haplotypes (chromosomes) in the sample.
    
      // Current values.
      em_haplotype_map  my_haplotypes;
      vector<em_phenotype>  my_phenotypes;
      mutable set<hap_freq, greater<hap_freq> >  my_current_frequencies;      
      
      // Final values
      em_haplotype_map  my_best_haplotypes;
      vector<em_phenotype>  my_best_phenotypes;
      mutable set<hap_freq, greater<hap_freq> >  my_final_frequencies;
      double  my_max_ln_likelihood;
    
    bool  maximized;
};


//----------------------------------------------------------------------------
//  Class:    member_em_phenotype_map
//                                                                         
//  Purpose:  base class for em_phenotype_maps utilizing members.
//
//  Note:     derived classes, unrelated_em_phenotype_map and family_em
//            phenotype map, are very similar and could be further
//            consolidated.  This would be at the expense of performance
//            for the unrelated class, however, which does not have to build a new
//            phenotype when members have identical marker values.            
//                                                                          
//----------------------------------------------------------------------------
//
class member_em_phenotype_map : public base_em_phenotype_map
{
  public:
    typedef SAGE::DECIPHER::member  member;
  
    member_em_phenotype_map(output_state& ostate, APP::Output_Streams& streams, const vector<member>& members, const locus_group& loci, 
                            const string& inner_sub_pop_name = "", 
                            const string& outer_sub_pop_name = "");  
    ~member_em_phenotype_map();
    
    // - Returns current phenotype, NOT final phenotype.
    //
    virtual const em_phenotype&  operator [](member m) const = 0;
    
    virtual set<member, member_order<member> >  members() const = 0;
    
    virtual const em_phenotype&  final_phenotype(member m) const = 0;
    
  protected:
    
    // Data members.
    const vector<member>&  my_member_pool;       // Potential members.
};

//----------------------------------------------------------------------------
//  Class:    pop_freq_writer
//                                                                         
//  Purpose:  write population frequency estimates.
//            
//----------------------------------------------------------------------------
//
class pop_freq_writer
{
  public:
    pop_freq_writer(ostream& out, double cutoff);
    
    static size_t  max_seq_width(const set<hap_freq, greater<hap_freq> >& haps);
    static void  cutoff_note(ostream& out);
                    
    void  set_frequencies(const set<hap_freq, greater<hap_freq> >* frequencies);
    void  set_ln_likelihood(double);                    
    void  write(size_t precision = 6, bool show_likelihood = true);
    
  private:
  
    // Data members.
    ostream&  my_out;
    const set<hap_freq, greater<hap_freq> >*  my_frequencies;
    double  my_cutoff;
    double  my_ln_likelihood;
};

#include "decipher/em.ipp"

}
} 

#endif

