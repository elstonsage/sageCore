#ifndef DECIPHER_OUTPUT_H
#define DECIPHER_OUTPUT_H

//============================================================================
// File:      output.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   11/15/4 - created.                                   djb
//                                                                          
// Notes:     Output writing classes.
//
//
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/family.h"

namespace SAGE
{

namespace DECIPHER
{

class order_by_sub_pop
{
  public:
    bool operator ()(const base_em_phenotype_map* map1, 
                     const base_em_phenotype_map* map2 );
};


//----------------------------------------------------------------------------
//  Class:    writer
//                                                                         
//  Purpose:  creates program output.
//
//  Note:     templatized on type of em_phenotype_map.
//            
//----------------------------------------------------------------------------
//
template<typename T> class writer
{
  public:

    // Constructor/destructor.
    writer(const output_state& os, ostream& summary, ostream& detail, const instructions& instr,
           const vector<p_values>& pvs, const vector<ld_record>& lds);
    writer(const output_state& os, ostream& summary, ostream& detail, const instructions& instr, 
           const vector<p_values>& pvs, const vector<const T*>& phenotype_maps, const vector<ld_record>& lds);
    ~writer();
    
    void  init_output();    
    void  write(const locus_group& loci, const string& region = "", size_t block_number = 0);    
    void  write(const string& region = "", size_t block_number = 0);
    
    static void  write_sub_pop_name(ostream& out, const string& sub_pop_name);
    static void  write_blank_lines(ostream& out, size_t count);    
    
  private:

    void  set_precision();

    void  write_title(ostream& out);
    void  write_options();
    void  write_prelims(ostream& out, const string& region, size_t block_number);    
    void  write_markers(ostream& out) const;
    void  write_ld() const;
      size_t  find_max_marker_width() const;
    
    // Population haplotype.
    //
    void  write_freq_est(const base_em_phenotype_map& phenotypes, const string& sub_pop_name = "");  
    void  write_all_possible_haplotypes(const base_em_phenotype_map& phenotypes, const string& sub_pop_name = "");
      
    // Most likely diplotypes for individuals.  
    //
    struct ml_widths
    {
      ml_widths();
    
      size_t  space;
      size_t  ped;
      size_t  mem;
      size_t  dip;        
    };    
    
    void write_most_likely();
    void  write_most_likely_pop(const member_em_phenotype_map& phenotypes);
    ml_widths  write_most_likely_header(const member_em_phenotype_map& phenotypes);
    
      void calc_max_widths(const member_em_phenotype_map& phenotypes, ml_widths& widths);
      
    // Write possible diplotypes for individuals.
    //
    enum  format { TABULAR, NON_TABULAR };    
    void  write_possible();
      void  write_possible_dispatch(format f);
      void  write_synteny_note();
        void  write_possible_long_pop(const member_em_phenotype_map* phenotypes);
          size_t  diplotype_col_width(const std::map<string, std::map<member, bool> >& table) const;
          size_t  member_col_width(const set<member, member_order<member> >& members) const;
          void  write_possible_table_header(const set<member, member_order<member> >& members, 
                                                                size_t offset, size_t col_width);
                                                                
                                                                
                                                                
        /*    Does not work with missing data.  NOT USED.
        
        void  write_possible_short_pop(const member_em_phenotype_map* phenotypes);
          size_t  max_member_width(const member_em_phenotype_map* phenotypes);
            string  member_name(member m);
          void  write_member_possibilities(const em_phenotype& phenotype, 
                                           const em_haplotype_map& haplotypes,
                                           const locus_group& loci,
                                           size_t member_width);
          
          static void  build_new_combs(const em_phenotype& phenotype,
                                       const em_haplotype_map& haplotypes,
                                       const locus_group& loci,
                                       vector<pair<vector<size_t>, vector<size_t> > >& new_combs);
            static bool  multiple_allele_pairs(size_t l, 
                                               const vector<pair<vector<size_t>, vector<size_t> > >& new_combs);
              static bool  same_alleles(const pair<size_t, size_t>& p1, const pair<size_t, size_t>& p2);
            static void  insert_missing_index(size_t l, vector<pair<vector<size_t>, vector<size_t> > >& new_combs);
                                       
                                       
          size_t  calc_indicator(const vector<size_t>& indicators,
                                 const vector<pair<vector<size_t>, vector<size_t> > >& new_combs);
                                        
          struct two_locus_diplotype
          {
            two_locus_diplotype(size_t l1, size_t l2, size_t h1l1, size_t h1l2, size_t h2l1, size_t h2l2);
            
            bool  operator ==(const two_locus_diplotype& other);
            bool  operator !=(const two_locus_diplotype& other);
          
            // Loci.
            size_t  locus1;
            size_t  locus2;
            
            // Alleles.
            size_t  hap1_locus1;
            size_t  hap1_locus2;
            size_t  hap2_locus1;
            size_t  hap2_locus2;
          };
          
          static void  cull_phase_unknown(vector<string>& indicators);
          static void  write_two_locus_diplotype(ostream& out, const two_locus_diplotype& dip,
                                                   const locus_group& loci);
          void  write_possible_lines(const vector<vector<string> >& lines, size_t member_width);
            void  calc_column_widths(const vector<vector<string> >& lines, vector<size_t>& col_widths);
        */
        
        
            
    // - Likelihood ratio and permutation tests.
    //
    void  write_p_values();        

    // Data members.
    const output_state&  my_output_state;
    ostream&  my_summary;
    ostream&  my_detail;
    const instructions&  my_instructions;
    const vector<p_values>&  my_p_values;
    vector<const T*>  my_phenotype_maps;
    const vector<ld_record>&  my_lds;
    
    ios::fmtflags  orig_summary_fmt;
    ios::fmtflags  orig_detail_fmt;
};


#include "output.ipp"

}
} 

#endif

