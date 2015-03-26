#ifndef DECIPHER_PARTITIONING_H
#define DECIPHER_PARTITIONING_H

//============================================================================
// File:      partitioning.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   4/18/7 - created.                                   djb
//                                                                          
// Notes:     declaration of a partitioner class. 
//                                                                          
// Copyright (c) 2007 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/founders.h"

namespace SAGE
{

namespace DECIPHER
{

//----------------------------------------------------------------------------
//  Class:    partitioner
//                                                                          
//  Purpose:  Divide pedigree data file records into sub-populations per 
//            user instructions.
//                                                                          
//----------------------------------------------------------------------------
//
class partitioner
{
  public:
  
    typedef instructions::partition_data  partition_data;
    typedef vector<partition_data>  partition_vector;
  
    // - map<outer partition name, map<inner partition name, vector<member (rep or singleton)> > >
    //
    typedef std::map<string, std::map<string, vector<member> > >  sub_pop_directory;
    
    // - map<outer partition name, map<inner partition name, 
    //      pair<vector<member (analyze subpedigree founder pool)>, vector<member (rep or singleton)> > > >
    //
    typedef std::map<string, std::map<string, pair<vector<member>, vector<member> > > >  founders_sub_pop_directory;    
  
    partitioner(const instructions& instr,
                const FilteredMultipedigree& filtered_mped, cerrorstream& errors);
    
    void  partition_population(sub_pop_directory& members, 
                               const locus_group& loci, const partition_vector& partitions) const;
                                   
    void  partition_population(founders_sub_pop_directory& members, 
                               const locus_group& loci, const partition_vector& partitions) const;

  private:

    //   M = member container type
    //
    template<typename M> 
    void  init_sub_pop_directory(std::map<string, std::map<string, M> >& sub_pop_members,
                                 const partition_vector& partitions) const;
                                 
    void  find_singletons(sub_pop_directory& sub_pop_members, 
                          pedigree_iterator& ped, const locus_group& loci,
                          const partition_vector& partitions) const;
                          
    void  find_singletons(founders_sub_pop_directory& sub_pop_members,
                          pedigree_iterator& ped, const locus_group& loci,
                          const partition_vector& partitions) const; 

    pair<bool, pair<string, string> >  
    partition_data_consistent(subpedigree_iterator& sped,
                              const partition_vector& partitions) const;

    member  get_informative_founder(subpedigree_iterator& sped, const locus_group& loci) const;    

    bool      designated_by_user(member_iterator m_iter) const;
    size_t    get_genotyped_count(member_iterator m_iter, const locus_group& loci) const;

    void      get_genotype_cnt_sorted_members(std::map<size_t, vector<member_iterator>, greater<size_t> >& g_members,
                                              subpedigree_iterator sped, const locus_group& loci) const;

    pair<member_iterator, pair<string, string> > 
    get_best_family_rep(subpedigree_iterator sped, 
                        const locus_group& loci,
                        const partition_vector& partitions) const;
    
    bool  member_missing_partition_value(const member_iterator& m_iter, 
                                         const partition_data& partition) const;
    
    pair< bool, pair<string, string> > member_belongs(member_iterator m_iter, 
                                                      const partition_vector& partitions) const;
    pair<bool, string> member_belongs_to_partition(member_iterator m_iter, 
                                                   const partition_data& partition) const;

    pair<bool, string> check_double_partition(double trait_value, const std::map<string, instructions::value>& sub_pops) const;
    pair<bool, string> check_string_partition(string string_value, const std::map<string, instructions::value>& sub_pops) const;

    void  dump_sub_pops_map(const sub_pop_directory& sub_pops) const;
    void  dump_sub_pops_map(const founders_sub_pop_directory& sub_pops) const;    

    
    // Data members
    const instructions&  my_instructions;
    const FilteredMultipedigree&  my_filtered_mped;
    cerrorstream&  my_errors;
};

#include "decipher/partitioning.ipp"

}
}

#endif

