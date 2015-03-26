//============================================================================
// File:      founders.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 10/11/5                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/founders.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

/* For forcing bad_alloc() exception to test exception handling.
size_t  founders_comb_generator::chunk_size = 10000;
*/

ostream&
operator <<(ostream& out, const ags_seq& seq)
{
  for(size_t l = 0; l < seq.size(); ++l)
  {
    out << seq[l].first.name() << "  ";
    
    set<allele_group>::const_iterator  ags_iter     = seq[l].second.begin();
    set<allele_group>::const_iterator  ags_end_iter = seq[l].second.end();
    for(; ags_iter != ags_end_iter; ++ags_iter)
    {
      out << "(";
      allele_group::const_iterator  ag_iter     = ags_iter->begin();
      allele_group::const_iterator  ag_end_iter = ags_iter->end();
      for(; ag_iter != ag_end_iter; ++ag_iter)    
      {
    
        out << (*ag_iter == MLOCUS::NPOS ? "no info" : seq[l].first.get_allele(*ag_iter).name());
        out << ", ";
      }
      
      out << ")  ";
    }
    
    out << endl;
  }
  
  return  out;
}

void
write_allele_group(ostream& out, MLOCUS::inheritance_model* locus, allele_group ag)
{
  allele_group::const_iterator  ag_iter     = ag.begin();
  allele_group::const_iterator  ag_end_iter = ag.end();
  for(; ag_iter != ag_end_iter; ++ag_iter)    
  {
    out << locus->get_allele(*ag_iter).name();
    out << ", ";
  }
}


//============================================================================
// IMPLEMENTATION:  founders_generator::iterator
//============================================================================
//
void
founders_generator::iterator::dump_phenotypes(ostream& out) const
{
  FPED::FilteredMultipedigree::member_const_iterator  m_iter     = my_subped_ptr->member_begin();
  FPED::FilteredMultipedigree::member_const_iterator  m_end_iter = my_subped_ptr->member_end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    out << "member " << m_iter->name() << endl;
    const locus_group&  loci = my_generator->loci();
    size_t  locus_count = loci.size();
    for(size_t l = 0; l < locus_count; ++l)
    {
      out << "locus " << loci[l].second->name() << endl;
      out << "genotype " 
          << loci[l].second->get_phenotype(m_iter->info().phenotype(loci[l].first)).name() << endl;
    }
  }
}

// - Build sequence of allele group sets.
//
bool
founders_generator::iterator::build_combs()
{
  if(my_member_iter == my_member_end_iter)
  {
    return  false;
  }
  else
  {
    bool& interrupted = my_generator->messages_interrupted();
    my_hap_seq_combs.clear();
    set_hap_count();  
    init_ags_seq();
    
    //cout << "\nallele group set sequence after initialization:\n" << my_ags_seq << endl;
    
    size_t  marker_count = my_ags_seq.size();

    #if 0
      cout << "member " << (*my_member_iter)->name() << endl;
      cout << "pedigree " << (*my_member_iter)->pedigree()->name() << endl;    
      cout << "allele pair sequence after initialization\n" << my_ags_seq << endl;    
      //dump_phenotypes(cout);
    #endif
    
    // - Create founder allele graphs.
    //
    vector<boost::shared_ptr<FounderAlleleGraph> >  graphs;
    McmcMeiosisMap  mm(*my_subped_ptr);
        
    size_t  bit_count = mm.get_meiosis_count();
    
    graphs.reserve(marker_count);
    for(size_t l = 0; l < marker_count; ++l)
    {
      if(my_ags_seq[l].second.size())   // Mendelian inconsistency was found during genotype elimination.
      {
        const string&  locus_name = my_ags_seq[l].first.name();
        const string&  pedigree_name = (*my_member_iter)->pedigree()->name();
        const string&  member_name = (*my_member_iter)->name();
        
        if(! interrupted)
        {
          interrupted = true;
          my_generator->messages() << endl;
        }
        
        my_generator->errors() << priority(information) << "A Mendelian inconsistency was found at locus, "
                               << locus_name << ", of the constituent pedigree in pedigree, " << pedigree_name
                               << ", containing member, " << member_name << ".  Substituting missing values "
                               << "at this locus for all members of the constituent pedigree ..." << endl;
                               
        graphs.push_back(SKIP_MARKER);
      }
      else
      {
        graphs.push_back(boost::shared_ptr<FounderAlleleGraph>(new FounderAlleleGraph(mm, my_ags_seq[l].first)));
      }
    }
    
    bool  valid_pattern_found = false;
    
    // - Examine each inheritance pattern.
    //
    size_t  pattern_count = 1 << bit_count;
    for(__U32 p = 0; p < pattern_count; ++p)
    {
      //cout << "pattern" << p << endl;
      bool  valid = true;
      bit_field  pattern(bit_count, p);
      size_t  graph_count = graphs.size();      
      for(size_t l = 0; l < graph_count; ++l)
      {
        //cout << "locus" << l << endl;
        if(graphs[l] == SKIP_MARKER)
        {
          continue;
        }
        
        graphs[l]->set_pattern(pattern);
        //graphs[l]->dump_graph(cout);        
        valid &= graphs[l]->is_pattern_valid();
      }
      
      //cout << boolalpha << valid << endl;
      
      if(valid)
      {
        valid_pattern_found = true;
      
        // - Remove data pertaining to previous valid pattern!
        //
        reset_ags_seq(graphs);
        //cout << my_ags_seq << endl;
           
        for(size_t l = 0; l < graphs.size(); ++l)
        {
          //cout << "locus " << l << endl; 
          if(graphs[l] == SKIP_MARKER)
          {
            continue;
          }
        
          for(size_t vas = 0; vas < graphs[l]->get_valid_allele_set_count(); ++vas)
          {
            allele_group  ag;
            member_iterator  m_iter     = my_subped_ptr->member_begin();
            member_iterator  m_end_iter = my_subped_ptr->member_end();
            for(; m_iter != m_end_iter; ++m_iter)
            {
              if(MPED::mp_utilities::is_founder(m_iter))     
              {
                //cout << "member " << m_iter->name() << " in pedigree " << m_iter->pedigree()->name() << endl;              
                std::pair<MLOCUS::allele, MLOCUS::allele> ap = 
                        graphs[l]->get_individual_allele_pattern(vas, *m_iter);

                if(ap.first.is_valid())
                {
                  assert(ap.first.name() != "~remap");
                }
                
                if(ap.second.is_valid())
                {
                  assert(ap.second.name() != "~remap");
                }
            
                //cout << " allele pair ";
                //write_allele_pair(cout, &(my_ags_seq[l].first), ap);
                //cout << endl;
              
                ag.push_back(ap.first.id());
                ag.push_back(ap.second.id());
              }
            }
            
            //write_allele_group(cout, &(my_ags_seq[l].first), ag);
            //cout << endl;
 
            my_ags_seq[l].second.insert(ag);
          }
          
          //cout << endl;
        }
      
        // - Consider haplotype combinations pertaining to one inheritance vector (pattern)
        //   at a time.
        //
        //cout << my_ags_seq << endl;
        build_combinations(interrupted);
      }
    }
    
    return  valid_pattern_found;
  }
}

// - Generate haplotype combinations for corresponding to current allele group
//   set sequence.
//
void
founders_generator::iterator::build_combinations(bool& messages_interrupted) //throw(bad_alloc)
{
  founders_comb_generator  comb_gen;
  
  comb_gen.generate(my_ags_seq, my_hap_seq_combs, my_hap_count);
}

void
founders_generator::iterator::set_hap_count()
{
  my_hap_count = 0;

  member_iterator  m_iter     = my_subped_ptr->member_begin();
  member_iterator  m_end_iter = my_subped_ptr->member_end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    if(MPED::mp_utilities::is_founder(m_iter))    
    {
      // - Founder allele graph does not yet handle x-linked markers.
      //
      my_hap_count += 2;    
    }            
  }
}

// - Check for inconsistent loci and supply subpedigree specific inheritance models.
//
void
founders_generator::iterator::init_ags_seq()
{
  my_inconsistent_loci.clear();
  my_inconsistent_group_loci.clear();

  const locus_group&  loci = my_generator->loci();
  size_t  locus_count = loci.size();
  my_ags_seq.resize(locus_count);
  
  // - Check for Mendelian inconsistencies.
  //
  pedigree_imodel_generator  ig;
  for(size_t l = 0; l < locus_count; ++l)
  {
    // - Check for mendelian inconsistencies.
    //
    ig.set_prior_remap(false);
    ig.set_genotype_elimination(true);
    ig.set_post_remap(false);
    
    my_ags_seq[l].first = ig(*my_subped_ptr, loci[l].first);
    my_ags_seq[l].second.clear();    
    
    if(ig.inconsistent())
    {
      // - Use original global inheritance model.  Will be skipped by build_combs().
      //
      my_ags_seq[l].first = *(loci[l].second);
      
      allele_group  dummy(my_hap_count, MLOCUS::NPOS);
      my_ags_seq[l].second.insert(dummy);
      my_inconsistent_loci.push_back(loci[l].first);
      my_inconsistent_group_loci.push_back(l);      
    }
  }
}

void
founders_generator::iterator::reset_ags_seq(const vector<boost::shared_ptr<FounderAlleleGraph> >& graphs)
{
  ags_seq::iterator  iter     = my_ags_seq.begin();
  ags_seq::iterator  end_iter = my_ags_seq.end();
  for(; iter != end_iter; ++iter)
  {
    iter->second.clear();
  }
  
  // - Restore 'inconsistenty tag'.
  //
  size_t inconsistent_count = my_inconsistent_group_loci.size();
  for(size_t l = 0; l < inconsistent_count; ++l)
  {
    allele_group  dummy(my_hap_count, MLOCUS::NPOS);
    my_ags_seq[my_inconsistent_group_loci[l]].second.insert(dummy);
  }
}

//============================================================================
// IMPLEMENTATION:  founders_comb_generator
//============================================================================
//
void
founders_comb_generator::init_group_seq(const ags_seq& a_seq, group_seq& grp_seq)
{
  size_t  locus_count = a_seq.size();
  grp_seq.assign(locus_count, group_set());
  
  for(size_t l = 0; l < locus_count; ++l)
  {
    set<allele_group>  working_set;
    expand_unknown_alleles(a_seq[l], working_set);  
    grp_seq[l].assign(working_set.begin(), working_set.end());
  }
}

// - Replace allele groups containing MLOCUS::NPOS (no information) at any position(s)
//   with groups containing all possible alleles at that position(s) eliminating any
//   resultant duplicate allele groups.
// 
void
founders_comb_generator::expand_unknown_alleles(const ags& grp_set, set<allele_group>& working_set)
{
  set<allele_group>::const_iterator  grp_iter     = grp_set.second.begin();
  set<allele_group>::const_iterator  grp_end_iter = grp_set.second.end();
  for(; grp_iter != grp_end_iter; ++grp_iter)
  {
    expand_group(*grp_iter, grp_set.first, working_set);
  }
}

void  
founders_comb_generator::expand_group(const allele_group& ag,
                                      const MLOCUS::inheritance_model& locus,
                                      set<allele_group>& working_set   )
{
  allele_group::const_iterator  ag_begin_iter = ag.begin();
  allele_group::const_iterator  ag_end_iter   = ag.end();
  allele_group::const_iterator  ag_iter = find(ag_begin_iter, ag_end_iter, MLOCUS::NPOS);
  if(ag_iter != ag_end_iter)
  {
    allele_iterator  allele_iter     = locus.allele_begin();
    allele_iterator  allele_end_iter = locus.allele_end();
    for(; allele_iter != allele_end_iter; ++allele_iter)
    {
      allele_group  new_group;
      populate_new_group(new_group, ag_begin_iter, ag_end_iter, ag_iter, allele_iter);
      expand_group(new_group, locus, working_set);
    }            
  }
  else
  {
    working_set.insert(ag);
  }
}

// - Populate new group with elements of old group substituting a new allele in the
//   position where allele information is missing.
//
void  
founders_comb_generator::populate_new_group(allele_group& new_group,
                                            allele_group::const_iterator old_begin_iter,
                                            allele_group::const_iterator old_end_iter,
                                            allele_group::const_iterator npos_position,
                                            allele_iterator allele_iter                 )
{
  allele_group::const_iterator  old_iter = old_begin_iter;
  for(; old_iter != old_end_iter; ++old_iter)
  {
    if(old_iter != npos_position)
    {
      new_group.push_back(*old_iter);
    }
    else
    {
      new_group.push_back(allele_iter->id());
    }
  }
}                                      

void
founders_comb_generator::generate(const ags_seq& a_seq, set<hap_seq_comb>& combs, size_t hap_count)
{
  group_seq  grp_seq;
  init_group_seq(a_seq, grp_seq);
  generate(a_seq, grp_seq, combs, hap_count);
}

// - Generate haplotypes recursively by breaking first set of genotypes
//   encountered into a set with one genotype and a set with the remaining
//   genotypes.  The process is repeated until there is one genotype at
//   every locus, ie, the haplotypes are unambiguous.
//
void
founders_comb_generator::generate(const ags_seq& a_seq,
                         group_seq& grp_seq, set<hap_seq_comb>& combs, size_t hap_count)
{
  /* For forcing bad_alloc() exception to test exception handling.
  vector<double>*  leak_ptr = new vector<double>(chunk_size);
  cout << "chunk_size " << chunk_size << " leak_ptr " << leak_ptr << endl;
  */

  for(size_t l = 0; l < grp_seq.size(); ++l)
  {
    while(grp_seq[l].size() > 1)   
    {
      group_seq  new_grp_seq(grp_seq);
      group_set  new_grp_set;
      
      new_grp_set.push_back(grp_seq[l][grp_seq[l].size() - 1]);
      new_grp_seq[l] = new_grp_set;
      
      grp_seq[l].pop_back();
      generate(a_seq, new_grp_seq, combs, hap_count);
    }
  }
  
  // - Note combinations are multsets, hence elements are sorted.
  //   combs is a set and therefore will not allow a duplicate combinations
  //   to be inserted regardless of the haplotype order used to construct them.
  //  
  combs.insert(comb(a_seq, grp_seq, hap_count));
}

// - Create a combination of haplotyes from a sequence of allele group sets
//   containing exactly one allele group at every locus.
//
hap_seq_comb
founders_comb_generator::comb(const ags_seq& a_seq, const group_seq& grp_seq, size_t hap_count)
{
  size_t  locus_count = a_seq.size();
  assert(grp_seq.size() == locus_count);

  vector<hap_seq>  haps;
  haps.assign(hap_count, hap_seq());
  
  for(size_t l = 0; l < locus_count; ++l)
  {
    assert(! grp_seq[l].empty());
    for(size_t h = 0; h < hap_count; ++h)
    haps[h].push_back(grp_seq[l][0][h]);
  }
  
  hap_seq_comb  c;
  
  c.insert(haps.begin(), haps.end());
  
  return  c;
}

void  
founders_comb_generator::write_group_set(ostream& out, const group_set& grp)
{
  size_t  group_count = grp.size();

  for(size_t g = 0; g < group_count; ++g)
  {
    size_t  allele_count = grp[g].size();
    for(size_t a = 0; a < allele_count; ++a)
    {
      out << grp[g][a] << "/";
    }
  }
  
  out << endl;
}

void  
founders_comb_generator::write_group_seq(ostream& out, const group_seq& seq)
{
  size_t  locus_count = seq.size();

  for(size_t l = 0; l < locus_count; ++l)
  {
    out << "locus " << l << "  ";
    write_group_set(out, seq[l]);
  }
  
  out << endl;
}

ostream&
operator <<(ostream& out, founders_comb_generator::group_set grp)
{
  founders_comb_generator::write_group_set(out, grp);

  return  out;
}


//============================================================================
// IMPLEMENTATION:  founders_em_phenotype_map
//============================================================================
//
// - Last argument to the parent class (family_em_phenotype_map) constructor 
//   indicates that phenotypes should not be built.
//
founders_em_phenotype_map::founders_em_phenotype_map(output_state& ostate,
                                                 APP::Output_Streams& streams,
                                                 const FilteredMultipedigree& filtered_mped,
                                                 const pair<vector<member>, vector<member> >& members, 
                                                 const instructions& instr, 
                                                 const locus_group& loci,
                                                 const string& inner_sub_pop_name,
                                                 const string& outer_sub_pop_name)
      : family_em_phenotype_map(ostate, streams, filtered_mped, members.second, instr, loci, inner_sub_pop_name, 
                                outer_sub_pop_name, false), 
        my_alt_member_pool(members.first)
{
  build();
}

// - NOTE:  code commented 'supressed in parent class', should actually be executed in 
//   family_em_phenotype_map, but that would cause progress messages to be duplicated
//   and error messages to be out of order.
//
void
founders_em_phenotype_map::build()
{
  if(my_output_state.allow_progress_msg)
  {
    if(my_sub_pop_name.empty())
    {
      my_messages << "    Total population. " << endl;
    }
    else
    {
      my_messages << "    Subpopulation: " << my_sub_pop_name << endl;
    }

    my_messages << "      Determining possible haplotype combinations ................" << flush; 
  }
  
  build_filtered_member_pools();

  // - Suppressed in parent class.
  //
  vector<member>  family_members;
  vector<member>  singletons;
  divide_member_pool(family_members, singletons);

  // - Founder pools.
  //
  founders_generator  f_generator(my_alt_filtered_member_pool, my_loci, my_errors, my_messages,
                                  my_output_state);
  related_sub_build(f_generator.begin(), f_generator.end());
  
  // - Reps.  Suppressed in parent class.
  //
  family_generator  ps_generator(family_members, my_loci, my_errors, my_messages, my_output_state);
  related_sub_build(ps_generator.begin(), ps_generator.end());
  
  // - Singletons.  Suppressed in parent class
  //
  pheno_seq_generator  singleton_ps_generator(singletons, my_loci);
  comb_generator  singleton_hc_generator;
  singleton_sub_build(singleton_ps_generator.begin(), singleton_ps_generator.end(), 
                                        singleton_hc_generator);

  my_haplotypes.update_counts(*this);   
  
  if(my_output_state.allow_progress_msg)
  {
    if(! my_output_state.msg_interrupted)
    {
      my_messages << "... done." << endl;
    }
    else
    {
      my_messages << "                                                                  ... done." << endl;  
      my_output_state.msg_interrupted = false;
    }
  }
}

// - Members in my_member_pool and my_member_alt_pool are from a FilteredMultipedigree in which no one is actually
//   filtered.  This function produces a vector of members belonging to a FilteredMultipedigree
//   in which members who are 'uninformative leaves' have been eliminated.
//
void
founders_em_phenotype_map::build_filtered_member_pools()
{
  // - Suppressed in parent class.
  //
  build_filtered_member_pool();

  vector<member>::const_iterator  m_iter     = my_alt_member_pool.begin();
  vector<member>::const_iterator  m_end_iter = my_alt_member_pool.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    member  m = my_filtered_mped.member_find((*m_iter)->pedigree()->name(), (*m_iter)->name());
    if(m)
    {
      my_alt_filtered_member_pool.push_back(m);
    }
    else
    {
      assert(false);   // No uninformative members should be in the potential pool.
    }
  }
}

void
founders_em_phenotype_map::dump_filtered_member_pools(ostream& out) const
{
  dump_filtered_member_pool(out);

  vector<member>::const_iterator  m_iter     = my_alt_filtered_member_pool.begin();
  vector<member>::const_iterator  m_end_iter = my_alt_filtered_member_pool.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    out << (*m_iter)->pedigree()->name() << ":" << (*m_iter)->name()
        << "  subped ";
    FilteredMultipedigree::member_const_iterator  sp_iter     = (*m_iter)->subpedigree()->member_begin();
    FilteredMultipedigree::member_const_iterator  sp_end_iter = (*m_iter)->subpedigree()->member_end();
    for(; sp_iter != sp_end_iter; ++sp_iter)
    {
      out << sp_iter->name() << "  ";
    }
    
    out << endl;    
  }
}

}
}
