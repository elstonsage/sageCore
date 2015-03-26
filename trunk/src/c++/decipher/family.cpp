//============================================================================
// File:      family.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 8/11/4                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/family.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

ostream&
operator <<(ostream& out, const aps_seq& seq)
{
  for(size_t l = 0; l < seq.size(); ++l)
  {
    out << seq[l].first.name() << "  ";
    
    set<allele_pair>::const_iterator  iter = seq[l].second.begin();
    set<allele_pair>::const_iterator  end_iter = seq[l].second.end();
    for(; iter != end_iter; ++iter)
    {
      out << "(" << (iter->first == MLOCUS::NPOS ? "no info" : seq[l].first.get_allele(iter->first).name());
      out << ", ";
      out << (iter->second == MLOCUS::NPOS ? "no info" : seq[l].first.get_allele(iter->second).name()) << ")  " ;
    }
    
    out << endl;
  }
  
  return  out;
}

void
write_allele_pair(ostream& out, MLOCUS::inheritance_model* locus, allele_pair ap)
{
  string  first_allele  = ap.first  == MLOCUS::NPOS ? "no info" : locus->get_allele(ap.first).name();
  string  second_allele = ap.second == MLOCUS::NPOS ? "no info" : locus->get_allele(ap.second).name();  

  out << first_allele << ", " << second_allele;
}


//============================================================================
// IMPLEMENTATION:  family_generator::iterator
//============================================================================
//
void
family_generator::iterator::dump_phenotypes(ostream& out) const
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

// - Build sequence of allele pair sets.
//
bool
family_generator::iterator::build_combs()
{
  if(my_member_iter == my_member_end_iter)
  {
    return  false;
  }
  else
  {
    bool& interrupted = my_generator->messages_interrupted();
    my_hap_seq_combs.clear();  
    init_aps_seq();
    size_t  marker_count = my_aps_seq.size();

    #if 0
      cout << "member " << (*my_member_iter)->name() << endl;
      cout << "pedigree " << (*my_member_iter)->pedigree()->name() << endl;    
      cout << "allele pair sequence after initialization\n" << my_aps_seq << endl;    
      dump_phenotypes(cout);
    #endif
    
    // - Create founder allele graphs.
    //
    vector<boost::shared_ptr<FounderAlleleGraph> >  graphs;
    McmcMeiosisMap  mm(*my_subped_ptr);
        
    size_t  bit_count = mm.get_meiosis_count();
    
    graphs.reserve(marker_count);
    for(size_t l = 0; l < marker_count; ++l)
    {
      // - inheritance model
      //
      #if 0
        my_aps_seq[l].first.print_info_sparse_matrix();
      #endif
    
      if(my_aps_seq[l].second.size())   // Mendelian inconsistency was found during genotype elimination.
      {
        const string&  locus_name = my_aps_seq[l].first.name();
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
        graphs.push_back(boost::shared_ptr<FounderAlleleGraph>(new FounderAlleleGraph(mm, my_aps_seq[l].first)));
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
        reset_aps_seq();
        
        for(size_t l = 0; l < graphs.size(); ++l)
        {
          //cout << "locus " << l << endl; 
          if(graphs[l] == SKIP_MARKER)
          {
            continue;
          }
        
          for(size_t vas = 0; vas < graphs[l]->get_valid_allele_set_count(); ++vas)
          {
            std::pair<MLOCUS::allele, MLOCUS::allele> ap = 
                    graphs[l]->get_individual_allele_pattern(vas, **my_member_iter);
            
            if(ap.first.is_valid())
            {
              assert(ap.first.name() != "~remap");
            }
            
            if(ap.second.is_valid())
            {
              assert(ap.second.name() != "~remap");
            }            
            
            //cout << " allele pair ";
            //write_allele_pair(cout, &(my_aps_seq[l].first), ap);
            allele_pair ap1 = std::make_pair(ap.first.id(), ap.second.id());
            my_aps_seq[l].second.insert(ap1);
          }
          
          //cout << endl;
        }
      
        // - Consider haplotype combinations pertaining to one inheritance vector (pattern)
        //   at a time.
        //
        build_combinations(interrupted);
      }
    }
    
    return  valid_pattern_found;
  }
}

// - Generate haplotype combinations for corresponding to current allele pair
//   set sequence.
//
void
family_generator::iterator::build_combinations(bool& messages_interrupted) throw(bad_alloc)
{
  family_comb_generator  comb_gen;

  try
  {
    comb_gen.generate(my_aps_seq, my_hap_seq_combs);
  }
  catch(const bad_alloc&)
  {
    if(! messages_interrupted)
    {
      my_generator->messages() << endl;
      messages_interrupted = true;
    }
    
    my_generator->errors() << priority(error) << "Not enough memory available to process member "
              << (*my_member_iter)->name() << " in pedigree "
              << (*my_member_iter)->pedigree()->name() << ".  Skipping member "
              << "and its constituent pedigree ... " << endl;
    
    throw;
  }
}

// - Initialize allele pair sequence.
//
void
family_generator::iterator::init_aps_seq()
{
  my_inconsistent_loci.clear();
  my_inconsistent_group_loci.clear();

  const locus_group&  loci = my_generator->loci();
  size_t  locus_count = loci.size();
  my_aps_seq.resize(locus_count);
  
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
    
    my_aps_seq[l].first = ig(*my_subped_ptr, loci[l].first);
    my_aps_seq[l].second.clear();    
    
    if(ig.inconsistent())
    {
      // - Use original global inheritance model.  Will be skipped by build_combs().
      //
      my_aps_seq[l].first = *(loci[l].second);
      my_aps_seq[l].second.insert(make_pair(MLOCUS::NPOS, MLOCUS::NPOS));
      my_inconsistent_loci.push_back(loci[l].first);
      my_inconsistent_group_loci.push_back(l);      
    }
  }
}

void
family_generator::iterator::reset_aps_seq()
{
  aps_seq::iterator  iter     = my_aps_seq.begin();
  aps_seq::iterator  end_iter = my_aps_seq.end();
  for(; iter != end_iter; ++iter)
  {
    iter->second.clear();
  }
  
  // - Restore 'inconsistency tag'.
  //
  size_t inconsistent_count = my_inconsistent_group_loci.size();
  for(size_t l = 0; l < inconsistent_count; ++l)
  {
    my_aps_seq[my_inconsistent_group_loci[l]].second.insert(make_pair(MLOCUS::NPOS, MLOCUS::NPOS));
  }
}

//============================================================================
// IMPLEMENTATION:  family_comb_generator
//============================================================================
//
void
family_comb_generator::init_pair_seq(const aps_seq& a_seq, pair_seq& pr_seq)
{
  size_t  locus_count = a_seq.size();
  pr_seq.assign(locus_count, pair_set());
  
  for(size_t l = 0; l < locus_count; ++l)
  {
    pr_seq[l].assign(a_seq[l].second.begin(), a_seq[l].second.end());
    expand_unknown_alleles(pr_seq[l], a_seq[l].first);
  }
}

// - Replace pairs containing MLOCUS::NPOS (no information) as one or both members with
//   with all possible pairs and eliminate duplicate pairs.
// 
void
family_comb_generator::expand_unknown_alleles(pair_set& pr_set, const penetrance_model& locus)
{
  pair_set  new_pr_set;

  for(size_t p = 0; p < pr_set.size(); ++p)
  {
    size_t  pair_unknowns = unknowns(pr_set[p]);
    switch(pair_unknowns)
    {
      case none:
      if(is_new_pair(new_pr_set, pr_set[p]))
      {
        new_pr_set.push_back(pr_set[p]);
      }
        
      break;
        
      case first_allele:
      case second_allele:
      {
        allele_iterator  allele_iter = locus.allele_begin();
        allele_iterator  allele_end_iter = locus.allele_end();
        for(; allele_iter != allele_end_iter; ++allele_iter)
        {
          if(pair_unknowns == first_allele)
          {
            allele_pair  pr = make_pair(allele_iter->id(), pr_set[p].second);
            if(is_new_pair(new_pr_set, pr))
            {
              new_pr_set.push_back(pr);
            }
          }
          else
          {
            allele_pair  pr = make_pair(pr_set[p].first, allele_iter->id());
            if(is_new_pair(new_pr_set, pr))
            {
              new_pr_set.push_back(pr);
            }          
          }          
        }        
        
        break;
      }
        
      case both:
      {
        allele_iterator  first_allele_iter = locus.allele_begin();
        allele_iterator  first_allele_end_iter = locus.allele_end();
        for(; first_allele_iter != first_allele_end_iter; ++first_allele_iter)
        {
          allele_iterator  second_allele_iter = locus.allele_begin();
          allele_iterator  second_allele_end_iter = locus.allele_end();
          for(; second_allele_iter != second_allele_end_iter; ++second_allele_iter)
          {
            allele_pair  pr = make_pair(first_allele_iter->id(), second_allele_iter->id());
            if(is_new_pair(new_pr_set, pr))
            {
              new_pr_set.push_back(pr);
            }
          }        
        }
      
        break;
      }
        
      default:
        assert(false);
    }
  }
  
  pr_set = new_pr_set;
}

family_comb_generator::UNKNOWNS  
family_comb_generator::unknowns(const allele_pair& pair)
{
  size_t  unks = none;

  if(pair.first == MLOCUS::NPOS)
  {
    unks += first_allele;
  }
  
  if(pair.second == MLOCUS::NPOS)
  {
    unks += second_allele;
  }
  
  return  static_cast<UNKNOWNS>(unks);
}

bool  
family_comb_generator::is_new_pair(const pair_set& pr_set, const allele_pair& pr)
{
  pair_set::const_iterator  begin_iter = pr_set.begin();
  pair_set::const_iterator  end_iter   = pr_set.end();
  
  return  find(begin_iter, end_iter, pr) == end_iter;
}

void
family_comb_generator::generate(const aps_seq& a_seq, set<hap_seq_comb>& combs)
{
  pair_seq  pr_seq;
  init_pair_seq(a_seq, pr_seq);
  generate(a_seq, pr_seq, combs);
}

// - Generate haplotypes recursively by breaking first set of genotypes
//   encountered into a set with one genotype and a set with the remaining
//   genotypes.  The process is repeated until there is one genotype at
//   every locus, ie, the haplotypes are unambiguous.
//
void
family_comb_generator::generate(const aps_seq& a_seq,
                         pair_seq& pr_seq, set<hap_seq_comb>& combs)
{
  for(size_t l = 0; l < pr_seq.size(); ++l)
  {
    while(pr_seq[l].size() > 1)   
    {
      pair_seq  new_pr_seq(pr_seq);
      pair_set  new_pr_set;
      
      new_pr_set.push_back(pr_seq[l][pr_seq[l].size() - 1]);
      new_pr_seq[l] = new_pr_set;
      
      pr_seq[l].pop_back();
      generate(a_seq, new_pr_seq, combs);
    }
  }
  
  // - Note combinations are multsets, hence elements are sorted.
  //   combs is a set and therefore will not allow a duplicate combinations
  //  to be inserted regardless of the haplotype order used to construct them.
  //  
  combs.insert(comb(a_seq, pr_seq));
}

// - Create a combination of haplotypes from a sequence of allele pair sets
//   containing exactly one allele pair at every locus.
//
hap_seq_comb
family_comb_generator::comb(const aps_seq& a_seq, const pair_seq& pr_seq)
{
  size_t  locus_count = a_seq.size();
  assert(pr_seq.size() == locus_count);

  hap_seq  hap1;
  hap_seq  hap2;
  
  for(size_t l = 0; l < locus_count; ++l)
  {
    assert(! pr_seq[l].empty());
    hap1.push_back(pr_seq[l][0].first);
    hap2.push_back(pr_seq[l][0].second);    
  }
  
  hap_seq_comb  c;
  
  c.insert(hap1);
  c.insert(hap2);
  
  return  c;
}

void  
family_comb_generator::write_pair_set(ostream& out, const pair_set& ps)
{
  for(size_t i = 0; i < ps.size(); ++i)
  {
    out << ps[i].first << "/" << ps[i].second << " ";
  }
  
  out << endl;
}

void  
family_comb_generator::write_pair_seq(ostream& out, const pair_seq& ps)
{
  for(size_t l = 0; l < ps.size(); ++l)
  {
    out << "locus " << l << "  ";
    write_pair_set(out, ps[l]);
  }
  
  out << endl;
}

ostream&
operator <<(ostream& out, family_comb_generator::pair_set ps)
{
  size_t  pair_count = ps.size();
  for(size_t p = 0; p < pair_count; ++p)
  {
    out << ps[p].first << ",  " << ps[p].second << endl;
  }
  
  out << endl;
  
  return  out;
}

//============================================================================
// IMPLEMENTATION:  family_em_phenotype_map
//============================================================================
//
family_em_phenotype_map::family_em_phenotype_map(output_state& ostate,
                                                 APP::Output_Streams& streams,
                                                 const FilteredMultipedigree& filtered_mped,
                                                 const vector<member>& members, 
                                                 const instructions& instr, 
                                                 const locus_group& loci,
                                                 const string& inner_sub_pop_name,
                                                 const string& outer_sub_pop_name,
                                                 bool build_phenotypes)
      : member_em_phenotype_map(ostate, streams, members, loci, inner_sub_pop_name, outer_sub_pop_name), 
        my_filtered_mped(filtered_mped), my_instructions(instr)
{
  // - Remove uninformative "leaves" from data (unless they are only children!).  This could more
  //   efficiently be done at a higher level, but there is not a clean way to do it.  Don't want to
  //   filter before partitioning as user may designate anyone as a family rep and everyone can 
  //   potentially have partitioning information.  Could do it after parttitioning in analysis::do_tasks(),
  //   but this is templatized on phenotype map type and filtration is only needed for family and founder
  //   phenotype maps.
  //
  FPED::has_informative_loci<FilteredMultipedigree::member_type>  base_filter(filtered_mped, false);
  init_filter(base_filter);
  not_uninformative_leaf  filter(base_filter);

  FPED::FilterResults  results = FPED::MPFilterer::add_multipedigree_filtered_by_members(my_filtered_mped, 
                                                                                          filtered_mped, filter);
  
  my_filtered_mped.construct();
  //dump_filtered_mped(cout);

  if(build_phenotypes)
  {
    build();
  }
  
  //dump(cout);
}

// - Add loci to has_informative_loci filter as needed.
//
void
family_em_phenotype_map::init_filter(FPED::has_informative_loci<FPED::FilteredMultipedigree::member_type>& base_filter) const
{
  size_t  locus_count = my_loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    base_filter.set_check_status_for_locus(my_loci[l].first, true);
  }            
}

void
family_em_phenotype_map::build()
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
 
  build_filtered_member_pool();
  
  vector<member>  family_members;
  vector<member>  singletons;
  divide_member_pool(family_members, singletons);

  family_generator  ps_generator(family_members, my_loci, my_errors, my_messages, my_output_state);
  related_sub_build(ps_generator.begin(), ps_generator.end());
  
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

// - Members in my_member_pool are from a FilteredMultipedigree in which no one is actually
//   filtered.  This function produces a vector of members belonging to a FilteredMultipedigree
//   in which members who are 'uninformative leaves' have been eliminated.
//
void
family_em_phenotype_map::build_filtered_member_pool()
{
  vector<member>::const_iterator  m_iter     = my_member_pool.begin();
  vector<member>::const_iterator  m_end_iter = my_member_pool.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    member  m = my_filtered_mped.member_find((*m_iter)->pedigree()->name(), (*m_iter)->name());
    if(m)
    {
      my_filtered_member_pool.push_back(m);
    }
    else
    {
      assert(false);   // No uninformative members should be in the potential pool.
    }
  }
}

void
family_em_phenotype_map::divide_member_pool(vector<member>& family_members, 
                                            vector<member>& singletons     ) const
{
  size_t  member_count = my_filtered_member_pool.size();
  for(size_t m = 0; m < member_count; ++m)
  {
    if(MPED::mp_utilities::nuclear_family_count(my_filtered_member_pool[m]) == 0)
    {
      singletons.push_back(my_filtered_member_pool[m]);
    }
    else
    {
      family_members.push_back(my_filtered_member_pool[m]);
    }
  }
  
  /*
  for(size_t m = 0; m < family_members.size(); ++m)
  {
    cout << family_members[m]->pedigree()->name() << ":" << family_members[m]->name() << endl;
  }
  
  for(size_t m = 0; m < singletons.size(); ++m)
  {
    cout << singletons[m]->pedigree()->name() << ":" << singletons[m]->name() << endl;
  } 
  */ 
}

void  
family_em_phenotype_map::singleton_sub_build(pheno_seq_generator::iterator iter, 
                                             pheno_seq_generator::iterator end_iter, 
                                             comb_generator& comb_gen)
{
  for(; iter != end_iter; ++iter)
  {
    set<hap_seq_comb>  combs;
    
    try
    {
      comb_gen.generate(iter, combs);
    }
    catch(const bad_alloc&)
    {
      if(! my_output_state.msg_interrupted)
      {
        my_messages << endl;
        my_output_state.msg_interrupted = true;
      }
      
      my_errors << priority(error) << "Not enough memory available to process member "
                << iter.get_member()->name() << " in pedigree "
                << iter.get_member()->pedigree()->name() << ".  Skipping member ..." << endl;
      continue;    
    }

    try
    {
      assert(! combs.empty());
      my_total_hap_count += combs.begin()->size();
      my_phenotypes.push_back(em_phenotype(&my_haplotypes, combs));
      my_phenotype_directory.insert(make_pair(iter.get_member(), my_phenotypes.size() - 1));
    }
    catch(const bad_alloc&)
    {
      if(! my_output_state.msg_interrupted)
      {
        my_messages << endl;
        my_output_state.msg_interrupted = true;
      }
      
      my_errors << priority(fatal) << "Not enough memory available to process "
                << "the data.  Currently processing member " << iter.get_member()->name() 
                << " in pedigree " << iter.get_member()->pedigree()->name() 
                << ".  Exiting program ..." <<  endl;
      exit(1);
    }            
  }
}

set<member, member_order<member> >
family_em_phenotype_map::members() const
{
  set<member, member_order<member> >  ordered_members;
  map<member, size_t>::const_iterator  m_iter     = my_phenotype_directory.begin();
  map<member, size_t>::const_iterator  m_end_iter = my_phenotype_directory.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    ordered_members.insert(m_iter->first);
  }
  
  return  ordered_members;  
}

void
family_em_phenotype_map::dump_filtered_member_pool(ostream& out) const
{
  vector<member>::const_iterator  m_iter     = my_filtered_member_pool.begin();
  vector<member>::const_iterator  m_end_iter = my_filtered_member_pool.end();
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

void
family_em_phenotype_map::dump(ostream& out) const
{
  out << "Haplotypes:" << endl;
  my_haplotypes.dump(out, my_total_hap_count, my_loci);

  map<member, size_t>::const_iterator  p_iter     = my_phenotype_directory.begin();
  map<member, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    out << "pedigree   " << p_iter->first->pedigree()->name()
        << "  member  " << p_iter->first->name() << endl;
        
    size_t  phenotype_index = p_iter->second;
    assert(phenotype_index < my_phenotypes.size());
        
    my_phenotypes[phenotype_index].dump(out, my_loci, my_haplotypes);
  }
}

void
family_em_phenotype_map::dump_filtered_mped(ostream& out) const
{
  out << "\nPed\tName\tP1\tP2\tSex" << endl;
  
  FilteredMultipedigree::pedigree_const_iterator  p_iter     = my_filtered_mped.pedigree_begin();
  FilteredMultipedigree::pedigree_const_iterator  p_end_iter = my_filtered_mped.pedigree_end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    string  pedigree_name = p_iter->name();
  
    // - Families
    //
    FilteredMultipedigree::subpedigree_const_iterator  subp_iter     = p_iter->subpedigree_begin();
    FilteredMultipedigree::subpedigree_const_iterator  subp_end_iter = p_iter->subpedigree_end();
    for(; subp_iter != subp_end_iter; ++subp_iter)
    {
      out << endl;
      FilteredMultipedigree::member_const_iterator  memb_iter     = subp_iter->member_begin();
      FilteredMultipedigree::member_const_iterator  memb_end_iter = subp_iter->member_end();
      for(; memb_iter != memb_end_iter; ++memb_iter)
      {
        string  member_name = memb_iter->name();
        FilteredMultipedigree::member_const_pointer  p1 = memb_iter->parent1();
        FilteredMultipedigree::member_const_pointer  p2 = memb_iter->parent2();        
        string  p1_name = p1 ? p1->name() : "0";
        string  p2_name = p2 ? p2->name() : "0";
        size_t  sex = memb_iter->get_effective_sex();
        out << pedigree_name << "\t" << member_name << "\t" << p1_name << "\t"
            << p2_name << "\t" << sex << endl;
      }      
    }

    // - Singletons.
    //
    FilteredMultipedigree::member_const_iterator  memb_iter     = p_iter->member_begin();
    FilteredMultipedigree::member_const_iterator  memb_end_iter = p_iter->member_end();
    for(; memb_iter != memb_end_iter; ++memb_iter)
    {
      if(! (memb_iter->parent1() || memb_iter->offspring_count()))
      {
        out << endl;
        
        string  member_name = memb_iter->name();
        string  p1_name = "0";
        string  p2_name = "0";
        size_t  sex = memb_iter->get_effective_sex();
        out << pedigree_name << "\t" << member_name << "\t" << p1_name << "\t"
            << p2_name << "\t" << sex << endl;        
      }
    }        
  }  
}

}
}
