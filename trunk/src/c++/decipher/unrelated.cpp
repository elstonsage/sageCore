//============================================================================
// File:      unrelated.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 2/26/4                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/unrelated.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

ostream&
operator <<(ostream& out, const pheno_seq& seq)
{
  pheno_seq::const_iterator  seq_iter;
  for(seq_iter = seq.begin(); seq_iter != seq.end(); seq_iter++)
  {
    if(seq_iter->second != seq_iter->first->get_missing_phenotype_id())
    {
      out << seq_iter->first->unphased_penetrance_begin(seq_iter->second).unphased_geno().name();
    }
    else
    { 
      string  allele_name = seq_iter->first->missing_allele_name();
      out << allele_name << "/" << allele_name;
    }
    
    out << " ";
  }
  
  return  out;
}

//============================================================================
// IMPLEMENTATION:  pheno_seq_generator::iterator
//============================================================================
//
// - Build sequence of phenotypes over given genome region.
//
void
pheno_seq_generator::iterator::build_seq()
{
  if(my_member_iter == my_member_end_iter)
  {
    return;
  }
  else
  {
    const locus_group&  loci = my_generator->loci();
    size_t  locus_count = loci.size();
    for(size_t l = 0; l < locus_count; ++l)
    {
      my_pheno_seq[l].second = (*my_member_iter)->info().phenotype(loci[l].first);
    }
  }
}


//============================================================================
// IMPLEMENTATION:  comb_generator
//============================================================================
//
// - Populate a data structure, geno_seq, which consists of a group of 
//   phased genotypes corresponding to the marker phenotype at each locus.
//
// - Determine all of the haplotype combinations (pairs) consistent
//   with sequence of marker phenotypes.
//
void
comb_generator::generate(pheno_seq_generator::iterator& p_seq_iter, set<hap_seq_comb>& combs)
{
  if(x_linked(p_seq_iter))
  {
    if(male(p_seq_iter))
    {
      generate_x_linked_male(p_seq_iter, combs);
      return;
    }
    else
    {
      assert(female(p_seq_iter));        
    }
  }
  
  // - X-linked females or autosomal case.
  //
  const pheno_seq&  p_seq = *p_seq_iter;
  
  geno_seq  g_seq;
  init_geno_seq(p_seq, g_seq);
  generate(p_seq, g_seq, combs);
}

void
comb_generator::init_geno_seq(const pheno_seq& p_seq, geno_seq& g_seq)
{
  size_t  p_seq_count = p_seq.size();
  for(size_t l = 0; l < p_seq_count; ++l)
  {
    geno_set  g_set;
    const penetrance_model*  locus = p_seq[l].first;
    pen_iter iter = locus->phased_penetrance_begin(p_seq[l].second);
    pen_iter iter_end = locus->phased_penetrance_end(p_seq[l].second);
    for(; iter != iter_end; ++iter)
    {
      g_set.push_back(iter.geno_id());
    }
    
    g_seq.push_back(g_set);
  }
}

bool
comb_generator::x_linked(pheno_seq_generator::iterator& p_seq_iter)
{
  assert(! p_seq_iter->empty());

  return  (*p_seq_iter)[0].first->is_x_linked();
}

bool
comb_generator::male(pheno_seq_generator::iterator& p_seq_iter)
{
  assert(! p_seq_iter->empty());

  return  p_seq_iter.get_member()->is_male();
}

bool
comb_generator::female(pheno_seq_generator::iterator& p_seq_iter)
{
  assert(! p_seq_iter->empty());

  return  p_seq_iter.get_member()->is_female();
}

void
comb_generator::generate_x_linked_male(pheno_seq_generator::iterator& p_seq_iter, set<hap_seq_comb>& combs)
{
  const pheno_seq&  p_seq = *p_seq_iter;
  x_allele_seq  x_seq;

  size_t  locus_count = p_seq.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    x_allele_set  x_set;
    
    const penetrance_model*  locus = p_seq[l].first;
    size_t                   phenotype_index = p_seq[l].second;
    if(phenotype_index == locus->get_missing_phenotype_id())
    {
      const string&  missing_allele_name = locus->missing_allele_name();
    
      allele_iterator  a_iter     = locus->allele_begin();
      allele_iterator  a_end_iter = locus->allele_end();
      for(; a_iter != a_end_iter; ++a_iter)
      {
        const string&  allele_name = a_iter->name();
        if(allele_name != "~Y" && allele_name != missing_allele_name)
        {
          x_set.push_back(a_iter->id());
        }
      }
    }
    else
    {
      assert(locus->unphased_penetrance_count(phenotype_index) == 1);
      
      unphased_genotype  geno = locus->unphased_penetrance_begin(phenotype_index).unphased_geno();

      allele  a1 = geno.allele1();
      allele  a2 = geno.allele2();
      const string  a1_name = a1.name();
      const string  a2_name = a2.name();
      size_t  a1_id = a1.id();
      size_t  a2_id = a2.id();
        
      if(a1_name != "~Y")
      {
        assert(a2_name == "~Y" || a2_name == a1_name);
        x_set.push_back(a1_id);
      }
      else
      {
        assert(a2_name != "~Y");
        x_set.push_back(a2_id);
      }
    }
    
    /*
    for(size_t i = 0; i < x_set.size(); ++i)
    {
      cout << "allele set " << (p_seq[l].first)->get_allele(x_set[i]).name() << " ";
    }
    cout << endl;
    */
    
    x_seq.push_back(x_set);
  }
  
  generate_x_linked_male(x_seq, combs);
}

// - Generate haplotype combinations, "monotypes", recursively by breaking 
//   first set of alleles encountered into a set with one allele and a set 
//   with the remaining alleles.  Repeat the process until there is one allele
//   at every locus, ie, the haplotype is unambiguous.
//
void  
comb_generator::generate_x_linked_male(x_allele_seq& a_seq, set<hap_seq_comb>& combs)
{
  size_t  locus_count = a_seq.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    while(a_seq[l].size() > 1)
    {
      x_allele_seq  new_a_seq(a_seq);
      x_allele_set  new_a_set;
      
      new_a_set.push_back(a_seq[l][a_seq[l].size() - 1]);
      new_a_seq[l] = new_a_set;
      
      a_seq[l].pop_back();
      generate_x_linked_male(new_a_seq, combs);
    }
  }
  
  hap_seq  h_seq;
  for(size_t l = 0; l < locus_count; ++l)
  {
    assert(a_seq[l].size() == 1);
    h_seq.push_back(a_seq[l][0]);
  }
  
  hap_seq_comb  comb;
  comb.insert(h_seq);
  combs.insert(comb);
}

// - Generate haplotype combinations recursively by breaking first set of genotypes
//   encountered into a set with one genotype and a set with the remaining
//   genotypes.  Repeat the process until there is one genotype at
//   every locus, ie, the haplotypes are unambiguous.
//
void
comb_generator::generate(const pheno_seq& p_seq,
                         geno_seq& g_seq, set<hap_seq_comb>& combs)
{
  size_t  g_seq_count = g_seq.size();
  for(size_t l = 0; l < g_seq_count; ++l)
  {
    while(g_seq[l].size() > 1)   
    {
      geno_seq  new_g_seq(g_seq);
      geno_set  new_g_set;
      
      new_g_set.push_back(g_seq[l][g_seq[l].size() - 1]);
      new_g_seq[l] = new_g_set;
      
      g_seq[l].pop_back();
      generate(p_seq, new_g_seq, combs);
    }
  }
  
  // - Note combinations are multsets, hence elements are sorted.
  //   combs is a set and therefore will not allow a duplicate combinations
  //   to be inserted regardless of the haplotype order used to construct them.
  //
  combs.insert(comb(p_seq, g_seq));
}

// - Create a combination of haplotypes from a sequence of genotypes sets
//   containing exactly one phased genotype at every locus.
//
hap_seq_comb
comb_generator::comb(const pheno_seq& p_seq, const geno_seq& g_seq)
{
  hap_seq  hap1;
  hap_seq  hap2;
  
  size_t  p_seq_count = p_seq.size();
  for(size_t l = 0; l < p_seq_count; ++l)
  {
    phased_genotype  g = p_seq[l].first->get_phased_genotype(g_seq[l][0]);
    hap1.push_back(g.allele1().id());
    hap2.push_back(g.allele2().id());    
  }
  
  hap_seq_comb  c;
  
  c.insert(hap1);
  c.insert(hap2);
  
  return  c;
}


//============================================================================
// IMPLEMENTATION:  unrelated_em_phenotype_map
//============================================================================
//
void
unrelated_em_phenotype_map::build()
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
  
  pheno_seq_generator  ps_generator(my_member_pool, my_loci);
  comb_generator  hc_generator;
  
  pheno_seq_generator::iterator  iter = ps_generator.begin();
  pheno_seq_generator::iterator  iter_end = ps_generator.end();
  for(; iter != iter_end; ++iter)
  {
    if(my_phenotype_directory.find(*iter) == my_phenotype_directory.end())
    {
      set<hap_seq_comb>  combs;
      
      try
      {
        hc_generator.generate(iter, combs);
      }
      catch(const bad_alloc&)
      {
        if(! my_output_state.msg_interrupted)
        {
          my_messages << endl;
          my_output_state.msg_interrupted = true;
        }
        
        my_errors << priority(error) << "Not enough memory available to process "
                  << "member " << iter.get_member()->name() << " in pedigree "
                  << iter.get_member()->pedigree()->name() << ".  Skipping this individual ..."
                  << endl;
        continue;
      }
      
      try
      {
        my_member_directory[iter.get_member()] = *iter; 
                 
        assert(! combs.empty());
        my_total_hap_count += combs.begin()->size();
        
        my_phenotypes.push_back(em_phenotype(&my_haplotypes, combs));
        my_phenotype_directory[*iter] = my_phenotypes.size() - 1;
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
    else
    {
      my_member_directory[iter.get_member()] = *iter; 
             
      my_phenotypes[my_phenotype_directory[*iter]].incr();
      my_total_hap_count += my_phenotypes[my_phenotype_directory[*iter]].comb_size();
    }
  }
   
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

set<member, member_order<member> >
unrelated_em_phenotype_map::members() const
{
  set<member, member_order<member> >  ordered_members;
  map<member, pheno_seq>::const_iterator  m_iter     = my_member_directory.begin();
  map<member, pheno_seq>::const_iterator  m_end_iter = my_member_directory.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    ordered_members.insert(m_iter->first);
  }
  
  return  ordered_members;
}

void
unrelated_em_phenotype_map::dump(ostream& out) const
{
  out << "Haplotypes:" << endl;
  my_haplotypes.dump(out, my_total_hap_count, my_loci);

  map<pheno_seq, size_t>::const_iterator  p_iter     = my_phenotype_directory.begin();
  map<pheno_seq, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    out << "phenotype  " << p_iter->first << endl;
    
    size_t  phenotype_index = p_iter->second;
    assert(phenotype_index < my_phenotypes.size());
    my_phenotypes[phenotype_index].dump(out, my_loci, my_haplotypes);
  }
}


}
}
