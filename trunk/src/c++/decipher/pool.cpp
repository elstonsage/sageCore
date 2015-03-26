//============================================================================
// File:      pool.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 1/9/6                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2006 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/pool.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

pool_pheno_seq_generator::iterator::reconcile  RECONCILIATION_METHOD = 
                                      pool_pheno_seq_generator::iterator::OVER_UNDER;    

ostream&  
operator <<(ostream& out, const pool_pheno pheno)
{
  out << "{";
  
  string  allele_name = "";
  size_t  allele_count = (size_t)(-1);
  
  allele_counts  counts = pheno.second;
  const MLOCUS::inheritance_model*  locus = pheno.first;
  
  allele_counts::const_iterator  c_iter     = counts.begin();
  allele_counts::const_iterator  c_end_iter = counts.end();
  for(; c_iter != c_end_iter; ++c_iter)
  {
    if(! allele_name.empty())
    {
      out << allele_name << ":" << allele_count << ", ";
    }
  
    allele_name = locus->get_allele(c_iter->first).name();
    allele_count = c_iter->second;
  }
 
  if(! allele_name.empty())
  {
    out << allele_name << ":" << allele_count;
  }
  
  out << "}";
  
  return  out;
}

ostream&
operator <<(ostream& out, const pool_pheno_seq& seq)
{
  pool_pheno_seq::const_iterator  seq_iter     = seq.begin();
  pool_pheno_seq::const_iterator  seq_end_iter = seq.end();  
  for(; seq_iter != seq_end_iter; ++seq_iter)
  {
    out << *seq_iter << " ";
  }
  
  return  out;
}

void
write_allele_weights(ostream& out, const allele_weights& weights, const MLOCUS::inheritance_model* locus)
{
  allele_weights::const_iterator  w_iter     = weights.begin();
  allele_weights::const_iterator  w_end_iter = weights.end();
  for(; w_iter != w_end_iter; ++w_iter)
  {
    out << "allele " << locus->get_allele(w_iter->first).name() << ", weight "
        << w_iter->second << endl;
  }
}

void
write_allele_counts(ostream& out, const allele_counts& counts, const MLOCUS::inheritance_model* locus)
{ 
  out << make_pair(locus, counts);
}

//============================================================================
// IMPLEMENTATION:  pool_pheno_seq_generator::iterator
//============================================================================
//
// - Build sequence of phenotypes over the haplotyping region.
//
bool
pool_pheno_seq_generator::iterator::build_seq()
{
  assert(my_member_iter != my_member_end_iter);
  
  my_expanded_pheno_seq_count = 1;
  bool  any_info = false;

  const locus_group&  loci = my_generator->loci();
  size_t  pool_size = get_pool_size();
  
  size_t  locus_count = loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    allele_counts&  counts = my_pool_pheno_seq[l].second;
    
    counts.clear();
    build_pool_phenotype(l, counts, pool_size);
    
    bool  informative_locus = ! counts.empty();
    
    any_info |= informative_locus;
    
    if(! informative_locus)
    {
      my_expanded_pheno_seq_count *= genotype_count(loci[l].second->allele_count(), pool_size); 
    }
  }
    
  // cout << "\nPool " << pool_name() << " has " << my_expanded_pheno_seq_count 
  //      << " possible phenotypes." << endl; 
  
  if(any_info  /* && my_expanded_pheno_seq_count <= MAX_EXPANDED_PHENO_SEQ_COUNT */)
  {
    return  true;
  }
  else
  {
    cerrorstream&  errors = my_generator->errors();

    bool&  interrupted = my_generator->messages_interrupted();
    if(! interrupted)
    {
      my_generator->messages() << endl;
      interrupted = true;
    }
    
    string  error_message;
    
    errors << priority(information) << "Pool, " << pool_name()
           << ", contains " << /* (any_info ? "insufficient" : */ "no" /*)*/ 
           << " information.  Skipping this pool ..." << endl; 
    
    return  false;
  }
}

// - Convert trait values to allele counts for the current record at the given
//   locus.  If the locus contains no usable information, leave it empty.
//
void  
pool_pheno_seq_generator::iterator::build_pool_phenotype(size_t l, allele_counts& counts, size_t pool_size)
{
  allele_weights  weights;
  
  get_allele_weights(l, weights);
  
  //const MLOCUS::inheritance_model*  locus = my_generator->instr().loci[l].second;
  //cout << "locus " << locus->name() << endl; 
  //write_allele_weights(cout, weights, locus);
  
  get_allele_counts(l, weights, counts, pool_size);
  //write_allele_counts(cout, counts, my_generator->instr().loci[l].second);  
}

// - 'Allele weights' are the values of traits that the user designated in the analysis
//   block as containing allele information for the given locus.  If information is 
//   missing for any allele or sum of weights is greater than one, leave weights
//   vector empty.  Weight of last allele is 1 - sum of weights of all other alleles.
//   
// - Issue warning if sum of weights is greater than one.
//
void  
pool_pheno_seq_generator::iterator::get_allele_weights(size_t l, allele_weights& weights)
{
  cerrorstream&  errors = my_generator->errors();

  const FPED::FilteredMemberInfo&  member_info = (*my_member_iter)->info();
  const MLOCUS::inheritance_model*  locus = my_generator->instr().loci[l].second;
  const set<instructions::pool_locus::allele>&  alleles = my_generator->instr().pool_loci[l].alleles;
  
  double  allele_weight_total = 0.0;
  string  last_allele_name = "";
  
  set<instructions::pool_locus::allele>::const_iterator  a_iter     = alleles.begin();
  set<instructions::pool_locus::allele>::const_iterator  a_end_iter = alleles.end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    size_t  trait_index = a_iter->index;
    if(trait_index == (size_t)(-1))        // Last allele has no associated trait;
    {
      assert(last_allele_name.empty());     // Only one last allele per locus!
      last_allele_name = a_iter->name;
      continue;
    }
    
    double  trait_value = member_info.trait(trait_index);
    if(! SAGE::isnan(trait_value))
    {
      if(! (0.0 <= trait_value && trait_value <= 1.0))
      {
        bool&  interrupted = my_generator->messages_interrupted();
        if(! interrupted)
        {
          my_generator->messages() << endl;
          interrupted = true;
        }  
      
        errors << priority(warning) << "An allele weight at locus, " << locus->name() 
               << ", for pool, " << pool_name() << ", is greater than 1 or less than 0. "
               << "Setting values to missing at this locus for this pool ..." << endl;
        weights.clear();
        
        return;
      } 
    
      allele_weight_total += trait_value;
      
      size_t  allele_id = locus->get_allele(a_iter->name).id();
      assert(allele_id != MLOCUS::NPOS);
      weights[allele_id] = trait_value; 
    }
    else       // All traits at a locus must have values for information to be meaningful.
    {
      weights.clear();
      
      return;
    } 
  }
  
  if(allele_weight_total <= 1.0)
  {
    assert(! last_allele_name.empty());
    size_t  allele_id = locus->get_allele(last_allele_name).id();
    assert(allele_id != MLOCUS::NPOS);
    weights[allele_id] = 1.0 - allele_weight_total; 
  }
  else
  {
    bool&  interrupted = my_generator->messages_interrupted();
    if(! interrupted)
    {
      my_generator->messages() << endl;
      interrupted = true;
    }  
  
    errors << priority(warning) << "Sum of allele weights at locus, " << locus->name() 
           << ", for pool, " << pool_name() << ", greater than 1. "
           << "Setting values to missing at this locus for this pool ..." << endl;
    weights.clear();
  }
}

string
pool_pheno_seq_generator::iterator::pool_name() const
{
  return  (*my_member_iter)->pedigree()->name() + ":" + (*my_member_iter)->name();
}
 
// - Knowing pool size, convert allele weights to allele counts for the given locus.
//   
void  
pool_pheno_seq_generator::iterator::get_allele_counts(size_t l, const allele_weights& weights, 
                                                         allele_counts& counts, size_t pool_size)
{
  assert(counts.empty());
  
  if(weights.empty())
  {
    return;                // Nothing to do.
  }
  
  size_t  total_allele_count = 0;

  allele_weights::const_iterator  w_iter     = weights.begin();
  allele_weights::const_iterator  w_end_iter = weights.end();
  for(; w_iter != w_end_iter; ++w_iter)
  {
    size_t  allele_count = round_2_integer(w_iter->second * pool_size);
    counts[w_iter->first] = allele_count;
    total_allele_count += allele_count;
  }
  
  if(total_allele_count != pool_size)
  {
    reconcile_counts(l, counts, pool_size, RECONCILIATION_METHOD);

    /*
    cerrorstream&  errors = my_generator->errors();
    errors << priority(information) << "Unable to determine allele counts for pool, "
           << pool_name() << ", at locus, " << my_generator->instr().loci[l].second->name()
           << ".  Setting allele counts to missing for this pool at this locus ..." << endl;
    counts.clear();
    */  
  }
}

void
pool_pheno_seq_generator::iterator::reconcile_counts(size_t l, allele_counts& counts, 
                                                     size_t pool_size, reconcile method)
{
  size_t  allele_count = count_alleles(counts);

  while(allele_count != pool_size)
  {
    //cout << "\npre-reconciliation   ";
    //write_allele_counts(cout, counts, my_generator->instr().loci[l].second);
    //cout << endl;
  
    allele_counts::iterator  adjustee;
    
    switch(method)
    {
      case GREATEST:
        adjustee = greatest(counts);
        break;
        
      case LEAST:
        adjustee = least(counts);
        break;
        
      case RANDOM:
        adjustee = random(counts);
        break;
        
      case OVER_UNDER:
        adjustee = over_under(counts, allele_count > pool_size);
          break;
          
      default:
        assert(false);
    }    
    
    if(allele_count > pool_size)
    {
      --(adjustee->second);
      --allele_count;
    }
    else
    {
      ++(adjustee->second);
      ++allele_count;
    }
    
    //cout << "post-reconciliation  ";
    //write_allele_counts(cout, counts, my_generator->instr().loci[l].second);
    //cout << endl;    
    
  }
}

size_t  
pool_pheno_seq_generator::iterator::count_alleles(const allele_counts& counts)
{
  size_t  allele_total = 0;

  allele_counts::const_iterator  a_iter     = counts.begin();
  allele_counts::const_iterator  a_end_iter = counts.end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    allele_total += a_iter->second;
  }
  
  return  allele_total;
}

// - Randomly choose an allele from among those with the greatest count.
//
allele_counts::iterator  
pool_pheno_seq_generator::iterator::greatest(allele_counts& counts)
{
  assert(! counts.empty());

  vector<allele_counts::iterator>  candidates;
  
  allele_counts::iterator  ac_iter     = counts.begin();
  allele_counts::iterator  ac_end_iter = counts.end();
  
  candidates.push_back(ac_iter);
  ++ac_iter;
  
  for(; ac_iter != ac_end_iter; ++ac_iter)
  {
    size_t  current_count = ac_iter->second;
    size_t  max_count = candidates.back()->second;
    
    if(current_count >= max_count)
    {
      if(current_count > max_count)
      {
        candidates.clear();
      }
      
      candidates.push_back(ac_iter);
    }
  }
  
  return  random_pick(candidates);
}

// - Randomly choose an allele from among those with the least count.
//
allele_counts::iterator  
pool_pheno_seq_generator::iterator::least(allele_counts& counts)
{
  assert(! counts.empty());

  vector<allele_counts::iterator>  candidates;
  
  allele_counts::iterator  ac_iter     = counts.begin();
  allele_counts::iterator  ac_end_iter = counts.end();
  
  candidates.push_back(ac_iter);
  ++ac_iter;
  
  for(; ac_iter != ac_end_iter; ++ac_iter)
  {
    size_t  current_count = ac_iter->second;
    size_t  min_count = candidates.back()->second;
    
    if(current_count <= min_count)
    {
      if(current_count < min_count)
      {
        candidates.clear();
      }
      
      candidates.push_back(ac_iter);
    }
  }
  
  return  random_pick(candidates);
}

// - Randomly choose an allele from among all alleles.
//
allele_counts::iterator  
pool_pheno_seq_generator::iterator::random(allele_counts& counts)
{
  assert(! counts.empty());

  vector<allele_counts::iterator>  candidates;
  
  allele_counts::iterator  ac_iter     = counts.begin();
  allele_counts::iterator  ac_end_iter = counts.end();
  for(; ac_iter != ac_end_iter; ++ac_iter)
  {
    candidates.push_back(ac_iter);
  }

  return  random_pick(candidates);
}

// - Randomly choose an allele from among those with greatest count if 
//   criteria is 'OVER' and from among those with least count if criteria
//   is 'UNDER'.
//
allele_counts::iterator  
pool_pheno_seq_generator::iterator::over_under(allele_counts& counts, bool over)
{
  if(over)
  {
    return  greatest(counts);
  }
  else
  {
    return  least(counts);
  }
}

const allele_counts::iterator&
pool_pheno_seq_generator::iterator::random_pick(const vector<allele_counts::iterator>& candidates)
{
  return  candidates[rand_src.uniform_integer(candidates.size())];
}
      
// - Value of instructions.pool_size_trait takes precedence over 
//   instructions.pool_size.
//
size_t
pool_pheno_seq_generator::iterator::get_pool_size() const
{
  const instructions&  user_instructions = my_generator->instr();
  size_t  trait_index = user_instructions.pool_size_trait;  
  size_t  pool_size = user_instructions.pool_size;

  if(trait_index != (size_t)(-1))
  {
    const FPED::FilteredMemberInfo&  member_info = (*my_member_iter)->info();
    double  trait_value = member_info.trait(trait_index);
    if(! SAGE::isnan(trait_value))
    {
      if(trait_value > 0)
      {
        size_t  integral_trait_value = round_2_integer(trait_value);
        pool_size = integral_trait_value;
      }
      else
      {
        cerrorstream&  errors = my_generator->errors();
        errors << priority(error) << "Size given for pool, " << pool_name() 
               << " is less than 1.  Setting pool size to " << pool_size  
               << " for this pool ..." << endl;        
      }
    }      
  }
  
  return  pool_size;
}

// - 'Statistician's method of rounding' from http://en.wikipedia.org/wiki/Rounding.  
//   Modified per conversation with Katrina Goddard so that 'round to even' rule 
//   is always invoked if the digit in the tenths place is a five.
//
size_t
pool_pheno_seq_generator::iterator::round_2_integer(double num)
{
  assert(num >= 0);
  
  size_t  whole_num = static_cast<size_t>(floor(num));
  size_t  tenths = static_cast<size_t>(floor(num * 10)) % 10;
  
  if((tenths > 5) || (tenths == 5 && (whole_num + 1) % 2 == 0))
  {
    ++whole_num;
  }
  
  return  whole_num; 
}


//============================================================================
// IMPLEMENTATION:  comb_generator
//============================================================================
//
ostream&  
pool_comb_generator::write_allele_group(ostream& out, const allele_group& group,
                                              const MLOCUS::inheritance_model* locus)
{
  out << "<";
  
  string  allele_name = "";
  
  allele_group::const_iterator  a_iter     = group.begin();
  allele_group::const_iterator  a_end_iter = group.end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    if(! allele_name.empty())
    {
      out << allele_name << ", ";
    }
  
    allele_name = locus->get_allele(*a_iter).name();
  }
 
  if(! allele_name.empty())
  {
    out << allele_name;
  }
  
  out << ">";
  
  return  out;
}

pool_comb_generator::allele_group
pool_comb_generator::counts_2_group(const allele_counts& counts)
{
  allele_group  group;

  allele_counts::const_iterator  c_iter     = counts.begin();
  allele_counts::const_iterator  c_end_iter = counts.end();
  for(; c_iter != c_end_iter; ++c_iter)
  {
    size_t  allele = c_iter->first;
    size_t  count  = c_iter->second; 
    for(size_t r = 0; r < count; ++r)
    {
      group.push_back(allele);
    }
  }

  return  group;  
}

allele_counts
pool_comb_generator::group_2_counts(const MLOCUS::inheritance_model* locus, const allele_group& group)
{
  allele_counts  counts;
  allele_counts::iterator  ac_end_iter = counts.end();
  
  MLOCUS::allele_iterator  a_iter     = locus->allele_begin();
  MLOCUS::allele_iterator  a_end_iter = locus->allele_end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    counts.insert(make_pair(a_iter->id(), 0));
  }  

  allele_group::const_iterator  g_iter     = group.begin();
  allele_group::const_iterator  g_end_iter = group.end();
  for(; g_iter != g_end_iter; ++g_iter)
  {
    assert(counts.find(*g_iter) != ac_end_iter);
    ++counts[*g_iter];
  }

  return  counts;  
}

void  
pool_comb_generator::write_comb(ostream& out, const pool_comb& comb, const pool_pheno_seq& seq)
{
  out << endl;

  pool_comb::const_iterator  c_iter     = comb.begin();
  pool_comb::const_iterator  c_end_iter = comb.end();
  for(; c_iter != c_end_iter; ++c_iter)
  {
    hap_seq  comb_elem = *c_iter;
    size_t  locus_count = comb_elem.size();
    
    assert(locus_count >= 1);
    assert(locus_count <= seq.size());
    
    string  current_allele = seq[0].first->get_allele(comb_elem[0]).name();
    for(size_t l = 1; l < locus_count; ++l)
    {
      out << current_allele << "-";
      current_allele = seq[l].first->get_allele(comb_elem[l]).name();
    }
    
    cout << current_allele << endl;
  }
}

void
pool_comb_generator::init_comb(pool_comb& comb, const allele_counts& counts)
{
  allele_counts::const_iterator  a_iter     = counts.begin();
  allele_counts::const_iterator  a_end_iter = counts.end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    size_t  allele_count = a_iter->second;
    for(size_t allele_copy = 0; allele_copy < allele_count; ++allele_copy)
    {
      hap_seq  new_seq(1, a_iter->first);
      comb.push_back(new_seq);
    }
  }  
}

// - Determine all of the haplotype combinations consistent
//   with sequence of marker phenotypes.
//
void
pool_comb_generator::generate(pool_pheno_seq_generator::iterator& seq_iter, vector<pool_comb>& combs)
{
  pool_pheno_seq  seq = *seq_iter;
  size_t  phenotype_count = seq_iter.phenotype_count();
  
  assert(! seq.empty());
  
  if(phenotype_count == 1)        // No missing information.
  {
    pool_comb  root_comb;
    init_comb(root_comb, seq[0].second);
    generate(root_comb, seq, combs);
    
    /*
    vector<pool_comb>::const_iterator  c_iter     = combs.begin();
    vector<pool_comb>::const_iterator  c_end_iter = combs.end();
    for(; c_iter != c_end_iter; ++c_iter)
    {
      write_comb(cout, *c_iter, *seq_iter);
      cout << endl;
    }
    */
  
  }
  else                                       // Missing information.
  {
    multi_pool_pheno_seq  expanded_pheno_seq; 
    expand_pheno_seq(seq_iter, expanded_pheno_seq);
    
#if 0
    multi_pool_pheno_seq::const_iterator  mpp_iter     = expanded_pheno_seq.begin();
    multi_pool_pheno_seq::const_iterator  mpp_end_iter = expanded_pheno_seq.end();
    for(; mpp_iter != mpp_end_iter; ++mpp_iter)
    {
      cout << "\nlocus " << mpp_iter->first->name() << endl;
      vector<allele_counts>::const_iterator  ac_iter     = mpp_iter->second.begin();
      vector<allele_counts>::const_iterator  ac_end_iter = mpp_iter->second.end();
      for(; ac_iter != ac_end_iter; ++ac_iter)
      {
        cout << "  " << make_pair(mpp_iter->first, *ac_iter) << endl;
      }
    }
#endif
    
    // - Generate all possible phenotype sequences and combinations for each.
    //
    vector<vector<pool_pheno_ptrs> >  sequences;
    
    multi_pool_pheno_seq::const_iterator  mp_iter     = expanded_pheno_seq.begin();
    multi_pool_pheno_seq::const_iterator  mp_end_iter = expanded_pheno_seq.end();
    assert(mp_iter != mp_end_iter);
    
    vector<allele_counts>::const_iterator  ac_iter     = mp_iter->second.begin();
    vector<allele_counts>::const_iterator  ac_end_iter = mp_iter->second.end();
    for(; ac_iter != ac_end_iter; ++ac_iter)
    {
      vector<pool_pheno_ptrs>  sequence;
      sequence.push_back(make_pair(mp_iter->first, &(*ac_iter)));
      
      multi_pool_pheno_seq::const_iterator  new_iter(mp_iter);
      enumerate_sequences(new_iter, mp_end_iter, sequence, sequences);
    }

    // - Determine haplotype combinations for each sequence.
    //
    size_t  sequence_count = sequences.size();
    for(size_t s = 0; s < sequence_count; ++s)
    {
      pool_pheno_seq  seq = make_dereferenced_sequence(sequences[s]);
      pool_comb  root_comb;
      init_comb(root_comb, seq[0].second);
      generate(root_comb, seq, combs);      
    }
  }
}

void
pool_comb_generator::enumerate_sequences(multi_pool_pheno_seq::const_iterator mp_iter,
                                         multi_pool_pheno_seq::const_iterator& mp_end_iter,
                                         const vector<pool_pheno_ptrs>& sequence,
                                         vector<vector<pool_pheno_ptrs> >& sequences)
{
  ++mp_iter;

  if(mp_iter == mp_end_iter)
  {
    sequences.push_back(sequence);
  }
  else
  {
    vector<allele_counts>::const_iterator  ac_iter     = mp_iter->second.begin();
    vector<allele_counts>::const_iterator  ac_end_iter = mp_iter->second.end();
    for(; ac_iter != ac_end_iter; ++ac_iter)
    {
      vector<pool_pheno_ptrs>  new_sequence(sequence);
      new_sequence.push_back(make_pair(mp_iter->first, &(*ac_iter)));
      
      multi_pool_pheno_seq::const_iterator  new_iter(mp_iter);      
      enumerate_sequences(new_iter, mp_end_iter, new_sequence, sequences);
    }              
  }
}

pool_pheno_seq
pool_comb_generator::make_dereferenced_sequence(const vector<pool_pheno_ptrs>& sequence)
{
  pool_pheno_seq  deref_seq;
   
  size_t  locus_count = sequence.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    deref_seq.push_back(make_pair(sequence[l].first, *(sequence[l].second)));
  }
  
  return  deref_seq;
}

void  
pool_comb_generator::expand_pheno_seq(pool_pheno_seq_generator::iterator& seq_iter, 
                                                    multi_pool_pheno_seq& expanded_seq)
{
  pool_pheno_seq::const_iterator  pheno_iter     = seq_iter->begin();
  pool_pheno_seq::const_iterator  pheno_end_iter = seq_iter->end();
  for(; pheno_iter != pheno_end_iter; ++pheno_iter)
  {
    vector<allele_counts>  genotypes;
    if(pheno_iter->second.empty())
    {
      enumerate_genotypes(pheno_iter->first, seq_iter.pool_size(), genotypes);        
    }
    else
    {
      genotypes.push_back(pheno_iter->second);
    }
    
    expanded_seq.push_back(make_pair(pheno_iter->first, genotypes));
  }
}

// - Determine all possible pool genotypes at the given locus.
//
void  
pool_comb_generator::enumerate_genotypes(const MLOCUS::inheritance_model* locus, size_t pool_size, 
                                                                   vector<allele_counts>& genotypes)
{
  MLOCUS::allele_iterator  a_iter     = locus->allele_begin();
  MLOCUS::allele_iterator  a_end_iter = locus->allele_end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    allele_group  genotype(1, a_iter->id());
    enumerate_genotypes(locus, pool_size, genotype, genotypes);
  }
}

// - Determine all possible pool genotypes at the given locus starting with the
//   partial genotype given.
//
void  
pool_comb_generator::enumerate_genotypes(const MLOCUS::inheritance_model* locus, size_t pool_size,
                                                                     const allele_group& genotype, 
                                                                    vector<allele_counts>& genotypes)
{
  if(genotype.size() == pool_size)
  {
    genotypes.push_back(group_2_counts(locus, genotype));
  }
  else
  {
    size_t  last_allele = genotype.back();
  
    MLOCUS::allele_iterator  a_iter     = locus->allele_begin();
    MLOCUS::allele_iterator  a_end_iter = locus->allele_end();
    for(; a_iter != a_end_iter; ++a_iter)
    {
      size_t  allele = a_iter->id();
      if(allele >= last_allele)
      {
        allele_group  new_genotype(genotype);
        new_genotype.push_back(allele);
        enumerate_genotypes(locus, pool_size, new_genotype, genotypes);
      }
    }  
  }
}

// - Given an incomplete haplotype combination, this function recursively generates
//   all possible combinations that are "descended" from it.
//
void  
pool_comb_generator::generate(pool_comb& comb, 
                              const pool_pheno_seq& seq, vector<pool_comb>& combs)
{
  size_t  locus_count = seq.size();
  
  assert(! comb.empty());
  
  if(comb.back().size() == locus_count)
  {
    combs.push_back(comb);
  }
  else
  {
    size_t  start = find_start(comb);
    //cout << "start " << start << endl;
    
    size_t  current_locus = comb[start].size();
    //cout << "current locus " << current_locus << endl;
    
    size_t  extension_size = get_extension_size(comb, start, locus_count);
    //cout << "extension size " << extension_size << endl;
    
    allele_counts  remaining_alleles = get_remaining_alleles(comb, seq, 
                                                             start, current_locus);
    //cout << "remaining alleles " << make_pair(seq[current_locus].first, remaining_alleles) << endl;
    
    set<allele_group>  groups;
    generate_allele_groups(remaining_alleles, extension_size, groups);
    
    /*
    set<allele_group>::const_iterator  g_iter     = groups.begin();
    set<allele_group>::const_iterator  g_end_iter = groups.end();
    for(; g_iter != g_end_iter; ++g_iter)
    {
      write_allele_group(cout, *g_iter, seq[current_locus].first);
      cout << " ";
    }
    
    cout << endl;
    */
    
    set<allele_group>::iterator  gr_iter     = groups.begin();
    set<allele_group>::iterator  gr_end_iter = groups.end();
    for(; gr_iter != gr_end_iter; ++gr_iter)
    {
      pool_comb  new_comb(comb);
      
      size_t  allele_count = gr_iter->size();
      for(size_t a = 0; a < allele_count; ++a)
      {
        size_t  hap = start + a;
        assert(hap < new_comb.size());
        new_comb[start + a].push_back((*gr_iter)[a]);
      }
      
      generate(new_comb, seq, combs);
    }
  }
}

// - Return index of first partial haplotype in the combination that
//   is shorter than the first one.  Return 0 if none exist.
//
size_t
pool_comb_generator::find_start(const pool_comb& comb)
{
  size_t  hap_count = comb.size();
  assert(hap_count > 0);

  size_t  start = 0;
  size_t  first_length = comb[0].size();
  
  for(size_t h = 1; h < hap_count; ++h)
  { 
    if(comb[h].size() < first_length)
    {
      start = h;
      break;
    }
  }
  
  return  start;
}

// - Determine number of partial haplotypes in the combination beginning with start
//   which are the same as start.  These will each be extended by one 
//   allele in the next branching.  A new combination(s) will be spawned for each
//   extension.  For example,
//
//   AAaa      AAaa      AAaa       AAaa
//         =>  bb   and  Bb    and  BB   if available alleles are bbBB.
//
//   In this case start is 0 (the first partial haplotype), the extension size is 2 and
//   3 new combinations are spawned.
//
size_t  
pool_comb_generator::get_extension_size(const pool_comb& comb, size_t start, size_t locus_count)
{
  assert(comb.back().size() < locus_count);
  
  size_t  haplotype_count = comb.size();
  
  size_t  extension_size = 1;
  for(size_t h = start + 1; h < haplotype_count; ++h)
  {
    if(comb[h] == comb[start])
    {
      ++extension_size;
    }
    else
    {
      break;
    }
  }
  
  return  extension_size;
}

// - Return alleles at the current locus that are "not yet used" in constructing
//   haplotypes for this combination.
//
allele_counts
pool_comb_generator::get_remaining_alleles(const pool_comb& comb, const pool_pheno_seq& seq,
                                                           size_t start, size_t current_locus)
{
  assert(start < comb.size());

  size_t  locus_count = seq.size();
  assert(0 < current_locus && current_locus < locus_count);
  
  allele_counts  remaining_alleles = seq[current_locus].second;
  
  for(size_t h = 0; h < start; ++h)
  {
    assert(comb[h].size() > current_locus);
    assert(remaining_alleles.count(comb[h][current_locus]) == 1);
    assert(remaining_alleles[comb[h][current_locus]] > 0);
    --remaining_alleles[comb[h][current_locus]];
  }
  
  return  remaining_alleles;
}               

// - Find all unique, unordered groups of the specified size using the elements of 
//   counts w/o replacement.
//
void 
pool_comb_generator::generate_allele_groups(const allele_counts& counts, size_t group_size, 
                                                                    set<allele_group>& results)
{
  allele_group  group = counts_2_group(counts);

  while(! group.empty())
  {
    size_t  allele = group.back();
    group.pop_back();
    
    allele_group  new_group(group);
    allele_group  starter(1, allele);
    
    generate_allele_groups(new_group, group_size, starter, results);
  }
}

// - Find all unique, unordered groups of the specified size using the elements of 
//   alleles w/o replacement using partial_group as a starting point.
//
void
pool_comb_generator::generate_allele_groups(allele_group& group, size_t group_size,
                                            allele_group& partial_group, set<allele_group>& results)
{
  if(partial_group.size() < group_size)
  {
    while(! group.empty())
    {
      size_t allele = group.back();
      group.pop_back();
      
      allele_group  new_group(group);
      allele_group  new_partial_group(partial_group);
      new_partial_group.push_back(allele);

      generate_allele_groups(new_group, group_size, new_partial_group, results);      
    }
  }
  else
  {
    sort(partial_group.begin(), partial_group.end());
    results.insert(partial_group);
  }
}                                            


//============================================================================
// IMPLEMENTATION:  pool_em_phenotype_map
//============================================================================
//
void
pool_em_phenotype_map::build()
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
  my_output_state.msg_interrupted = false;

  pool_pheno_seq_generator  pps_generator(my_member_pool, my_loci, my_errors,
                                         my_messages, my_instructions, my_output_state);
  pool_comb_generator  hc_generator;
  
  pool_pheno_seq_generator::iterator  iter = pps_generator.begin();
  pool_pheno_seq_generator::iterator  iter_end = pps_generator.end();
  for(; iter != iter_end; ++iter)
  {
    if(my_phenotype_directory.find(*iter) == my_phenotype_directory.end())
    {
      vector<pool_comb>  combs;
      
      try   // Generate possible combinations of haplotypes.
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
                  << "pool " << iter.pool_name() << ".  Skipping this pool ..."
                  << endl;
        continue;
      }
      
      try   // Make phenotypes and update global collection of haplotypes
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
                  << "the data.  Currently processing pool " << iter.pool_name()
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

set<member, member_order<member> >
pool_em_phenotype_map::members() const
{
  set<member, member_order<member> >  ordered_members;
  map<member, pool_pheno_seq>::const_iterator  m_iter     = my_member_directory.begin();
  map<member, pool_pheno_seq>::const_iterator  m_end_iter = my_member_directory.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    ordered_members.insert(m_iter->first);
  }
  
  return  ordered_members;
}

void
pool_em_phenotype_map::dump(ostream& out) const
{
  out << "Haplotypes:" << endl;
  my_haplotypes.dump(out, my_total_hap_count, my_loci);

  map<pool_pheno_seq, size_t>::const_iterator  p_iter     = my_phenotype_directory.begin();
  map<pool_pheno_seq, size_t>::const_iterator  p_end_iter = my_phenotype_directory.end();
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


