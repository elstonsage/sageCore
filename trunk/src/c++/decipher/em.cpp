//============================================================================
// File:      em.cpp
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

#include "decipher/em.h"
#include "numerics/print_util.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

std::string  
hap_seq_string(const hap_seq& seq, const locus_group& loci)
{
  ostringstream  out;
  write_hap_seq(out, seq, loci);
  
  return  out.str();
}

void
write_hap_seq(ostream& out, const hap_seq& seq, const locus_group& loci)
{
  size_t  seq_size = seq.size();
  for(size_t l = 0; l < seq_size; l++)
  {
    string  allele_str = "";
    if(seq[l] == MLOCUS::NPOS)
    {
      allele_str = loci[l].second->missing_allele_name();
    }
    else
    {
      allele_str = loci[l].second->get_allele(seq[l]).name();
    }
    
    out << allele_str;
    if(l < seq_size - 1)
    {
      out << "-";
    }
  }
}

void
write_markers(ostream& out, const locus_group& loci)
{
  out << "Markers in order:" << endl;
    
  out << "  ";
  size_t  locus_count = loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    out << loci[l].second->name() << " ";
  }
                     
  out << "\n\n" << endl;
}

void
dump_locus(ostream& out, const inheritance_model* locus)
{
  out << "BEGIN LOCUS DUMP" << endl;
  

  out << "missing allele name " << locus->missing_allele_name() << endl;
  out << "missing allele id " << locus->get_allele(locus->missing_allele_name()).id() << endl;
  
  /* This causes undefined behavior.  MLOCUS::NPOS is out-of-bounds for allele lookup.
  out << "MLOCUS::NPOS name " << locus->get_allele(MLOCUS::NPOS).name() << endl;
  */
  
  out << endl;
  
  allele_iterator  a_iter     = locus->allele_begin();
  allele_iterator  a_end_iter = locus->allele_end();
  for(; a_iter != a_end_iter; ++a_iter)
  {
    out << "allele id " << a_iter->id() << "  name " << a_iter->name() << endl;
  }
  
  out << "END LOCUS DUMP\n\n" << endl;
}

//============================================================================
// IMPLEMENTATION:  em_haplotype_map
//============================================================================
//
// - Revise haplotype counts to reflect new weights for haplotype 
//   combinations.  This is the 'maximization' step.
//
void
em_haplotype_map::update_counts(const base_em_phenotype_map& pm)
{
  reset();
  
  base_em_phenotype_map::const_iterator  p_iter = pm.begin();
  base_em_phenotype_map::const_iterator  p_end_iter = pm.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    if((! counts_initialized) || p_iter->is_ambiguous())
    {
      const vector<pair<em_phenotype::combination, double> >&  combs = p_iter->combinations();
      size_t  comb_count = combs.size();      
      for(size_t c = 0; c < comb_count; c++)
      {
        double  amount = p_iter->count() * combs[c].second;
        
        size_t  hap_count = combs[c].first.size();
        for(size_t h = 0; h < hap_count; h++)
        {
          (*this)[combs[c].first[h]].incr(amount);
        }
      }        
    }
  } 
}

void
em_haplotype_map::dump(ostream& out, size_t total_hap_count, const locus_group& loci) const
{
  map<size_t, em_haplotype>::const_iterator  h_iter = my_haplotypes.begin();
  map<size_t, em_haplotype>::const_iterator  h_end_iter = my_haplotypes.end();
  for(; h_iter != h_end_iter; ++h_iter)
  {
    write_hap_seq(out, h_iter->second.sequence(), loci);
    out << "   " << h_iter->second.new_freq(total_hap_count) << endl;
  }
}

//============================================================================
// IMPLEMENTATION:  em_phenotype
//============================================================================
//
em_phenotype::em_phenotype(em_haplotype_map* haplotypes,
                           const set<hap_seq_comb>& combinations)
      : my_count(1), ambiguous(combinations.size() > 1)
{
  assert(haplotypes != 0);

  // - Convert haplotype sequences into indices into the haplotype_map.
  //
  set<hap_seq_comb>::const_iterator  c_iter = combinations.begin();
  set<hap_seq_comb>::const_iterator  c_end_iter = combinations.end();
  for(; c_iter != c_end_iter; c_iter++)
  {
    my_combinations.push_back(make_pair(combination(), 0.0));
    hap_seq_comb::const_iterator  hs_iter = c_iter->begin();
    hap_seq_comb::const_iterator  hs_end_iter = c_iter->end();
    for(; hs_iter != hs_end_iter; ++hs_iter)
    {
      size_t  index = haplotypes->add_haplotype(*hs_iter);
      my_combinations.back().first.push_back(index);
    }
  }
  
  init_weights();
}

em_phenotype::em_phenotype(em_haplotype_map* haplotypes,
                           const vector<vector<hap_seq> >& combinations)
      : my_count(1), ambiguous(combinations.size() > 1)
{
  assert(haplotypes != 0);

  // - Convert haplotype sequences into indices into the haplotype_map.
  //
  vector<vector<hap_seq> >::const_iterator  c_iter = combinations.begin();
  vector<vector<hap_seq> >::const_iterator  c_end_iter = combinations.end();
  for(; c_iter != c_end_iter; c_iter++)
  {
    my_combinations.push_back(make_pair(combination(), 0.0));
    vector<hap_seq>::const_iterator  hs_iter = c_iter->begin();
    vector<hap_seq>::const_iterator  hs_end_iter = c_iter->end();
    for(; hs_iter != hs_end_iter; ++hs_iter)
    {
      size_t  index = haplotypes->add_haplotype(*hs_iter);
      my_combinations.back().first.push_back(index);
    }
  }
  
  init_weights();
}

// - Revise haplotype combination weights to reflect new haplotype counts. 
//   This is the 'expectation' step. Ito et. al.  AJHG 72:384-398, 2003.  
//   Equation (1).
//
void
em_phenotype::update_weights(size_t chromosome_count, em_haplotype_map& haplotypes)
{
  if(ambiguous)
  {
    double  prior_sum = 0.0;
  
    size_t  comb_count = my_combinations.size();
    for(size_t c = 0; c < comb_count; c++)
    {
      double  prior = calc_prior(my_combinations[c].first, chromosome_count, haplotypes);
      my_combinations[c].second = prior;
      prior_sum += prior;  
    }
    
    for(size_t c = 0; c < comb_count; c++)
    {
      my_combinations[c].second /= prior_sum;
    }
  }
}

double
em_phenotype::calc_prior(const combination& comb, size_t chromosome_count,
                         em_haplotype_map& haplotypes) const
{
  double  prior = factorial(comb.size());
  
  // - Iterate over distinct haplotypes in the combination.
  //
  combination_iterator  iter(comb);
  for(; ! iter.at_end(); iter++)
  {
    prior *= pow(haplotypes[(*iter).first].new_freq(chromosome_count),
                 static_cast<double>((*iter).second)) /
             factorial((*iter).second);
  }
  
  return  prior;
}

// - Assign weights randomly to haplotype combinations.
//
void
em_phenotype::init_weights()
{
  if(ambiguous)
  {
    size_t  total_weight = 0;
    size_t  comb_count = my_combinations.size();
    for(size_t c = 0; c < comb_count; c++)
    {
      size_t  weight = rand_src.uniform_integer(100) + 1;  // Don't want weight of 0.
      my_combinations[c].second = weight;
      total_weight += weight;
    }
    
    for(size_t c = 0; c < comb_count; c++)
    {
      my_combinations[c].second /= total_weight;
    }
  }
  else
  {
    my_combinations[0].second = 1.0;
  }
}

const set<comb_prob, greater<comb_prob> >&
em_phenotype::probabilities(const locus_group& loci, const em_haplotype_map& haplotypes) const
{
  my_probabilities.clear();
  vector<pair<combination, double> >::const_iterator  c_iter = my_combinations.begin();
  vector<pair<combination, double> >::const_iterator  c_end_iter = my_combinations.end();
  for(; c_iter != c_end_iter; ++c_iter)
  {
    my_probabilities.insert(comb_prob(combination_string(c_iter, loci, haplotypes), c_iter->second));
  } 

  return  my_probabilities;
}

// - Return combinations (diplotypes) formatted as a single string whose
//   constituent haplotypes are in alphabetical order.
//
const vector<string>&
em_phenotype::comb_strs(const locus_group& loci, const em_haplotype_map& haplotypes) const
{
  if(my_comb_strs.empty())
  {
    vector<pair<combination, double> >::const_iterator  c_iter = my_combinations.begin();
    vector<pair<combination, double> >::const_iterator  c_end_iter = my_combinations.end();
    for(; c_iter != c_end_iter; ++c_iter)
    {
      my_comb_strs.push_back(combination_string(c_iter, loci, haplotypes));
    }  
  }

  return  my_comb_strs;
}

string
em_phenotype::combination_string(const vector<pair<combination, double> >::const_iterator& c_iter,
                                 const locus_group& loci, const em_haplotype_map& haps) const
{
  multiset<string>  haplotypes;
  
  vector<size_t>::const_iterator  hs_iter = c_iter->first.begin();
  vector<size_t>::const_iterator  hs_end_iter = c_iter->first.end();
  for(; hs_iter != hs_end_iter; ++hs_iter)
  {
    haplotypes.insert(hap_seq_string(haps.index_to_hap_seq(*hs_iter), loci));
  }
  
  string  comb_str;
  multiset<string>::const_iterator  h_iter = haplotypes.begin();
  multiset<string>::const_iterator  h_end_iter = haplotypes.end();
  for(; h_iter != h_end_iter; ++h_iter)
  {
    comb_str += *h_iter;
    comb_str += "  ";
  }
  
  assert(! comb_str.empty());
  comb_str.resize(comb_str.size() - 1);   // Remove the final "/".
  
  return  comb_str;
}

void
em_phenotype::dump(ostream& out, const locus_group& loci, const em_haplotype_map& haplotypes) const
{
  out << "count  " << my_count << endl;
  out << "haplotype combinations" << endl;
  vector<pair<combination, double> >::const_iterator  c_iter = my_combinations.begin();
  vector<pair<combination, double> >::const_iterator  c_end_iter = my_combinations.end();
  for(; c_iter != c_end_iter; ++c_iter)
  {
    vector<size_t>::const_iterator  hs_iter = c_iter->first.begin();
    vector<size_t>::const_iterator  hs_end_iter = c_iter->first.end();
    for(; hs_iter != hs_end_iter; ++hs_iter)
    {
      write_hap_seq(out, haplotypes.index_to_hap_seq(*hs_iter), loci);
      out << endl;
    }
    
    out << "weight " << c_iter->second << "\n" << endl;
  }
}

//============================================================================
// IMPLEMENTATION:  em_phenotype::combination_iterator
//============================================================================
//
em_phenotype::combination_iterator::combination_iterator(const combination& comb)
{     
  vector<size_t>::const_iterator  iter = comb.begin();
  vector<size_t>::const_iterator  end_iter = comb.end();
  for(; iter != end_iter; ++iter)
  {
    if(previously_found(iter))
    {
      my_distinct_haplotypes[*iter]++;
    }
    else
    {
      my_distinct_haplotypes[*iter] = 1;
    }
  }
  
  my_internal_iterator = my_distinct_haplotypes.begin();
}

bool
em_phenotype::combination_iterator::previously_found(vector<size_t>::const_iterator iter) const
{
  bool  found = false;
  
  map<size_t, size_t>::const_iterator  h_iter = my_distinct_haplotypes.begin();
  map<size_t, size_t>::const_iterator  h_end_iter = my_distinct_haplotypes.end();
  for(; h_iter != h_end_iter; h_iter++)
  {
    if(*iter == h_iter->first)
    {
      found = true;
      break;
    }
  }
  
  return  found;
}

//============================================================================
// IMPLEMENTATION:  base_em_phenotype_map
//============================================================================
//
// - The EM algorithm for haplotype frequency estimation as described in
//   Excoffier and Slatkin, Mol. Biol. Evol. 12(5):921-927. 1995.
//
void
base_em_phenotype_map::maximize(double epsilon, size_t total_runs, ostream& dump_file,
                                bool dump, double cutoff, bool silent)
{
  if(maximized)
  {
    return;
  }
  
  assert(! empty());
  
  if(! silent)
  {
    my_messages << "      Maximizing likelihood ......................................" << flush;
  }
  
  if(dump)
  {
    init_dump(dump_file);
    dump_file << endl;
  }
  
  pop_freq_writer  freq_writer(dump_file, cutoff);
  for(size_t run = 0; run < total_runs; ++run)
  {
    init_weights();
    my_haplotypes.zero();
    my_haplotypes.update_counts(*this);   
    
    if(dump)
    {
      dump_file << "run " << run + 1 << endl;
      dump_file << "start point" << endl;
      freq_writer.set_frequencies(&current_frequencies());
      freq_writer.write(10, false);
      dump_file << endl;      
    }
    
    while(! my_haplotypes.converged(my_total_hap_count, epsilon))
    {
      update_weights();
      my_haplotypes.update_counts(*this);
    }
    
    double  current_ln_likelihood = ln_likelihood();
  
    if(dump)
    {
      dump_file << "end point" << endl;
      freq_writer.set_frequencies(&current_frequencies());
      freq_writer.set_ln_likelihood(current_ln_likelihood);
      freq_writer.write(10);
      dump_file << endl;
    }
    
    if(SAGE::isnan(my_max_ln_likelihood) || current_ln_likelihood > my_max_ln_likelihood)
    {
      my_max_ln_likelihood = current_ln_likelihood;
      my_best_haplotypes = my_haplotypes;
      my_best_phenotypes = my_phenotypes;    
    }    
  }
  
  maximized = true;
  
  if(! silent)
  {
    my_messages << "... done. " << endl;
  }
}

void 
base_em_phenotype_map::init_dump(ostream& dump_file) const
{
  string  pop_name = my_sub_pop_name.empty() ? "total population" : my_sub_pop_name;
  dump_file << "Maximizing " << pop_name << "...\n" << endl;
  pop_freq_writer::cutoff_note(dump_file);
  write_markers(dump_file, my_loci);     
}

const set<hap_freq, greater<hap_freq> >&
base_em_phenotype_map::final_frequencies() const
{
  assert(maximized);

  if(my_final_frequencies.empty())
  {
    const map<size_t, em_haplotype>& haps = my_best_haplotypes.haplotypes();
    set_frequencies(my_final_frequencies, haps);
  }

  return  my_final_frequencies;
}

const set<hap_freq, greater<hap_freq> >&
base_em_phenotype_map::current_frequencies() const
{
  my_current_frequencies.clear();
  const map<size_t, em_haplotype>& haps = my_haplotypes.haplotypes();
  set_frequencies(my_current_frequencies, haps);

  return  my_current_frequencies;
}

void
base_em_phenotype_map::set_frequencies(set<hap_freq, greater<hap_freq> >& frequencies,
                                       const map<size_t, em_haplotype>& haps          ) const
{
  map<size_t, em_haplotype>::const_iterator  iter = haps.begin();
  map<size_t, em_haplotype>::const_iterator  end_iter = haps.end();
  for(; iter != end_iter; ++iter)
  {
    string  haplotype = hap_seq_string(iter->second.sequence(), my_loci);
    double  frequency = iter->second.new_freq(my_total_hap_count);
    frequencies.insert(hap_freq(haplotype, frequency));
  } 
}

// - Equation (4), Excoffier and Slatkin, 1995 and Ito et. al.  AJHG 72:384-398, 2003.
//
double
base_em_phenotype_map::ln_likelihood()
{
  log_double  like(1.0); 
  const_iterator  p_iter = begin();
  const_iterator  p_end_iter = end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    log_double  sum(0);
    vector<pair<em_phenotype::combination, double> >::const_iterator  c_iter = p_iter->combinations().begin();
    vector<pair<em_phenotype::combination, double> >::const_iterator  c_end_iter = p_iter->combinations().end();
    for(; c_iter != c_end_iter; c_iter++)
    {
      sum += log_double(p_iter->calc_prior(c_iter->first, my_total_hap_count, my_haplotypes));
    }
    
    like *= sum.pow(static_cast<double>(p_iter->count()));
  }
  
  return  like.get_log();
}

void
base_em_phenotype_map::dump(ostream& out) const
{
  out << "Haplotypes:" << endl;
  my_haplotypes.dump(out, my_total_hap_count, my_loci);

  vector<em_phenotype>::const_iterator  p_iter = my_phenotypes.begin();
  vector<em_phenotype>::const_iterator  p_end_iter = my_phenotypes.end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    p_iter->dump(out, my_loci, my_haplotypes);
  }
}


//============================================================================
// IMPLEMENTATION: pop_freq_writer
//============================================================================
//
size_t
pop_freq_writer::max_seq_width(const set<hap_freq, greater<hap_freq> >& haps)
{
  size_t  max_width = 0;
  
  set<hap_freq, greater<hap_freq> >::const_iterator  h_iter = haps.begin();
  set<hap_freq, greater<hap_freq> >::const_iterator  h_end_iter = haps.end();
  for(; h_iter != h_end_iter; h_iter++)
  {
    size_t  width = h_iter->hap.size();
    if(width > max_width)
    {
      max_width = width;
    }
  }  
  
  return  max_width;
}

void
pop_freq_writer::write(size_t precision, bool show_likelihood)
{
  assert(my_frequencies != 0);

  my_out.setf(ios::showpoint);

  size_t  hap_width = max_seq_width(*my_frequencies);
  const size_t  SPACE = 10;
  const string  HAP_LABEL = "Haplotype";
  size_t  hap_label_width = HAP_LABEL.size();
  size_t  greater_width = hap_width > hap_label_width ? hap_width : hap_label_width;

  ostringstream  underline;
  underline << setw(hap_label_width) << setfill('-') << "" << setfill(' ');
  my_out << setw(hap_width) << HAP_LABEL << setw(SPACE) << "" << "Frequency\n"
         << setw(hap_width) << underline.str() << setw(SPACE) << "" << "---------" << endl;

  // - Body.
  //
  size_t FREQ_WIDTH = 9;

  size_t  lines_written = 0;
  double  total = 0;
  bool  first_written = false;
  set<hap_freq, greater<hap_freq> >::const_iterator  f_iter = my_frequencies->begin();
  set<hap_freq, greater<hap_freq> >::const_iterator  f_end_iter = my_frequencies->end();
  for(; f_iter != f_end_iter; f_iter++)
  {
    double  frequency = f_iter->freq;
    if((! first_written) || frequency >= my_cutoff)
    {
      total += frequency;
      my_out << setw(greater_width) << f_iter->hap;
      my_out << setw(SPACE) << "";
      if( frequency < 1.0e-4 )
        my_out << fp_scientific(frequency, FREQ_WIDTH - 4, FREQ_WIDTH - 6) << endl;
      else
        my_out << fp(frequency, FREQ_WIDTH, FREQ_WIDTH - 2) << endl;
      
      first_written = true;
      ++lines_written;
    }
  }
  
  // - Footer.
  //
  if(lines_written > 1)
  {
    ios::fmtflags  orig_flags = my_out.flags();
    my_out << setw(greater_width + SPACE) << "" << "---------\n"
           << setw(greater_width + SPACE) << "Total " << std::fixed << fp(total, FREQ_WIDTH, FREQ_WIDTH - 2) << endl;
               
    my_out.flags(orig_flags);
    my_out << showpoint;               
  }
  
  my_out << endl;
  
  if(show_likelihood)
  {
    assert(! SAGE::isnan(my_ln_likelihood));
    
    /*  
    my_out << "Ln likelihood    " << setprecision(precision) << my_ln_likelihood 
           << setprecision(6) << endl;
    */
    
    my_out << "Ln likelihood    " << double_to_string(my_ln_likelihood) << endl;
  }  
  
  my_out.unsetf(ios::showpoint);
}

}
}
