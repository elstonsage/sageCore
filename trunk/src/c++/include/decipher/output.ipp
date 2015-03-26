//============================================================================
// File:      output.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   11/15/4 - created.                                   djb
//                                                                          
// Notes:      
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


//============================================================================
// IMPLEMENTATION:  writer
//============================================================================
//
template<typename T> inline
writer<T>::writer(const output_state& os, ostream& summary, ostream& detail,
                  const instructions& instr, const vector<p_values>& pvs, const vector<ld_record>& lds)
    : my_output_state(os), my_summary(summary), my_detail(detail), my_instructions(instr),
      my_p_values(pvs), my_lds(lds)
{}

template<typename T> inline
writer<T>::writer(const output_state& os, ostream& summary, ostream& detail, 
               const instructions& instr, const vector<p_values>& pvs,
               const vector<const T*>& phenotype_maps, const vector<ld_record>& lds)
    : my_output_state(os), my_summary(summary), my_detail(detail), my_instructions(instr), 
      my_p_values(pvs), my_phenotype_maps(phenotype_maps), my_lds(lds)
{
  orig_summary_fmt = my_summary.flags();
  orig_detail_fmt  = my_detail.flags();
  
  my_summary.unsetf(ios::showpos | ios::fixed | ios::scientific);
  my_detail.unsetf(ios::showpos | ios::fixed | ios::scientific);
  
  my_summary << showpoint << left;
  my_detail << showpoint << left;
}

template<typename T> inline
writer<T>::~writer()
{
  my_summary.flags(orig_summary_fmt);
  my_detail.flags(orig_detail_fmt);
}

template<typename T> inline void
writer<T>::write(const locus_group& loci, const string& region, size_t block_number)
{
  write_prelims(my_summary, region, block_number);

  if(my_output_state.first_block && ! my_lds.empty())
  {
    if(! region.empty())
    {
      my_detail << "Region -- " << region << "\n" << endl;  
    }
  
    write_ld();
  }
  
  DECIPHER::write_markers(my_summary, loci);  
}

template<typename T> inline void
writer<T>::write_prelims(ostream& out, const string& region, size_t block_number)
{
  if(block_number)
  {
    // - 'dec' should not be necessary, but without it, block_number is written
    //    in hex notation.
    //
    my_summary << "Block " << noshowpos << dec << block_number << endl;
  }
}

template<typename T> inline void
writer<T>::write(const string& region, size_t block_number)
{
  write_prelims(my_summary, region, block_number);
  write_markers(my_summary);
  write_blank_lines(my_summary, 1);
  
  // - Haplotype frequency estimates.
  //
  if(my_instructions.pop_freq)
  {
    my_summary << "                         Haplotype Frequency Estimates\n\n" << endl;
    pop_freq_writer::cutoff_note(my_summary);
  
    size_t  map_count = my_phenotype_maps.size();
  
    if(map_count == 1)
    {
      write_freq_est(*(my_phenotype_maps[0]));
      write_blank_lines(my_summary, 2);
    }
    else
    {
      for(size_t m = 0; m < map_count; ++m)
      {
        write_freq_est(*(my_phenotype_maps[m]), my_phenotype_maps[m]->sub_pop_name());
        write_blank_lines(my_summary, 2);
      }
    }
    
    write_blank_lines(my_summary, 2);
  }
  
  // - Detailed results.
  //
  if(! region.empty())
  {
    my_detail << "Region -- " << region << "\n" << endl;  
  }

  if(my_output_state.first_block && ! my_lds.empty())
  {
    write_ld();
  }  
  
  if(block_number)
  {
    my_detail << "Block " << noshowpos << dec << block_number << endl;
  }
  
  write_markers(my_detail);
  write_blank_lines(my_detail, 1);    
  
  // - Most likely haplotype combinations.
  //
  if(my_instructions.most_likely_diplotypes)
  {
    write_most_likely();
    write_blank_lines(my_detail, 2);
  }
  
  // - Possible diplotypes.
  //
  if(my_instructions.all_possible_diplotypes || my_instructions.all_possible_diplotypes_table)
  {
    write_possible();
    write_blank_lines(my_detail, 2);
  }
  
  // - Likelihood ratio and permutation tests.
  //
  if(my_instructions.likelihood_ratio_test || my_instructions.compute_empirical_pvalue)
  {
    write_p_values();
    write_blank_lines(my_summary, 2);
  }
  
  // - All possible haplotypes.  This is a hidden feature.
  //
  if(my_instructions.all_possible_haplotypes)
  {
    my_detail << "\n\n                         All Possible Haplotypes\n\n" << endl;
  
    size_t  map_count = my_phenotype_maps.size();
  
    if(map_count == 1)
    {
      write_all_possible_haplotypes(*(my_phenotype_maps[0]));
      write_blank_lines(my_summary, 2);
    }
    else
    {
      for(size_t m = 0; m < map_count; ++m)
      {
        write_all_possible_haplotypes(*(my_phenotype_maps[m]), my_phenotype_maps[m]->sub_pop_name());
        write_blank_lines(my_summary, 2);
      }
    }
    
    write_blank_lines(my_summary, 2);
  }    
}

template<typename T> inline void
writer<T>::write_blank_lines(ostream& out, size_t count)
{
  size_t  adj_count = count - 1;

  for(size_t l = 0; l < adj_count; ++l)
  {
    out << "\n";
  }
  
  out << endl;
}

template<typename T> inline void
writer<T>::write_sub_pop_name(ostream& out, const string& sub_pop_name)
{
  if(! sub_pop_name.empty())
  {
    out << "Population:  " << sub_pop_name << endl;
  }
}

template<typename T> inline void
writer<T>::set_precision()
{
  double  epsilon = my_instructions.epsilon;
  
  size_t  precision = 0;
  while(epsilon < 1)
  {
    epsilon *= 10;
    precision += 1;
  }

  my_summary << setprecision(precision);  
  my_detail << setprecision(precision);
}

template<typename T> inline void
writer<T>::write_title(ostream& out)
{
  size_t  title_width = my_instructions.title.size();
  out << setfill('=') << setw(title_width) << "" << "\n"
      << my_instructions.title << "\n"
      << setfill('=') << setw(title_width) << "" << endl;
             
  out << setfill(' ');
}

template<typename T> inline void
writer<T>::write_markers(ostream& out) const
{
  assert(! my_phenotype_maps.empty());

  const locus_group&  loci = my_phenotype_maps[0]->loci();
  DECIPHER::write_markers(out, loci);  
}

template<typename T> inline void
writer<T>::write_ld() const
{
  size_t  max_width = find_max_marker_width();
  my_detail << "\n" << setfill(' ') << setw(max_width) << "" << "Linkage Disequilibrium\n" << endl;
  my_detail << setw(max_width) << left << "Marker" << "  " << setw(max_width) << left << "Marker" << "    LD" << endl;
  my_detail << setfill('-') << setw(max_width) << "" << "  " << setw(max_width) << "" << "    ---------" << endl;
  my_detail << setfill(' ');

  size_t  ld_count = my_lds.size();
  for(size_t L = 0; L < ld_count; ++L)
  {
    my_detail << setw(max_width) << left << my_lds[L].marker1 << "  "
              << setw(max_width)  << left << my_lds[L].marker2 <<  "    "
              << my_lds[L].ld << endl;
  }

  my_detail << "\n\n\n" << endl;       
}

template<typename T> inline size_t
writer<T>::find_max_marker_width() const
{
  size_t  max_width = 0;

  size_t  ld_count = my_lds.size();
  for(size_t L = 0; L < ld_count; ++L)
  {
    size_t  width = my_lds[L].marker1.size();
    if(width > max_width)
    {
      max_width = width; 
    }
    
    width = my_lds[L].marker2.size();
    if(width > max_width)
    {
      max_width = width; 
    }    
  }
  
  return  max_width;
}

template<typename T> inline void
writer<T>::init_output()
{
  // - Summary file.
  //
  write_title(my_summary);
  write_blank_lines(my_summary, 2);
  write_options();
  
  // - Detail file.
  //
  write_title(my_detail);
  write_blank_lines(my_detail, 2);
}

// - Write analysis options.
//
template<typename T> inline void
writer<T>::write_options()
{
  assert(my_instructions.partitions.size() == 2);
  string  sliding_window                = my_instructions.sliding_window ? "yes" : "no";
  string  four_gamete_rule              = my_instructions.four_gamete_rule ? "yes" : "no";
  string  ld_blocks                     = my_instructions.ld_blocks ? "yes" : "no";
  string  maf_filter                    = my_instructions.maf_filter ? "yes" : "no";  
  string  dump                          = my_instructions.dump ? "yes" : "no";
  string  subpopulations                = my_instructions.partitions[1].valid() ? "yes" : "no";
  string  analysis_unit                 = instructions::unit_2_string(my_instructions.analysis_unit);
  string  family_reps_specified         = my_instructions.family_rep != (size_t)(-1) ? "yes" : "no";
  string  haplotype_frequencies         = my_instructions.pop_freq ? "yes" : "no";
  string  all_possible_diplotypes       = my_instructions.all_possible_diplotypes ? "yes" : "no";
  string  all_possible_diplotypes_table = my_instructions.all_possible_diplotypes_table ? "yes" : "no";
  string  most_likely_diplotypes        = my_instructions.most_likely_diplotypes ? "yes" : "no";
  string  likelihood_ratio_test         = my_instructions.likelihood_ratio_test ? "yes" : "no";
  string  compute_empirical_pvalue      = my_instructions.compute_empirical_pvalue ? "yes" : "no";
  
  bool  show_em_parameters = my_instructions.pop_freq                 ||
                             my_instructions.most_likely_diplotypes   ||
                             my_instructions.likelihood_ratio_test    ||
                             my_instructions.compute_empirical_pvalue   ;
                             
  bool  show_dump          = my_instructions.pop_freq               ||
                             my_instructions.most_likely_diplotypes   ;

  // - General.
  //
  my_summary << "Options Selected" << "\n"
             << "================" << "\n";
             
  if(my_instructions.analysis_unit != instructions::POOL)
  {        
    size_t  region_count = my_instructions.regions.size();
    
    assert(region_count > 0);
    string  region = region_count > 1 ? "multiple" : my_instructions.regions.front().name;
   
    my_summary << "Haplotype region                                   " << region << endl;
  }
  
  // - Filters sub-block.
  //
  my_summary << "Minor allele frequency filter                      " << maf_filter << endl;
  if(maf_filter == "yes")
  {
    my_summary << "  Threshold                                        " 
               << noshowpoint << my_instructions.maf_threshold << showpoint << endl;
  }              
  
  
  // - Blocks sub-block.
  //
  my_summary << "Sliding window                                     " << sliding_window << endl;
  if(sliding_window == "yes")
  {
    my_summary << "  Window width                                     " << my_instructions.window_width << endl;
  }            
  
  my_summary << "Four gamete rule blocks                            " << four_gamete_rule << endl;
  if(four_gamete_rule == "yes")
  {
    my_summary << "  Threshold                                        " 
               << noshowpoint << my_instructions.fg_threshold << showpoint << endl;
  }

  my_summary << "LD blocks                                          " << ld_blocks << endl;
  if(ld_blocks == "yes")
  {
    my_summary << "  Threshold                                        " 
               << noshowpoint << my_instructions.ld_threshold << showpoint << endl;
  }              
             
  if(show_em_parameters)
  {           
    my_summary  << "EM algorithm convergence criterion                 " 
                << noshowpoint << my_instructions.epsilon << showpoint << "\n"
                << "Number of EM algorithm starting states             " 
                << my_instructions.starting_points << endl;
    if(show_dump)
    { 
      my_summary << "  Dump                                             " << dump << endl;
      
      if(my_instructions.dump)
      {
        my_summary << "    Cutoff                                         " 
                   << noshowpoint << my_instructions.dump_cutoff << showpoint << endl;
      }
    }
  }
  
  my_summary << endl;
    
  // - Data sub-block.
  // 
  my_summary << "Subpopulations specified                           " << subpopulations << "\n"
             << "Analysis unit                                      " << analysis_unit << endl;
             
  if(my_instructions.analysis_unit == instructions::FAMILY_REP)
  {
    my_summary << "  Family representatives specified                 " << family_reps_specified << endl;
  }             
  
  my_summary << endl;
  
  // - Tasks sub-block.
  //
  my_summary << "Estimate haplotype frequencies                     " << haplotype_frequencies << endl;
  
  if(my_instructions.pop_freq)
  {
    my_summary << noshowpoint;
    my_summary << "  Cutoff                                           " << my_instructions.freq_cutoff << endl;
    my_summary << showpoint;
  }
  
  my_summary << "Show all possible haplotype combinations table     " << all_possible_diplotypes_table << "\n"
             << "List most likely haplotype combinations            " << most_likely_diplotypes << endl;

  if(my_instructions.most_likely_diplotypes)
  {
    my_summary << noshowpoint;
    my_summary << "  Cutoff                                           " << my_instructions.likely_cutoff << endl;
    my_summary << showpoint;
  }             
             
  my_summary << "Do likelihood ratio test                           " << likelihood_ratio_test << "\n"
             << "Compute empirical p-value                          " << compute_empirical_pvalue << endl;
             
  if(my_instructions.compute_empirical_pvalue)
  {
    if(my_instructions.permutations != 0)
    {
      my_summary << "  Number of permutations                           " << my_instructions.permutations << endl;
    }
    else
    {
      my_summary << noshowpoint;
      my_summary << "  Maximum number of permutations                   " << my_instructions.max_permutations << "\n"
                 << "  Width                                            " << my_instructions.width            << "\n"
                 << "  Confidence                                       " << my_instructions.confidence       << endl;
      my_summary << showpoint;
    }
  }
}

// - Write population frequency estimates.
//
template<typename T> inline void
writer<T>::write_freq_est(const base_em_phenotype_map& phenotypes, const string& sub_pop_name)
{
  assert(phenotypes.is_maximized());
  
  pop_freq_writer  fw(my_summary, my_instructions.freq_cutoff);
  fw.set_frequencies(&(phenotypes.final_frequencies()));
  fw.set_ln_likelihood(phenotypes.max_ln_likelihood());
  
  // - Header.
  //
  if(! sub_pop_name.empty())
  {
    writer<T>::write_sub_pop_name(my_summary, sub_pop_name);
    writer<T>::write_blank_lines(my_summary, 1);
  }
  
  fw.write();
}

// - Write all possible haplotypes for the given population.
//
template<typename T> inline void
writer<T>::write_all_possible_haplotypes(const base_em_phenotype_map& phenotypes, const string& sub_pop_name)
{
  // - Header.
  //
  if(! sub_pop_name.empty())
  {
    writer<T>::write_sub_pop_name(my_detail, sub_pop_name);
    writer<T>::write_blank_lines(my_detail, 1);
  }
  
  my_detail << "Haplotypes\n"
            << "----------" << endl;
  
  const std::map<size_t, em_haplotype>&  haplotypes = phenotypes.haplotypes().haplotypes();
  
  set<string>  hap_strings;
  std::map<size_t, em_haplotype>::const_iterator  h_iter     = haplotypes.begin();
  std::map<size_t, em_haplotype>::const_iterator  h_end_iter = haplotypes.end();
  for(; h_iter != h_end_iter; ++h_iter)
  {
    hap_strings.insert(hap_seq_string(h_iter->second.sequence(), my_instructions.loci));
  }
  
  set<string>::const_iterator  hs_iter     = hap_strings.begin();
  set<string>::const_iterator  hs_end_iter = hap_strings.end();
  for(; hs_iter != hs_end_iter; ++hs_iter)
  {
    my_detail << *hs_iter << endl;
  }
}

// - Write most likely haplotype combination(s) for individuals.
//
template<typename T> inline void
writer<T>::write_most_likely()
{
  size_t  map_count = my_phenotype_maps.size();
  assert(map_count > 0);
  
  my_detail << "                      Most Likely Haplotype Combinations\n\n" << endl;
  my_detail << "Note:  Haplotype combinations listed have estimated probabilities greater than or equal\n"
            << "       to the cutoff or have the greatest probability estimate." << endl;  

  if(map_count == 1)
  {
    write_most_likely_pop(*(my_phenotype_maps[0]));
  }
  else
  {
    size_t  last_map_index = map_count - 1;
  
    for(size_t m = 0; m < map_count; ++m)
    {
      write_most_likely_pop(*(my_phenotype_maps[m]));
      if(m != last_map_index)
      {
        write_blank_lines(my_detail, 2);
      }
    }
  }
  
  write_blank_lines(my_detail, 2);
}

// - Write most likely haplotype combination(s) for individuals in a sub-population
//
template<typename T> inline void
writer<T>::write_most_likely_pop(const member_em_phenotype_map& phenotypes)
{
  assert(phenotypes.is_maximized());

  write_blank_lines(my_detail, 1);
  
  // - Header.
  //
  const string&  sub_pop_name = phenotypes.sub_pop_name();  
  if(! sub_pop_name.empty())
  {
    write_sub_pop_name(my_detail, sub_pop_name);
    write_blank_lines(my_detail, 1);
  }
    
  ml_widths  widths = write_most_likely_header(phenotypes);

  // - Body.
  //
  set<member, member_order<member> >  ordered_members = phenotypes.members();
  
  set<member, member_order<member> >::const_iterator  m_iter = ordered_members.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = ordered_members.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    const em_phenotype&  phenotype = phenotypes.final_phenotype(*m_iter);

    const set<comb_prob, greater<comb_prob> >&  probabilities = 
                  phenotype.probabilities(phenotypes.loci(), phenotypes.final_haplotypes());
    
    size_t  ped_spc = widths.ped + widths.space;
    size_t  mem_spc = widths.mem + widths.space;
    size_t  dip_spc = widths.dip + widths.space;
    size_t  ped_spc_mem_spc = ped_spc + mem_spc;
    set<comb_prob, greater<comb_prob> >::const_iterator  p_iter = probabilities.begin();
    set<comb_prob, greater<comb_prob> >::const_iterator  p_end_iter = probabilities.end();
    bool  first_prob = true;
    for(; p_iter != p_end_iter; ++p_iter)
    {
      if(! (first_prob || p_iter->prob >= my_instructions.likely_cutoff))
      {
        break;
      }
    
      bool  pool = my_instructions.analysis_unit == instructions::POOL;    
    
      if(first_prob)
      {
        if(pool)
        {
          string  pool_name = (*m_iter)->pedigree()->name() + ":" + (*m_iter)->name();
          my_detail << setw(ped_spc) << pool_name;
        }
        else
        {
          my_detail << setw(ped_spc) << (*m_iter)->pedigree()->name()
                    << setw(mem_spc) << (*m_iter)->name();
        }
                  
        my_detail << setw(dip_spc) << p_iter->comb << p_iter->prob << endl;
        first_prob = false;
      }
      else
      {
        my_detail << setw(pool ? ped_spc : ped_spc_mem_spc) << ""
                  << setw(dip_spc) << p_iter->comb << p_iter->prob << endl;      
      }
    }
  }
}

template<typename T> inline typename writer<T>::ml_widths
writer<T>::write_most_likely_header(const member_em_phenotype_map& phenotypes)
{
  const size_t  SPACE = 5;
  
  ml_widths  widths;
  widths.space = SPACE;
  calc_max_widths(phenotypes, widths);
  
  string  PED_LABEL = "";
  string  MEM_LABEL = "";
  
  if(my_instructions.analysis_unit == instructions::POOL)
  {
    PED_LABEL = "Pool";    // Treat as pool label.
    MEM_LABEL = "";        // Not needed.
  }
  else
  {
    PED_LABEL = "Pedigree";
    MEM_LABEL = "Member";
  }
  
  const string  DIP_LABEL = "Combination";
  const string  PRB_LABEL = "Probability";
  size_t  ped_label_width = PED_LABEL.size();
  size_t  mem_label_width = MEM_LABEL.size();
  size_t  dip_label_width = DIP_LABEL.size();
  size_t  prb_label_width = PRB_LABEL.size();
  
  if(my_instructions.analysis_unit == instructions::POOL)
  {
    size_t  pool_width = widths.ped + widths.mem + 1;
    widths.ped = ped_label_width > pool_width ? ped_label_width : pool_width;
    widths.mem = 0;  
  }
  else
  {
    widths.ped = ped_label_width > widths.ped ? ped_label_width : widths.ped;
    widths.mem = mem_label_width > widths.mem ? mem_label_width : widths.mem;
  }
    
  widths.dip = dip_label_width > widths.dip ? dip_label_width : widths.dip;  
  
  my_detail << setw(widths.ped) << PED_LABEL << setw(SPACE) << "";

  if(my_instructions.analysis_unit != instructions::POOL)
  {
    my_detail << setw(widths.mem) << MEM_LABEL << setw(SPACE) << "";
  }

  my_detail << setw(widths.dip) << DIP_LABEL << setw(SPACE) << ""
            << PRB_LABEL << endl;
             
  ostringstream  underline;
  underline << setw(ped_label_width) << setfill('-') << "" << setfill(' ');
  my_detail << setw(widths.ped + SPACE) << underline.str();
  
  if(my_instructions.analysis_unit != instructions::POOL)
  {
    underline.str("");  
    underline << setw(mem_label_width) << setfill('-') << "" << setfill(' ');
    my_detail << setw(widths.mem + SPACE) << underline.str();  
  }
  
  underline.str("");
  underline << setw(dip_label_width) << setfill('-') << "" << setfill(' ');
  my_detail << setw(widths.dip + SPACE) << underline.str();
  
  my_detail << setw(prb_label_width) << setfill('-') << "" << setfill(' ') << endl;
  
  return  widths;
}

// - Return maximum pedigree and member diplotype widths.
//
template<typename T> inline void
writer<T>::calc_max_widths(const member_em_phenotype_map& phenotypes, ml_widths& widths)
{
  // - Pedigree and member name widths.
  //
  set<member, member_order<member> >  members = phenotypes.members();
  
  size_t  max_ped_width = 0;
  size_t  max_mem_width = 0;
  set<member, member_order<member> >::const_iterator  m_iter = members.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = members.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    size_t  ped_width = (*m_iter)->pedigree()->name().size();
    size_t  mem_width = (*m_iter)->name().size();
    
    if(ped_width > max_ped_width)
    {
      max_ped_width = ped_width;
    }
    
    if(mem_width > max_mem_width)
    {
      max_mem_width = mem_width;
    }    
  }
  
  widths.ped = max_ped_width;
  widths.mem = max_mem_width;
  
  // - Diplotype widths.
  //
  size_t  max_width = 0;
  
  base_em_phenotype_map::const_iterator  ph_iter     = phenotypes.final_begin();
  base_em_phenotype_map::const_iterator  ph_end_iter = phenotypes.final_end();
  for(; ph_iter != ph_end_iter; ++ph_iter)
  {
    const set<comb_prob, greater<comb_prob> >& probs = ph_iter->probabilities(phenotypes.loci(),
                                                          phenotypes.final_haplotypes());
    set<comb_prob, greater<comb_prob> >::const_iterator  cp_iter     = probs.begin();
    set<comb_prob, greater<comb_prob> >::const_iterator  cp_end_iter = probs.end();
    for(; cp_iter != cp_end_iter; ++cp_iter)
    {
      size_t  width = cp_iter->comb.size();
      if(width > max_width)
      {
        max_width = width;
      }
    }
  }
  
  widths.dip = max_width;
}

// - Write all possible diplotypes for each individual.
//
template<typename T> inline void
writer<T>::write_possible()
{
  my_detail << "                      All Possible Haplotype Combinations" << endl;
  write_blank_lines(my_detail, 1);
  
  assert(! my_instructions.all_possible_diplotypes);
  
  /* Does not work.  Not used.
  if(my_instructions.all_possible_diplotypes)
  {
    write_possible_dispatch(my_phenotype_maps, NON_TABULAR);
    write_blank_lines(my_detail, 1);
  }
  */
  
  if(my_instructions.all_possible_diplotypes_table)
  {
    write_possible_dispatch(TABULAR);
  }
}

template<typename T> inline void  
writer<T>::write_possible_dispatch(format f)
{
  size_t  map_count = my_phenotype_maps.size();
  assert(map_count > 0);
  
  if(map_count == 1)
  {
    switch(f)
    {
      case TABULAR:
        write_possible_long_pop(my_phenotype_maps[0]);
        break;
      
      /* Not working.  Not used.  
      case NON_TABULAR:
        write_synteny_note();
        write_possible_short_pop(my_phenotype_maps[0]);
        break;
      */
        
      default:
        assert(false);
    }
  }
  else
  {
    size_t  last_map_index = map_count - 1;
  
    for(size_t m = 0; m < map_count; ++m)
    {
      switch(f)
      {
        case TABULAR:
          write_possible_long_pop(my_phenotype_maps[m]);
          break;
          
        /* Not working.  Not used.
        case NON_TABULAR:
          if(m == 0)
          {
            write_synteny_note();
          }
          
          write_possible_short_pop(my_phenotype_maps[m]);
          break;
        */
          
        default:
          assert(false);
      }
      
      if(m != last_map_index)
      {
        write_blank_lines(my_detail, 2);
      }
    }
  }  
}

template<typename T> inline void
writer<T>::write_synteny_note()
{
  my_detail << "Note:  In what follows, loci having the same number between alleles\n"
            << "       are always on the same respective haplotypes.  Absence of a\n" 
            << "       number indicates that phase is unknown.  Unknown alleles are\n" 
            << "       designated with an allele missing value code." << endl;
}

// - Write possible haplotype combinations (using long format) for individuals
//   in a sub-population.
//
template<typename T> inline void
writer<T>::write_possible_long_pop(const member_em_phenotype_map* phenotypes)
{
  write_sub_pop_name(my_detail, phenotypes->sub_pop_name());
  write_blank_lines(my_detail, 1);
  
  map<string, map<member, bool> >  table;
  
  // - Build the diplotypes x members table.
  //
  set<member, member_order<member> >  member_set = phenotypes->members();    
  
  base_em_phenotype_map::const_iterator  p_iter = phenotypes->begin();
  base_em_phenotype_map::const_iterator  p_end_iter = phenotypes->end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    map<member, bool>  member_map;

    set<member, member_order<member> >::const_iterator  m_iter = member_set.begin();
    set<member, member_order<member> >::const_iterator  m_end_iter = member_set.end();
    for(; m_iter != m_end_iter; ++m_iter)
    {
      if(&(*phenotypes)[*m_iter] == &(*p_iter))
      {
        member_map[*m_iter] = true;
      }
    }
    
    const vector<string>& combinations = p_iter->comb_strs(phenotypes->loci(), phenotypes->haplotypes());
    
    size_t  cs_count = combinations.size();
    for(size_t  cs = 0; cs < cs_count; ++cs)
    {
      string   combination = combinations[cs];
      map<string, map<member, bool> >::iterator  comb_end_iter = table.end();
      map<string, map<member, bool> >::iterator  comb_result   = table.find(combination);
      if(comb_result == comb_end_iter)
      {
        table.insert(make_pair(combination, member_map));
      }
      else
      {
        map<member, bool>::const_iterator  m_iter     = member_map.begin();
        map<member, bool>::const_iterator  m_end_iter = member_map.end();
        for(; m_iter != m_end_iter; ++m_iter)
        {
          table[combination][(*m_iter).first] = true;
        }
      }
    }
  }
  
  // - Display the diplotypes x members table.
  //
  size_t  dip_col_width = diplotype_col_width(table);
  size_t  mem_col_width = member_col_width(member_set);
  write_possible_table_header(member_set, dip_col_width, mem_col_width);
  
  size_t  left_width = static_cast<size_t>(ceil(static_cast<double>(mem_col_width - 1) / 2.0));
  size_t  right_width = mem_col_width - left_width;
  
  map<string, map<member, bool> >::const_iterator  t_iter = table.begin();
  map<string, map<member, bool> >::const_iterator  t_end_iter = table.end();
  for(; t_iter != t_end_iter; ++t_iter)
  {
    my_detail << setw(dip_col_width) << t_iter->first;
    
    my_detail << right;
    
    map<member, bool>::const_iterator  member_not_found = t_iter->second.end();

    set<member, member_order<member> >::const_iterator  m_iter     = member_set.begin();
    set<member, member_order<member> >::const_iterator  m_end_iter = member_set.end();
    for(; m_iter != m_end_iter; ++m_iter)
    {
      my_detail << setw(left_width);
    
      if(t_iter->second.find(*m_iter) == member_not_found)
      {
        my_detail << "-";
      }
      else
      {
        my_detail << "x";
      }
      
      my_detail << setw(right_width) << "";
    }
    
    my_detail << left << endl;
  }
}

template<typename T> inline size_t
writer<T>::diplotype_col_width(const map<string, map<member, bool> >& table) const
{
  size_t  max_col_width = 0;
  
  map<string, map<member, bool> >::const_iterator t_iter = table.begin();
  map<string, map<member, bool> >::const_iterator t_end_iter = table.end();
  for(; t_iter != t_end_iter; ++t_iter)
  {
    size_t  col_width = t_iter->first.size() + 1;
  
    if(col_width > max_col_width)
    {
      max_col_width = col_width;
    }
  }
  
  return  max_col_width;
}

template<typename T> inline size_t
writer<T>::member_col_width(const set<member, member_order<member> >& members) const
{
  size_t  max_col_width = 0;
  
  set<member, member_order<member> >::const_iterator  m_iter = members.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = members.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    size_t  col_label_size = ((*m_iter)->pedigree()->name()).size() + ((*m_iter)->name()).size() + 2;
    if(col_label_size > max_col_width)
    {
      max_col_width = col_label_size;
    }
  }
  
  return  max_col_width;
}

template<typename T> inline void
writer<T>::write_possible_table_header(const set<member, member_order<member> >& members, size_t offset, size_t col_width)
{
  my_detail << setw(offset + 10) << "";
  
  if(my_instructions.analysis_unit == instructions::POOL)
  {
    my_detail << "Pool" << endl;
  }
  else
  {
    my_detail << "Pedigree:Member" << endl;
  }
  
  my_detail << setw(offset) << "";

  set<member, member_order<member> >::const_iterator  m_iter = members.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = members.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    my_detail << setw(col_width) << (*m_iter)->pedigree()->name() + ":" + (*m_iter)->name();
  }
  
  my_detail << endl;
}

/* Does not work.  Not used.
// - Write possible haplotype combinations (using short, i.e. vertical, format) for individuals
//   in a sub-population.
//
template<typename T> inline void  
writer<T>::write_possible_short_pop(const member_em_phenotype_map* phenotypes)
{
  const size_t  SPACER = 5;
  write_sub_pop_name(my_detail, phenotypes->sub_pop_name());
  write_blank_lines(my_detail, 1);
  
  size_t  member_width = max_member_width(phenotypes);
  member_width += SPACER;
  
  set<member, member_order<member> >  member_set = phenotypes->members();
  set<member, member_order<member> >::const_iterator  m_iter = member_set.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = member_set.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    my_detail << setw(member_width) << member_name(*m_iter);
    write_member_possibilities((*phenotypes)[*m_iter], phenotypes->haplotypes(),
                                phenotypes->loci(), member_width);
    write_blank_lines(my_detail, 1);
  }
}

template<typename T> inline size_t
writer<T>::max_member_width(const member_em_phenotype_map* phenotypes)
{
  size_t  max_width = 0;

  set<member, member_order<member> >  member_set = phenotypes->members();
  set<member, member_order<member> >::const_iterator  m_iter     = member_set.begin();
  set<member, member_order<member> >::const_iterator  m_end_iter = member_set.end();
  for(; m_iter != m_end_iter; ++m_iter)
  {
    size_t  width = member_name(*m_iter).size();
    if(width > max_width)
    {
      max_width = width;    
    }
  }
  
  return  max_width;  
}

template<typename T> inline string
writer<T>::member_name(member m)
{
  return "pedigree " + m->pedigree()->name() + ", " +
          "member " + m->name();
}

// - Construct haplotype combinations (DIPLOTYPES ONLY) for the given phenotype
//   in a compact format using an indicator variable to show phase
//   information.
//
template<typename T> inline void
writer<T>::write_member_possibilities(const em_phenotype& phenotype, 
                                   const em_haplotype_map& haplotypes,
                                   const locus_group& loci,
                                   size_t member_width)
{
  vector<pair<hap_seq, hap_seq> >  new_combs;
  build_new_combs(phenotype, haplotypes, loci, new_combs);
  
  vector<vector<string> >  lines(3);
  
  // - First combination.
  //
  assert(! new_combs.empty());

  const hap_seq&  hap_seq_00 = new_combs[0].first;
  const hap_seq&  hap_seq_01 = new_combs[0].second;
  
  size_t  locus_count = hap_seq_00.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    const string&  missing = loci[l].second->missing_allele_name();
    lines[0].push_back(hap_seq_00[l] == MLOCUS::NPOS ? missing : loci[l].second->get_allele(hap_seq_00[l]).name());
    lines[2].push_back(hap_seq_01[l] == MLOCUS::NPOS ? missing : loci[l].second->get_allele(hap_seq_01[l]).name());
  }
  
  // - Calculate phase indicators.
  //
  vector<size_t>  indicators;
  
  indicators.push_back(1);
  lines[1].push_back(long2str(1));
  for(size_t l = 1; l < locus_count; ++l)
  {
    indicators.push_back(calc_indicator(indicators, new_combs));
    lines[1].push_back(long2str(indicators.back())); 
  }
  
  cull_phase_unknown(lines[1]);
  write_possible_lines(lines, member_width);
}

// - Look at all pairs of alleles across diplotypes.  Where more than one unordered
//   pair exists, substitute MLOCUS::NPOS at one or both alleles in the haplotypes
//   which make up the diplotypes.
//
template<typename T> inline void
writer<T>::build_new_combs(const em_phenotype& phenotype,
                        const em_haplotype_map& haplotypes,
                        const locus_group& loci,
                        vector<pair<hap_seq, hap_seq> >& new_combs)
{
  const vector<pair<em_phenotype::combination, double> >&  combs = phenotype.combinations();
  size_t  i_model_count = loci.size();

  // - Get sequences.
  //
  size_t  combs_count = combs.size();
  for(size_t c = 0; c < combs_count; ++c)
  {
    assert(combs[c].first.size() == 2);           // Diplotype.
    new_combs.push_back(pair<vector<size_t>, vector<size_t> >());

    new_combs[c].first  = haplotypes.index_to_hap_seq(combs[c].first[0]);
    new_combs[c].second = haplotypes.index_to_hap_seq(combs[c].first[1]);
  }
  
  
  // - Dump new combinations.
  //
  for(size_t c = 0; c < combs_count; ++c)
  {
    size_t  locus_count = new_combs[0].first.size();
    assert(locus_count == i_model_count);
    for(size_t l = 0; l < locus_count; ++l)
    {
      cout << new_combs[c].first[l] << "-";
    }
    
    cout << endl;
    
    for(size_t l = 0; l < locus_count; ++l)
    {
      cout << new_combs[c].second[l] << "-";
    }
    
    cout << "\n" << endl;    
  }
  
  
  // - Find multiple allele pairs.
  //
  cout << endl;
  cout << boolalpha;
  for(size_t l = 0; l < i_model_count; ++l)
  {
    cout << "locus " << l << " has multiple allele pairs, " <<  multiple_allele_pairs(l, new_combs) << endl;
    
    if(multiple_allele_pairs(l, new_combs))
    {
      insert_missing_index(l, new_combs);
    }
  }
  
  //cout << boolalpha;
  
  
  // - Dump new combinations.
  //
  //cout << "new combinations" << endl;
  //for(size_t c = 0; c < combs_count; ++c)
  //{
  //  size_t  locus_count = new_combs[0].first.size();
  //  assert(locus_count == i_model_count);
  //  for(size_t l = 0; l < locus_count; ++l)
  //  {
  //    cout << new_combs[c].first[l] << "-";
  //  }
  //  
  //  cout << endl;
  //
  //  for(size_t l = 0; l < locus_count; ++l)
  //  {
  //    cout << new_combs[c].second[l] << "-";
  //  }
  //  
  //  cout << "\n" << endl;    
  //} 
   
}

// - Among the combinations, are there multiple unordered pairs of alleles?
//
template<typename T> inline bool
writer<T>::multiple_allele_pairs(size_t l,
                              const vector<pair<hap_seq, hap_seq> >& new_combs)
{
  bool  multiple_pairs = false;
  
  size_t  comb_count = new_combs.size();
  assert(comb_count > 0);
    
  pair<size_t, size_t>  alleles_l0 = make_pair(new_combs[0].first[l], new_combs[0].second[l]);
  for(size_t c = 1; c < comb_count && ! multiple_pairs; ++c)
  {
    pair<size_t, size_t>  alleles_lc = make_pair(new_combs[c].first[l], new_combs[c].second[l]);
    if(! same_alleles(alleles_l0, alleles_lc))
    {
      multiple_pairs = true;
    }
  }
  
  return  multiple_pairs;    
}

template<typename T> inline bool
writer<T>::same_alleles(const pair<size_t, size_t>& p1, const pair<size_t, size_t>& p2)
{
  return  p1.first == p2.first  && p1.second == p2.second ||
          p1.first == p2.second && p1.second == p2.first    ;
}

// - Substitute MLOCUS::NPOS for ambiguous alleles.
//
template<typename T> inline void
writer<T>::insert_missing_index(size_t l, vector<pair<hap_seq, hap_seq> >& new_combs)
{
  bool  replace_first = false;
  bool  replace_second = false;

  size_t  comb_count = new_combs.size();
  assert(comb_count > 0);
  
  size_t  allele_00l = new_combs[0].first[l];
  size_t  allele_01l = new_combs[0].second[l];
  for(size_t c = 1; c < comb_count && ! (replace_first && replace_second); ++c)
  {
    if(allele_00l != new_combs[c].first[l])
    {
      replace_first = true;
    } 
    
    if(allele_01l != new_combs[c].second[l])
    {
      replace_second = true;
    }     
  }
  
  for(size_t c = 0; c < comb_count; ++c)
  {
    if(replace_first)
    {
      new_combs[c].first[l] = MLOCUS::NPOS;
    }
    
    if(replace_second)
    {
      new_combs[c].second[l] = MLOCUS::NPOS;    
    }
  }
}

// - If an indicator appears only once, replace it with a blank.
//
template<typename T> inline void
writer<T>::cull_phase_unknown(vector<string>& indicators)
{
  vector<string>::iterator  begin_iter = indicators.begin();
  vector<string>::iterator  end_iter   = indicators.end();
  vector<string>::iterator  iter       = begin_iter;
  for(; iter != end_iter; ++iter)
  {
    if(count(begin_iter, end_iter, *iter) == 1)
    {
      *iter = " ";
    }
  }
}

// - Determine indicator for next locus.  Alleles at loci with same indicator
//   are always syntenic.
// 
template<typename T> inline size_t
writer<T>::calc_indicator(const vector<size_t>& indicators,
                       const vector<pair<hap_seq, hap_seq> >& new_combs)
{
  size_t  comb_count = new_combs.size();
  size_t  cl = indicators.size();                // Current locus.
  size_t  max_ind = *max_element(indicators.begin(), indicators.end());
  
  for(size_t ind = 1; ind <= max_ind; ++ind)     // Try to match an existing indicator.
  {
    bool  ind_good = true;
    for(size_t l = 0; l < cl && ind_good; ++l)   // Check all loci with the current indicator.
    {
      if(indicators[l] == ind)
      {
        const hap_seq&  hap_seq_00 = new_combs[0].first;
        const hap_seq&  hap_seq_01 = new_combs[0].second;        
        two_locus_diplotype  dip0(l, cl, hap_seq_00[l], hap_seq_00[cl], 
                                         hap_seq_01[l], hap_seq_01[cl] );
        
        // - Do all combinations indicate current locus syntenic with locus, l?
        //
        size_t  tld_count = 1;
        for(size_t c = 1; c < comb_count && tld_count == 1; ++c)
        {
          const hap_seq&  hap_seq_c0 = new_combs[c].first;
          const hap_seq&  hap_seq_c1 = new_combs[c].second;
          two_locus_diplotype  dip(l, cl, hap_seq_c0[l], hap_seq_c0[cl], 
                                          hap_seq_c1[l], hap_seq_c1[cl] );      
          
          if(dip != dip0)
          {
            ++tld_count;
          }
        }
        
        if(tld_count > 1)   // Not this indicator.  Try next one.
        {
          ind_good = false;
        } 
      }
    }
    
    if(ind_good)
    {
      return  ind;
    }
  } 
    
  return  ++max_ind;
}

// - Write so that haplotypes run vertically, not horizontaly as implied by
//   structure of the lines matrix.
//
template<typename T> inline void
writer<T>::write_possible_lines(const vector<vector<string> >& lines, size_t member_width)
{
  size_t  col_count = lines.size();
  size_t  row_count = lines[0].size();

  vector<size_t>  col_widths;
  calc_column_widths(lines, col_widths);
  
  for(size_t row = 0; row < row_count; ++row)
  {
    if(row != 0)
    {
      my_detail << setfill(' ') << setw(member_width) << ""; 
    }
  
    for(size_t col = 0; col < col_count; ++col)
    {
      my_detail << setw(col_widths[col]) << lines[col][row] << " ";
    }
    
    my_detail << endl;
  }
}

template<typename T> inline void
writer<T>::calc_column_widths(const vector<vector<string> >& lines, vector<size_t>& col_widths)
{
  size_t  col_count = lines.size();
  size_t  row_count = lines[0].size();
  
  for(size_t col = 0; col < col_count; ++col)
  {
    size_t  max_width = 0;
    for(size_t row = 0; row < row_count; ++row)
    {
      size_t  width = lines[col][row].size();
      if(width > max_width)
      {
        max_width = width;
      }
    }
    
    col_widths.push_back(max_width);
  }  
}
*/

// - Likelihood ratio and permutation tests.
//
template<typename T> inline void
writer<T>::write_p_values()
{
  size_t p_value_count = my_p_values.size();
  if(p_value_count)
  {
    my_summary << "                         Likelihood Ratio Tests" << endl;
    write_blank_lines(my_summary, 1);  
  
    for(size_t pv = 0; pv < p_value_count; ++pv)
    {
      my_summary << "Population:  " << my_p_values[pv].outer_sub_pop_name << endl;
      
      if(my_instructions.likelihood_ratio_test)
      {
        string  a_value   = double_to_string(my_p_values[pv].asymptotic);
        string  c_ln_like = double_to_string(my_p_values[pv].composite_ln_like);        
        string  w_ln_like = double_to_string(my_p_values[pv].whole_ln_like);        
        string  t_stat    = double_to_string(my_p_values[pv].test_statistic);
        
        my_summary << "Ln likelihood (haplotype frequencies estimated separately for each constituent subpopulation)  "
                   << c_ln_like << endl;
        my_summary << "Ln likelihood (one set of haplotype frequencies estimated for whole population)  "
                   << w_ln_like << endl;
        my_summary << "Test statistic  "<< t_stat << endl;
        my_summary << "Degrees of freedom  " << my_p_values[pv].degrees_of_freedom << endl;                   
        my_summary << "Asymptotic p-value  "  << a_value << endl;
      }
    
      if(my_instructions.compute_empirical_pvalue)
      {
        string  e_value = SAGE::isnan(my_p_values[pv].empirical) ? "  ---" : doub2str(my_p_values[pv].empirical);  
        my_summary << "Empirical p-value   "  << e_value 
                   << "  (" << my_p_values[pv].permutations << " permutations)" << endl;
      }
      
      write_blank_lines(my_summary, 2);
    }

    my_summary << "\nNote:  '---' means the value could not be calculated." << endl;      
  }
}

//============================================================================
// IMPLEMENTATION:  writer<T>::ml_widths
//============================================================================
//
template<typename T> inline
writer<T>::ml_widths::ml_widths()
{}

/* Client code does not work, so this code is not used.
//============================================================================
// IMPLEMENTATION:  writer::two_locus_diplotype
//============================================================================
//
inline
writer::two_locus_diplotype::two_locus_diplotype(size_t l1, size_t l2, 
                                                 size_t h1l1, size_t h1l2, size_t h2l1, size_t h2l2)
      : locus1(l1), locus2(l2), 
        hap1_locus1(h1l1), hap1_locus2(h1l2), hap2_locus1(h2l1), hap2_locus2(h2l2)
{}

inline bool
writer::two_locus_diplotype::operator ==(const two_locus_diplotype& other)
{
  return  locus1 == other.locus1 && locus2 == other.locus2 && 
          ((hap1_locus1 == other.hap1_locus1 && hap1_locus2 == other.hap1_locus2) &&
           (hap2_locus1 == other.hap2_locus1 && hap2_locus2 == other.hap2_locus2)   ) ||
          ((hap1_locus1 == other.hap2_locus1 && hap1_locus2 == other.hap2_locus2) &&
           (hap2_locus1 == other.hap1_locus1 && hap2_locus2 == other.hap1_locus2)   )   ;
}

inline bool
writer::two_locus_diplotype::operator !=(const two_locus_diplotype& other)
{
  return  ! operator ==(other);
}

// - Write so haplotypes are displayed vertically.
//
inline void
writer::write_two_locus_diplotype(ostream& out, const two_locus_diplotype& dip,
                                    const locus_group& loci)
{
  const string&  missing1 = loci[dip.locus1].second->missing_allele_name();
  const string&  missing2 = loci[dip.locus2].second->missing_allele_name();

  out << dip.locus1 << "  " 
      << (dip.hap1_locus1 == MLOCUS::NPOS ? missing1 : loci[dip.locus1].second->get_allele(dip.hap1_locus1).name()) << " "
      << (dip.hap2_locus1 == MLOCUS::NPOS ? missing1 : loci[dip.locus1].second->get_allele(dip.hap2_locus1).name()) << endl;
  out << dip.locus2 << "  " 
      << (dip.hap1_locus2 == MLOCUS::NPOS ? missing2 : loci[dip.locus2].second->get_allele(dip.hap1_locus2).name()) << " "
      << (dip.hap2_locus2 == MLOCUS::NPOS ? missing2 : loci[dip.locus2].second->get_allele(dip.hap2_locus2).name()) << endl;  
}
*/ 


