//============================================================================
// File:      analysis.cpp
//
// Author:    Dan Baechle
//
// History:   4/8/4 - created.                         djb
//
// Notes:     Implementation of analysis class.
//
//
// Copyright (c) 2004 R.C. Elston
// All Rights Reserved
//============================================================================

#include "decipher/analysis.h"

namespace SAGE
{

namespace DECIPHER
{

const double  LD_ERROR_TOLERANCE = .05;
const double  ALLELE_FREQUENCY_TOLERANCE = .20;

void
write_blocks(ostream& out, const block_vector& bv)
{
  out << endl;

  size_t  block_count = bv.size();
  for(size_t b = 0; b < block_count; ++b)
  {
    out << "(" << bv[b].first << ", " << bv[b].second << ")" << endl;
  }
}


//============================================================================
// IMPLEMENTATION:  analysis
//============================================================================
//
// - Perform calculations per user instructions and write results.
//
void
analysis::build_filtered_mped()
{
  FPED::MPFilterer::add_multipedigree_filtered_by_members(my_filtered_mped, my_mped, FPED::always_keep());
  my_filtered_mped.construct();
}

void
analysis::analyze()
{
  assert(my_instructions.valid);
  set_test_seed();

  my_messages << "\n"
              << "  Performing analysis: " << my_instructions.title << double_line << endl;

  write_user_options();

  if(nothing_to_do())
  {
    my_messages << "    Nothing to do for this analysis ...\n" << double_line << endl;
    return;
  }

  my_summary_file << "\n\nResults\n"
                  << "=======\n\n" << endl;

  if(my_instructions.analysis_unit != instructions::POOL)
  {
    do_non_pool_analysis();
  }
  else
  {
    analyze_block(my_current_region);
  }

  my_messages << double_line << endl;
}

// - For test purposes.  Undocumented feature.
//
void
analysis::set_test_seed()
{
  if(my_instructions.seed)
  {
    rand_src.reseed(my_instructions.seed);
  }
}

bool
analysis::maximization_needed() const
{
  return   my_instructions.pop_freq                  ||
           my_instructions.most_likely_diplotypes    ||
           my_instructions.likelihood_ratio_test     ||
           my_instructions.compute_empirical_pvalue    ;
}

bool
analysis::nothing_to_do() const
{
  return   ! (maximization_needed()                         ||
              my_instructions.all_possible_diplotypes_table ||
              my_instructions.four_gamete_rule              ||
              my_instructions.ld_blocks                       );
}

void
analysis::write_user_options()
{
  vector<p_values>  dummy;
  writer<base_em_phenotype_map>  wrtr(my_output_state, my_summary_file, my_detail_file, my_instructions, dummy, my_current_region_lds);
  wrtr.init_output();
}

void
analysis::do_non_pool_analysis()
{
  bool  multiple_regions = my_instructions.regions.size() > 1;
  bool  first_region = true;
  list<instructions::region_data>::const_iterator  r_iter     = my_instructions.regions.begin();
  list<instructions::region_data>::const_iterator  r_end_iter = my_instructions.regions.end();
  for(; r_iter != r_end_iter; ++r_iter, first_region = false)
  {
    my_current_region_name = r_iter->name;
    my_output_state.first_block = "true";
    if(multiple_regions)
    {
      my_messages << "\n  Region -- " << my_current_region_name << endl;
      my_summary_file << "Region -- " << my_current_region_name << "\n" << endl;
    }

    // - Filter markers if specified by the user.
    //
    build_current_region(*r_iter);

    if(no_block_options_specified())
    {
      analyze_block(my_current_region);
    }
    else
    {
      if(my_instructions.sliding_window)
      {
        analyze_by_sliding_window();
      }

      if(my_instructions.four_gamete_rule)
      {
        analyze_by_four_gamete_rule_blocks();
      }

      if(my_instructions.ld_blocks)
      {
        analyze_by_ld_blocks();
      }
    }
  }
}

bool
analysis::no_block_options_specified() const
{
  return  ! (my_instructions.sliding_window   ||
             my_instructions.four_gamete_rule ||
             my_instructions.ld_blocks          );
}

void
analysis::analyze_by_sliding_window()
{
  my_messages << "  Preparing sliding window ......................................." << flush;
  block_vector  sw_blocks = sliding_window_blocks(*this)();

  if(! sw_blocks.empty())
  {
    write_done();
    size_t  block_count = sw_blocks.size();
    assert(block_count > 0);

    for(size_t b = 0; b < block_count; ++b)
    {
      size_t  block_number = b + 1;
      my_messages << "\n  Window " << block_number << " of " << block_count << "." << endl;

      const locus_group  block_group = build_block(sw_blocks, b);
      analyze_block(block_group, b ? "" : my_current_region_name, block_number);
    }
  }
  else
  {
    my_messages << endl;
    my_errors << priority(error) << "Insufficient number of loci in region, "
              << my_current_region_name << ", for sliding window option.  "
              << "Skipping sliding window ..." << endl;

    write_done();
  }
}

void
analysis::analyze_by_four_gamete_rule_blocks()
{
  my_messages << "\n  Determining blocks with four gamete rule ......................." << flush;
  block_vector  fg_blocks = four_gamete_rule_blocks(*this)();
  write_done();

  if(! fg_blocks.empty())
  {
    size_t  block_count = fg_blocks.size();
    for(size_t b = 0; b < block_count; ++b)
    {
      if(b)
      {
        my_messages << "\n";
      }

      size_t  block_number = b + 1;
      my_messages << "  Four-gamete-rule block " << block_number << " of " << block_count << "." << endl;

      const locus_group  block_group = build_block(fg_blocks, b);
      analyze_block(block_group, b ? "" : my_current_region_name, block_number);
    }
  }
  else
  {
    my_messages <<  "\n  No blocks found." << endl;
  }
}

void
analysis::analyze_by_ld_blocks()
{
  my_messages << "\n  Determining blocks via linkage disequilibrium .................." << flush;

  linkage_disequilibrium_blocks  ld_blks(*this);
  block_vector  ld_blocks = ld_blks();
  my_current_region_lds = ld_blks.ld_data();
  write_done();

  if(! ld_blocks.empty())
  {
    size_t  block_count = ld_blocks.size();
    for(size_t b = 0; b < block_count; ++b)
    {
      if(b)
      {
        my_messages << "\n";
      }

      size_t  block_number = b + 1;
      my_messages << "  LD block " << block_number << " of " << block_count << "." << endl;

      const locus_group  block_group = build_block(ld_blocks, b);
      analyze_block(block_group, b ? "" : my_current_region_name, block_number);
    }
  }
  else
  {
    my_messages <<  "\n  No blocks found." << endl;
  }
}

void
analysis::build_current_region(const instructions::region_data& region)
{
  my_current_region.clear();

  if(my_instructions.maf_filter)
  {
    my_summary_file << "\nMarkers excluded from analysis on basis of minor allele frequency:\n  ";
  }

  bool  markers_excluded = false;
  const locus_group&  loci = region.loci;
  size_t  locus_count = loci.size();
  for(size_t L = 0; L < locus_count; ++L)
  {
    if(! my_instructions.maf_filter)
    {
      my_current_region.push_back(loci[L]);
    }
    else
    {
      if(low_maf_freq(loci[L].second))
      {
        markers_excluded = true;
        const string&  marker_name = (loci[L].second)->name();
        my_errors << priority(information) << "Marker, " << marker_name << ", omitted from the current "
                  << "region on the basis of low minor allele frequency ..." << endl;
        my_summary_file << marker_name << " ";
      }
      else
      {
        my_current_region.push_back(loci[L]);
      }
    }
  }

  if(my_instructions.maf_filter)
  {
    if(! markers_excluded)
    {
      my_summary_file << "None";
    }

    my_summary_file << "\n" << endl;
  }
}

bool
analysis::low_maf_freq(const inheritance_model* locus)
{
  bool  low_maf = false;

  size_t  allele_count = locus->allele_count();
  if(allele_count > 2)
  {
    my_errors << priority(warning) << "Marker '" << locus->name() << "' not considered for "
              << "filtration by minor allele frequency as it is not a SNP." << endl;
  }
  else
  {
    double  freq1 = locus->get_allele(0).frequency();
    double  freq2 = locus->get_allele(1).frequency();

    assert(! isnan(freq1));
    assert(! isnan(freq2));

    if(allele_count < 2 || std::min(freq1, freq2) < my_instructions.maf_threshold)
    {
      low_maf =  true;
    }
  }

  return  low_maf;
}

void
analysis::write_done()
{
  if(! my_output_state.msg_interrupted)
  {
    my_messages << "... done." << endl;
  }
  else
  {
    my_messages << "                                                                  ... done" << endl;
    my_output_state.msg_interrupted = false;
  }
}

const locus_group
analysis::build_block(const block_vector& bv, size_t index) const
{
  locus_group  under_construction;

  assert(index < bv.size());

  size_t  begin = bv[index].first;
  size_t  end   = bv[index].second;
  for(size_t l = begin; l <= end; ++l)
  {
    under_construction.push_back(my_current_region[l]);
  }

  return  under_construction;
}

// - Perform analysis for the given loci.
//
void
analysis::analyze_block(const locus_group& loci,
                        const string& region, size_t block_number)
{
  if( (! (maximization_needed() || my_instructions.all_possible_diplotypes_table)) &&
      (my_instructions.four_gamete_rule || my_instructions.ld_blocks))
  {
    // - Nothing to do, but write blocks.
    //
    vector<p_values>  dummy;
    writer<base_em_phenotype_map>  wrtr(my_output_state, my_summary_file, my_detail_file, my_instructions, dummy, my_current_region_lds);
    wrtr.write(loci, region, block_number);
    my_output_state.first_block = false;

    return;
  }

  switch(my_instructions.analysis_unit)
  {
    case instructions::EACH_INDIVIDUAL:
    {
      partitioner::sub_pop_directory  members;
      my_partitioner.partition_population(members, loci, my_instructions.partitions);
      do_tasks<unrelated_em_phenotype_map, vector<member> >(members, loci, region, block_number);
      break;
    }

    case instructions::FAMILY_REP:
    {
      partitioner::sub_pop_directory  members;
      my_partitioner.partition_population(members, loci, my_instructions.partitions);
      do_tasks<family_em_phenotype_map, vector<member> >(members, loci, region, block_number);
      break;
    }

    case instructions::FAMILY_FOUNDERS:
    {
      partitioner::founders_sub_pop_directory  members;
      my_partitioner.partition_population(members, loci, my_instructions.partitions);
      do_tasks<founders_em_phenotype_map, pair<vector<member>, vector<member> > >(members, loci, region, block_number);
      break;
    }

    case instructions::POOL:
    {
      partitioner::sub_pop_directory  members;
      my_partitioner.partition_population(members, my_instructions.loci, my_instructions.partitions);
      do_tasks<pool_em_phenotype_map, vector<member> >(members, my_instructions.loci, region);
      break;
    }

    default:
      assert(false);
  }

  my_output_state.first_block = false;
}

template<typename M, typename V> void
analysis::do_tasks(const std::map<string, std::map<string, V> >& members,
                   const locus_group& loci, const string& region, size_t block_number)
{
  vector<p_values>  test_results;

  vector<const M*>  phenotype_maps;

  typename std::map<string, std::map<string, V> >::const_iterator  outer_iter     = members.begin();
  typename std::map<string, std::map<string, V> >::const_iterator  outer_end_iter = members.end();
  for(; outer_iter != outer_end_iter; ++outer_iter)
  {
    vector<const member_em_phenotype_map*>  m_phenotype_maps;

    // - Create em objects and maximize if necessary.
    //
    typename std::map<string, V>::const_iterator  inner_iter     = outer_iter->second.begin();
    typename std::map<string, V>::const_iterator  inner_end_iter = outer_iter->second.end();
    for(; inner_iter != inner_end_iter; ++inner_iter)
    {
      M*  map_ptr = new M(my_output_state, my_streams, my_filtered_mped,
                          inner_iter->second,
                          my_instructions, loci,
                          inner_iter->first, outer_iter->first);

      if(! map_ptr->empty())
      {
        if(maximization_needed())
        {
          map_ptr->maximize(my_instructions.epsilon, my_instructions.starting_points, my_dump_file,
                            my_instructions.dump, my_instructions.dump_cutoff, false);
        }

        //map_ptr->dump(my_summary_file);

        phenotype_maps.push_back(map_ptr);
        m_phenotype_maps.push_back(map_ptr);
      }
      else
      {
        string  sub_pop_name = map_ptr->sub_pop_name();
        if(sub_pop_name.empty())
        {
          my_errors << priority(error) << "Total population contains no usable data.  Skipping "
                    << "analysis ..." << endl;
          delete map_ptr;
          return;
        }
        else
        {
          my_errors << priority(error) << "Subpopulation, " << sub_pop_name << " contains "
                    << "no usable data.  Ignoring this subpopulation ..." << endl;
          delete map_ptr;
        }
      }
    }

    // - Do likelihood ratio and permutation tests.
    //
    if(my_instructions.likelihood_ratio_test || my_instructions.compute_empirical_pvalue)
    {
      my_messages << "    Doing likelihood ratio test(s) ..............................." << flush;

      if(m_phenotype_maps.size() > 1)
      {
        try
        {
          test_results.push_back(do_tests(m_phenotype_maps, loci));
        }
        catch(const bad_alloc&)
        {
          if(! my_output_state.msg_interrupted)
          {
            my_messages << endl;
            my_output_state.msg_interrupted = true;
          }

          my_errors << priority(error) << "Not enough memory available to do "
                    << "likelihood ratio test(s) and/or to calculate empirical p-value(s). "
                    << "Skipping these tasks ..." << endl;
        }
      }
      else
      {
        if(! my_output_state.msg_interrupted)
        {
          my_messages << endl;
          my_output_state.msg_interrupted = true;
        }

        if(outer_iter->first == "")
        {
           my_errors << priority(error) << "Multiple subpopulations not found.  "
                     << "Skipping likelihood ratio test/empirical p-value computation ..." << endl;
        }
        else
        {
           my_errors << priority(error) << "subpopulation, " << outer_iter->first
                     << ", Does not contain multiple subpopulations.  "
                     << "Skipping likelihood ratio test/empirical p-value computation ..." << endl;
        }
      }

      write_done();
    }     // END INNER PARTITION LOOP
  }     // END OUTER PARTITION LOOP

  // - Write results.
  //
  if(! phenotype_maps.empty())
  {
    writer<M>  wrtr(my_output_state, my_summary_file,
                    my_detail_file, my_instructions, test_results, phenotype_maps, my_current_region_lds);
    wrtr.write(region, block_number);

    for(size_t m = 0; m < phenotype_maps.size(); ++m)
    {
      delete phenotype_maps[m];
    }
  }
  else
  {
    my_errors << priority(error) << "No valid subpopulations found.  Skipping analysis/window ..." << endl;

    // - Write user output to summary file.
    //
    writer<M>  wrtr(my_output_state, my_summary_file,
                    my_detail_file, my_instructions, test_results, phenotype_maps, my_current_region_lds);
    wrtr.write(region, block_number);

    return;
  }
}

// - Likelihood ratio test and empirical p-value.
//
p_values
analysis::do_tests(const vector<const member_em_phenotype_map*>& phenotype_maps, const locus_group& loci)
{
  p_values  result;

  result.outer_sub_pop_name = phenotype_maps[0]->outer_sub_pop_name();

  sub_pop_shuffler  shuffler;
  shuffler.set_sub_pops(phenotype_maps);

  // - Calculate denominator of likelihood ratio.
  //
  rebuilt_em_phenotype_map  whole_map(my_output_state, my_streams, loci, shuffler.get_whole_pop(), phenotype_maps);
  whole_map.maximize(my_instructions.epsilon, my_instructions.starting_points, my_dump_file);
  size_t  whole_ipc = whole_map.independent_param_count();
  double  whole_ln_likelihood = whole_map.max_ln_likelihood();

  result.whole_ln_like = whole_ln_likelihood;

  // - Calculate numerator of likelihood ratio.
  //
  double  composite_ln_likelihood = 0;
  size_t  map_count = phenotype_maps.size();

  assert(map_count > 1);

  for(size_t m = 0; m < map_count; ++m)
  {
    assert(phenotype_maps[m]->is_maximized());
    assert(phenotype_maps[m]->outer_sub_pop_name() == result.outer_sub_pop_name);
    double  sub_pop_ln_likelihood = phenotype_maps[m]->max_ln_likelihood();
    composite_ln_likelihood += sub_pop_ln_likelihood;
  }

  result.composite_ln_like = composite_ln_likelihood;

  // - Theoreticaly the test statistic will always be non-negative.  If it isn't
  //   (perhaps for numerical reasons), chdtrc() will return QNAN and the program
  //   will show "---" as the result.
  //
  size_t  degrees_of_freedom = (map_count - 1) * whole_ipc;
  double  test_statistic = 2 * (composite_ln_likelihood - whole_ln_likelihood);

  result.degrees_of_freedom = degrees_of_freedom;
  result.test_statistic = test_statistic;
  result.asymptotic = chdtrc(degrees_of_freedom, test_statistic);

  if(my_instructions.compute_empirical_pvalue)
  {
    result.empirical = do_permutation_test(test_statistic, whole_ln_likelihood,
                                           shuffler, phenotype_maps, loci, result.permutations);
  }

  return  result;
}

double
analysis::do_permutation_test(double test_statistic,
                              double whole_ln_likelihood, sub_pop_shuffler& suf,
                              const vector<const member_em_phenotype_map *>& phenotype_maps,
                              const locus_group& loci, size_t& total_permutations)
{
  size_t fixed_permutations = my_instructions.permutations;
  double max_permutations   = my_instructions.max_permutations;
  double min_permutations   = my_instructions.min_permutations;
  double alpha              = my_instructions.confidence;
  double width              = my_instructions.width;

  if(   !finite(width) || width < 100*std::numeric_limits<double>::epsilon()
     || !finite(alpha) || alpha < 100*std::numeric_limits<double>::epsilon() )
  {
    width = 0.2;
    alpha = 0.95;
  }

  // See binomial trails theory for explaination
  double precision = pow(inv_normal_cdf((1.+alpha)/2.)/width, 2);

  if( fixed_permutations )
  {
    min_permutations = fixed_permutations;
    max_permutations = fixed_permutations;
  }

#if 0
  cout << "precision:" << precision << endl;
  cout << "Original:" << endl;
  suf.dump(cout);
#endif

  double emp_p = 1.0;

  size_t significant_permutations = 1;
  total_permutations       = 1;

  for( size_t i = 0; i < max_permutations; ++i )
  {
    bool valid_shuffle = suf.do_shuffle(my_instructions.seed);

    if( !valid_shuffle )
      continue;

#if 0
  cout << "Permutation " << i+1 << ":" << endl;
  suf.dump(cout);
#endif

    // - Calculate numerator of likelihood ratio.
    //
    double  composite_ln_likelihood = 0.;

    const vector<vector<member> >& new_sub_pops = suf.get_new_sub_pops();

    for( size_t s = 0; s < new_sub_pops.size(); ++s )
    {
      const vector<member>& a_new_sub_pop = new_sub_pops[s];

      rebuilt_em_phenotype_map a_new_sub_map(my_output_state, my_streams, loci, a_new_sub_pop, phenotype_maps);

      a_new_sub_map.maximize(my_instructions.epsilon, my_instructions.starting_points, my_dump_file);

      double new_sub_ln_likelihood = a_new_sub_map.max_ln_likelihood();

      composite_ln_likelihood += new_sub_ln_likelihood;
    }

    double new_statistic = 2 * (composite_ln_likelihood - whole_ln_likelihood);

    if( new_statistic > test_statistic )
      ++significant_permutations;

    ++total_permutations;

    emp_p = double(significant_permutations) / double(total_permutations);

    double m = i * emp_p / (1.0 - emp_p);

#if 0
  cout << "test_statistic          = " << test_statistic << endl;
  cout << "composite_ln_likelihood = " << composite_ln_likelihood << endl;
  cout << "whole_ln_likelihood     = " << whole_ln_likelihood << endl;
  cout << "new_statistic           = " << new_statistic << endl;
  cout << "emp_p                   = " << emp_p << endl;
  cout << "m                       = " << m << endl << endl;
#endif

    if( i > min_permutations && m > precision )
      break;
  }

  --total_permutations;    // Since this starts with a value of 1.

#if 0  // Degrees of freedom experiment
  cout << "\nppp empirical pvalue " << emp_p << endl;
#endif

  return emp_p;
}


//============================================================================
// IMPLEMENTATION:  blocks
//============================================================================
//
blocks::blocks(analysis& a)
      : my_analysis(a)
{}

blocks::~blocks()
{}

const block_vector&
blocks::operator()() const
{
  return  my_blocks;
}


//============================================================================
// IMPLEMENTATION:  sliding_window_blocks
//============================================================================
//
sliding_window_blocks::sliding_window_blocks(analysis& a)
      : blocks(a)
{
  size_t  locus_count = my_analysis.my_current_region.size();

  size_t  begin = 0;
  size_t  end = my_analysis.my_instructions.window_width - 1;

  if(my_analysis.my_instructions.window_width <= locus_count)
  {
    while(end < locus_count)
    {
      my_blocks.push_back(make_pair(begin, end));
      ++begin;
      ++end;
    }
  }
}


//============================================================================
// IMPLEMENTATION:  four_gamete_rule_blocks
//============================================================================
//
four_gamete_rule_blocks::four_gamete_rule_blocks(analysis& a)
      : blocks(a)
{
  size_t  locus_count = my_analysis.my_current_region.size();
  if(locus_count > 1)
  {
    size_t begin = 0;
    size_t end   = 1;

    for(; end < locus_count; ++end)
    {
      if(recomb(begin, end))
      {
        size_t  last = end - 1;
        if(last > begin)         // No blocks of 1 locus.
        {
          my_blocks.push_back(make_pair(begin, last));
        }

        begin = end;
      }
    }

    // - This block terminated by end of region, not a recombination.
    //
    size_t  last = end - 1;
    if(last > begin)
    {
      my_blocks.push_back(make_pair(begin, last));
    }
  }
}

// - Estimate haplotype frequencies for pairs consisting of the end locus and each
//   locus in the range [begin, end).  If any of these pairs results in four
//   haplotype frequencies greater than the threshold, a recombination is deemed
//   to have taken place.
//
bool
four_gamete_rule_blocks::recomb(size_t begin, size_t end)
{
  my_analysis.my_output_state.allow_progress_msg = false;

  bool  recomb_found = false;
  partitioner::partition_vector  whole_pop(2);

  for(size_t l = begin; l < end; ++l)
  {
    base_em_phenotype_map*  two_locus_map_ptr = NULL;
    locus_group  duo;

    duo.push_back(my_analysis.my_current_region[l]);
    duo.push_back(my_analysis.my_current_region[end]);
    switch(my_analysis.my_instructions.analysis_unit)
    {
      case instructions::EACH_INDIVIDUAL:
      {
        partitioner::sub_pop_directory  members;
        my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

        assert(members.find("") != members.end());
        assert(members[""].find("") != members[""].end());

        two_locus_map_ptr = new unrelated_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                           members[""][""], my_analysis.my_instructions, duo,
                                                           "", "");
        break;
      }

      case instructions::FAMILY_REP:
      {
        partitioner::sub_pop_directory  members;
        my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

        assert(members.find("") != members.end());
        assert(members[""].find("") != members[""].end());

        two_locus_map_ptr = new family_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                        members[""][""], my_analysis.my_instructions, duo,
                                                        "", "", true);
        break;
      }

      case instructions::FAMILY_FOUNDERS:
      {
        partitioner::founders_sub_pop_directory  members;
        my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

        assert(members.find("") != members.end());
        assert(members[""].find("") != members[""].end());

        two_locus_map_ptr = new founders_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                          members[""][""], my_analysis.my_instructions, duo,
                                                          "", "");
        break;
      }

      default:
        assert(false);
    }

    // - recomb found is true only if # of haps w. (freqs > threshold) == 4.
    //
    two_locus_map_ptr->maximize(my_analysis.my_instructions.epsilon, my_analysis.my_instructions.starting_points, my_analysis.my_dump_file);

    size_t  hap_count = common_freqs(two_locus_map_ptr, my_analysis.my_instructions.fg_threshold);
    assert(hap_count <= 4);

    if(hap_count == 4)
    {
      recomb_found = true;
      break;
    }

    delete  two_locus_map_ptr;
  }

  my_analysis.my_output_state.allow_progress_msg = true;
  return  recomb_found;
}

size_t
four_gamete_rule_blocks::common_freqs(const base_em_phenotype_map* pheno_map, double threshold) const
{
  size_t  freq_count = 0;
  const set<hap_freq, greater<hap_freq> >&  freqs = pheno_map->final_frequencies();

  set<hap_freq, greater<hap_freq> >::const_iterator  f_iter     = freqs.begin();
  set<hap_freq, greater<hap_freq> >::const_iterator  f_end_iter = freqs.end();
  for(; f_iter != f_end_iter; ++f_iter)
  {
    if(f_iter->freq > threshold)
    {
      ++freq_count;
    }
  }

  return  freq_count;
}


//============================================================================
// IMPLEMENTATION:  linkage_disequilibrium_blocks
//============================================================================
//
linkage_disequilibrium_blocks::linkage_disequilibrium_blocks(analysis& a)
      : blocks(a)
{
  size_t  locus_count = my_analysis.my_current_region.size();
  if(locus_count > 1)
  {
    size_t begin = 0;         // first_valid_snp();
    size_t end   = 1;         // begin + 1;

    for(; end < locus_count; ++end)
    {
      if(! ld(end - 1, end))     // Break between blocks
      {
        size_t  last = end - 1;
        if(last > begin)         // Single locus blocks not allowed.
        {
          my_blocks.push_back(make_pair(begin, last));
        }

        begin = end;
      }
    }

    // - This block terminated by end of region, lack of ld.
    //
    size_t  last = end - 1;
    if(last > begin)
    {
      my_blocks.push_back(make_pair(begin, last));
    }
  }
}

const vector<ld_record>&
linkage_disequilibrium_blocks::ld_data() const
{
  return  my_lds;
}

// - Does LD between the given loci exceed the threshold?
//
bool
linkage_disequilibrium_blocks::ld(size_t locus1, size_t locus2)
{
  my_analysis.my_output_state.allow_progress_msg = false;

  // - Estimate gamete frequencies.
  //
  partitioner::partition_vector  whole_pop(2);

  base_em_phenotype_map*  two_locus_map_ptr = NULL;
  locus_group  duo;

  duo.push_back(my_analysis.my_current_region[locus1]);
  duo.push_back(my_analysis.my_current_region[locus2]);
  switch(my_analysis.my_instructions.analysis_unit)
  {
    case instructions::EACH_INDIVIDUAL:
    {
      partitioner::sub_pop_directory  members;
      my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

      assert(members.find("") != members.end());
      assert(members[""].find("") != members[""].end());

      two_locus_map_ptr = new unrelated_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                         members[""][""], my_analysis.my_instructions, duo,
                                                         "", "");
      break;
    }

    case instructions::FAMILY_REP:
    {
      partitioner::sub_pop_directory  members;
      my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

      assert(members.find("") != members.end());
      assert(members[""].find("") != members[""].end());

      two_locus_map_ptr = new family_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                      members[""][""], my_analysis.my_instructions, duo,
                                                      "", "", true);
      break;
    }

    case instructions::FAMILY_FOUNDERS:
    {
      partitioner::founders_sub_pop_directory  members;
      my_analysis.my_partitioner.partition_population(members, duo, whole_pop);

      assert(members.find("") != members.end());
      assert(members[""].find("") != members[""].end());

      two_locus_map_ptr = new founders_em_phenotype_map(my_analysis.my_output_state, my_analysis.my_streams, my_analysis.my_filtered_mped,
                                                        members[""][""], my_analysis.my_instructions, duo,
                                                        "", "");
      break;
    }

    default:
      assert(false);
  }

  two_locus_map_ptr->maximize(my_analysis.my_instructions.epsilon, my_analysis.my_instructions.starting_points, my_analysis.my_dump_file);

  double  lld = lewontins_ld(two_locus_map_ptr);
  my_lds.push_back(ld_record(duo[0].second->name(), duo[1].second->name(), lld));

  delete  two_locus_map_ptr;
  my_analysis.my_output_state.allow_progress_msg = true;

  return  lld >= my_analysis.my_instructions.ld_threshold;
}

// - Weir, p120.
//
double
linkage_disequilibrium_blocks::lewontins_ld(const base_em_phenotype_map* em_map)
{
  check_allele_frequencies(em_map);

  const locus_group&  loci = em_map->loci();
  assert(loci.size() == 2);
  assert(loci[0].second->allele_count() == 2);
  assert(loci[1].second->allele_count() == 2);

  const em_haplotype_map&  final_haplotype_map = em_map->final_haplotypes();

  assert(! final_haplotype_map.haplotypes().empty());

  hap_seq  seq_00;
  seq_00.push_back(0);
  seq_00.push_back(0);

  // - Arbitrarily make this the haplotype composed of the 0th allele at each locus.
  //
  double  freq_ab;
  size_t  index_00 = final_haplotype_map.hap_seq_to_index(seq_00);
  if(index_00 != (size_t)(-1))
  {
    const em_haplotype&  hap_ab = final_haplotype_map.get_haplotype(index_00);
    freq_ab = hap_ab.new_freq(em_map->total_hap_count());
  }
  else
  {
    freq_ab = 0.0;
  }

  double  freq_a       = loci[0].second->get_allele(0).frequency();
  double  freq_b       = loci[1].second->get_allele(0).frequency();
  double  freq_a_compl = 1 - freq_a;
  double  freq_b_compl = 1 - freq_b;

  double  d_ab = freq_ab - freq_a * freq_b;
  double  d_ab_prime;


  if(d_ab < 0)
  {
    d_ab_prime = d_ab / max(- freq_a * freq_b, - freq_a_compl * freq_b_compl);
  }
  else if(d_ab > 0)
  {
    d_ab_prime = d_ab / min(freq_a_compl * freq_b, freq_a * freq_b_compl);
  }
  else
  {
    d_ab_prime = 0;
  }

  if(d_ab_prime - 1.0  > LD_ERROR_TOLERANCE)
  {
    if(! my_analysis.my_output_state.msg_interrupted)
    {
      my_analysis.my_messages << endl;
      my_analysis.my_output_state.msg_interrupted = true;
    }

    my_analysis.my_errors << priority(error) << "Calculated D' between markers " << loci[0].second->name() << " "
                          << loci[1].second->name() << " is greater than 1.  Allele frequecies in the "
                          << "marker locus description file may differ substantially from those in the pedigree data. "
                          << "See detail file for D' values." << endl;
  }

  return  d_ab_prime < 0 ? - d_ab_prime : d_ab_prime;
}

// - Compare allele frequencies implied by estimated gamete frequencies with
//   those from the marker locus description file and issue a warning if they
//   differ substantially.
//
void
linkage_disequilibrium_blocks::check_allele_frequencies(const base_em_phenotype_map* em_map)
{
  // - Get gamete frequencies.
  //
  const locus_group&  loci = em_map->loci();
  assert(loci.size() == 2);

  const em_haplotype_map&  gametes = em_map->final_haplotypes();
  size_t  total_hap_count = em_map->total_hap_count();

  hap_seq seq_00;
  seq_00.push_back(0);
  seq_00.push_back(0);
  hap_seq seq_01;
  seq_01.push_back(0);
  seq_01.push_back(1);
  hap_seq seq_11;
  seq_11.push_back(1);
  seq_11.push_back(1);

  size_t  index_00 = gametes.hap_seq_to_index(seq_00);
  size_t  index_01 = gametes.hap_seq_to_index(seq_01);
  size_t  index_11 = gametes.hap_seq_to_index(seq_11);

  double  freq_00 = index_00 != (size_t)(-1) ? gametes.get_haplotype(index_00).new_freq(total_hap_count) : 0.0;
  double  freq_01 = index_01 != (size_t)(-1) ? gametes.get_haplotype(index_01).new_freq(total_hap_count) : 0.0;
  double  freq_11 = index_11 != (size_t)(-1) ? gametes.get_haplotype(index_11).new_freq(total_hap_count) : 0.0;

  // - Calculate implied allele frequencies.  Arbitrarily choose allele with index
  //   of 0.
  //
  double  locus0_allele0 = freq_00 + freq_01;
  double  locus1_allele0 = 1.0 - (freq_01 + freq_11);

  // - Test frequencies against those provided in the marker locus description file.
  //
  test_allele_frequency(loci[0].second, locus0_allele0);
  test_allele_frequency(loci[1].second, locus1_allele0);
}

void
linkage_disequilibrium_blocks::test_allele_frequency(const MLOCUS::inheritance_model* locus, double allele_frequency)
{
  double  freq_diff = locus->get_allele(0).frequency() - allele_frequency;
  double  abs_freq_diff = freq_diff > 0 ? freq_diff : - freq_diff;

  string  locus_name = locus->name();

  if(abs_freq_diff > ALLELE_FREQUENCY_TOLERANCE * allele_frequency &&
     locus_name != my_last_freq_warning                               )
  {
    if(! my_analysis.my_output_state.msg_interrupted)
    {
      my_analysis.my_messages << endl;
      my_analysis.my_output_state.msg_interrupted = true;
    }

    my_analysis.my_errors << priority(warning) << "Estimated allele frequencies for marker, " << locus->name()
                          << ", deviates by more than " << ALLELE_FREQUENCY_TOLERANCE * 100 << " percent "
                          << "from frequencies supplied in the marker locus description file." << endl;
    my_last_freq_warning = locus_name;
  }
}

}
}
