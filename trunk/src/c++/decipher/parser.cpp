//============================================================================
// File:      parser.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   4/5/4 - created.                         djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              parser 
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/parser.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
parser::parser(const RefMultiPedigree* mp, const decipher_data& data, 
               ostream& m, cerrorstream& errors, genome_description* gd)
    : BasicParser(errors), my_ref_mped(mp), my_data(data), my_messages(m), my_genome(gd)
{
  my_messages << endl << "Parsing DECIPHER analyses ...\n" << endl;

  reset();
  my_analysis_id = 0;
}

void
parser::parse_symbols(const SymbolTable* syms)
{}

void 
parser::parse_parameter(const LSFBase* param)
{}

void  
parser::parse_test_parameter(const LSFBase* param)
{}

void
parser::parse(const LSFBase* params)
{
  parse_test_parameter_section(params);
}

// - Top level.  Iterate through parameters of decipher_analysis
//   block.  Client code must insure that params->name() is 'decipher'
//   or 'decipher_analysis'.
//
void  
parser::parse_test_parameter_section(const LSFBase* params)
{
  assert(params && my_ref_mped);
  
  reset();  
  init_parse();
  print_header();
  parse_out_attr(params);
  
  vector<LSFBase*>  region_ptrs;  

  if(params->List() != 0)
  {
    LSFList::const_iterator  iter;
    for(iter = params->List()->begin(); 
        iter != params->List()->end() && my_instructions.valid; 
        ++iter)
    {
      string  parameter = toUpper((*iter)->name());
    
      if(parameter == "TITLE") 
      {
        parse_title(*iter); 
      }
      else if(parameter == "REGION")
      {
        region_ptrs.push_back(*iter);
      }
      else if(parameter == "EPSILON")
      {
        parse_epsilon(*iter);
      }
      else if(parameter == "STARTING_POINTS")
      {
        parse_starting_points(*iter);
      }
      else if(parameter == "DUMP")
      {
        parse_dump(*iter);
      }    
      else if(parameter == "SEED")
      {
        parse_seed(*iter);
      }
      else if(parameter == "FILTERS")
      {
        parse_filters(*iter);
      }      
      else if(parameter == "BLOCKS")
      {
        parse_blocks(*iter);
      }        
      else if(parameter == "DATA")
      {
        parse_data(*iter);
      }
      else if(parameter == "TASKS")
      {
        parse_tasks(*iter);
      }
      else
      {
        errors << priority(error) << "parameter '" << parameter << "' not recognized.  "
               << "Ignoring parameter ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "No parameters read for Decipher analysis block." << endl;
  }
  
  size_t  region_count = region_ptrs.size();
  if(my_instructions.analysis_unit != instructions::POOL)
  {
    if(region_count)
    {
      for(size_t r = 0; r < region_count; ++r)
      {
        parse_region(region_ptrs[r]);
      }
    }
    else
    {
      errors << priority(critical) << "No parameter, region, found with parameter, "
             << "analysis_unit, set to '" << instructions::unit_2_string(my_instructions.analysis_unit)
             << "'.  Ignoring decipher block ..." << endl;
      my_instructions.valid = false;
    }
  }
  else
  {
    if(region_count)
    {
      errors << priority(warning) << "parameter, region, given with parameter, analysis_unit, "
             << "set to 'pool'.  Ignoring parameter, region ... " << endl;
    }
  }
  
  // - Constraints.
  //
  check_for_required_sub_pops();
  check_dump_conditions();
  
  if(my_instructions.analysis_unit != instructions::POOL)
  {
    if(my_instructions.valid)                // critical
    {
      check_for_regions();
    }
    
    if(my_instructions.valid)                // critical
    {
      validate_regions();
    }
    
    if(my_instructions.valid)
    {
      check_four_gamete_threshold();
      check_for_mld_file();
    }
  }
  else
  {
    if(my_instructions.valid)
    {
      check_pool_locus_count();              // critical
    }
    
    if(my_instructions.valid)
    {
      check_sliding_window();
      check_maf_filter();
      check_four_gamete_rule();
      check_ld();
      build_loci_for_pools();
    }
  }

  print_footer();
}

void
parser::reset()
{
  my_instructions.reset();
  my_pool_loci.clear();
  my_trait_registry.clear();
  my_string_registry.clear();
}

void parser::init_parse()
{
  // Create the analysis identifier string.
  //
  ++my_analysis_id;
  string  id = long2str(my_analysis_id);

  // Set the title and output to defaults.
  //
  my_instructions.title          = "Analysis " + id;
  my_instructions.file_name_root = "decipher_analysis" + id; 
}

void
parser::check_for_required_sub_pops()
{
  assert(my_instructions.partitions.size() == PARTITION_COUNT);
  if(my_instructions.partitions[1].sub_pops.size() < 2)
  {
    if(my_instructions.likelihood_ratio_test)
    {
      errors << priority(error) << "likelihood ratio test requires multiple subpopulations.  "
                                << "Skipping likelihood ratio test ..." << endl;
      my_instructions.likelihood_ratio_test = false;
    }
    
    if(my_instructions.compute_empirical_pvalue)
    {
      errors << priority(error) << "computation of an empirical p-value requires multiple subpopulations.  "
                                << "Skipping computation of empirical p-value ..." << endl;
      my_instructions.compute_empirical_pvalue = false;
    }    
  }
}

void
parser::check_dump_conditions()
{
  bool  applicable = my_instructions.pop_freq || my_instructions.most_likely_diplotypes;
  if(my_instructions.dump && ! applicable)
  {
    errors << priority(error) << "Parameter, dump, applicable only if population frequencies "
                              << "or most likely diplotypes are selected.  Ignoring 'dump' ..." << endl;
    my_instructions.dump = false;
  }
}

// - Invalidate analysis if there is not at least one region.
//
void  
parser::check_for_regions()
{
  if(my_instructions.regions.empty())
  {
    errors << priority(critical) << "No valid region given.  Skipping analysis ..." << endl;
    my_instructions.valid = false;    
  }
}

// - Build and check loci for each region.  Remove those that are invalid.
//
void
parser::validate_regions()
{
  list<instructions::region_data>::iterator  r_iter     = my_instructions.regions.begin();
  list<instructions::region_data>::iterator  r_end_iter = my_instructions.regions.end();
  for(; r_iter != r_end_iter; ++r_iter)
  {
    build_loci(*r_iter);
    
    if(r_iter->valid)
    {
      check_locus_count(*r_iter);
    }
    else
    {
      continue;
    }
    
    if(r_iter->valid)
    {
      check_for_codominance(*r_iter);
    }
    else
    {
      continue;
    }    
    
    if(r_iter->valid)
    {
      check_x_linkage(*r_iter);
    }
    else
    {
      continue;
    }        
    
    if(r_iter->valid)
    {
      check_related_x_linkage(*r_iter);
    }
    
    if(r_iter->valid)
    {
      check_snp_requirement(*r_iter);
    }
  }
  
  // - Remove invalid regions.
  //
  r_iter = my_instructions.regions.begin();
  my_instructions.regions.remove_if(instructions::invalid_region);
  
  check_for_regions();
}

// - Loci may be built using the genome description or directly from
//   first and last indices if they are given.
//
void
parser::build_loci(instructions::region_data& region)
{
  assert(region.valid);
  assert((region.first == MAX_SIZE && region.last == MAX_SIZE) ||
         (region.first != MAX_SIZE && region.last != MAX_SIZE)   );  

  if(region.first == MAX_SIZE)
  {
    genome_description::region_type  reg =
          const_cast<genome_description*>(my_genome)->region(region.name);

    genome_description::region_type::iterator  l_iter     = reg.begin();
    genome_description::region_type::iterator  l_end_iter = reg.end();
    for(; l_iter != l_end_iter; ++l_iter)
    {
      region.loci.push_back(make_pair((*l_iter).marker_index(), (*l_iter).locus()));
    }
  }
  else
  {
    region.valid = check_first_last_order(region);
    if(region.valid)    
    {
      const marker_order&  original_marker_order = my_data.original_marker_order();    
      APP::marker_order::iterator  first  = original_marker_order.find(region.first);
      APP::marker_order::iterator  last   = original_marker_order.find(region.last);
      APP::marker_order::iterator  end    = original_marker_order.end();    
      assert(first != end && last != end);
    
      ++last;               // marker range is [first, last]
      for(; first != last; ++first)
      {
        region.loci.push_back(make_pair(first.imodel_index(), *first));
      }
    }
  }
}

bool
parser::check_first_last_order(const instructions::region_data& region)
{
  const marker_order&  original_marker_order = my_data.original_marker_order();
  APP::marker_order::iterator  first  = original_marker_order.find(region.first);
  APP::marker_order::iterator  last   = original_marker_order.find(region.last);
  APP::marker_order::iterator  end    = original_marker_order.end();    
  assert(first != end && last != end);

  bool  order_valid = last > first;
  
  if(! order_valid)
  {
    errors << priority(error) << "'first' marker is not before 'last' marker in region '"
           << region.name << "  Skipping region ..." << endl;
  }
  
  return  order_valid;
}

// - If there are fewer than two loci, issue error message and invalidate region.
//
void
parser::check_locus_count(instructions::region_data& region)
{
  if(region.loci.size() < 2)
  {
    errors << priority(error) << "Fewer that two valid loci specified for haplotype "
           << "region '" << region.name << "'.  Skipping region ..." << endl;
    region.valid = false;
  }
}

// - If there are any non-codominant markers in the haplotyping region, issue an
//   error message and invalidate region.
//
void
parser::check_for_codominance(instructions::region_data& region)
{
  size_t  locus_count = region.loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    if(! (region.loci[l].second)->codominant())
    {
      errors << priority(error) << "One or more markers in haplotyping region '"
             << region.name << "' are not codominant.  Skipping region ..." << endl;
      region.valid = false;
      break;      
    }
  }
}

// - Determine if markers are x-linked.  If some are and some aren't, issue and error
//   and invalidate region.
//
void
parser::check_x_linkage(instructions::region_data& region)
{
  size_t  locus_count = region.loci.size();
  for(size_t l = 0; l < locus_count; ++l)
  {
    if(l == 0)
    {
      region.x_linked = (region.loci[l].second)->is_x_linked();
    }
    else 
    {
      if(region.x_linked != (region.loci[l].second)->is_x_linked())
      {
        errors << priority(error) << "Some markers in haplotyping region '"
               << region.name << "' are x-linked and some "
               << "are not.  Skipping region ..." << endl;
        region.valid = false;
        break;
      }
    }
  }
}

void
parser::check_related_x_linkage(instructions::region_data& region)
{
  if(region.x_linked && (! my_instructions.analysis_unit == instructions::EACH_INDIVIDUAL))
  {
    errors << priority(error) << "Use of x-linked markers only allowed with analysis_unit "
           << "set to 'each_individual'.  "
           << "Skipping x-linked haplotyping region '" << region.name << "' ..." << endl;
    region.valid = false;
  }
}

void  
parser::check_snp_requirement(instructions::region_data& region)
{
  if(my_instructions.four_gamete_rule || my_instructions.ld_blocks)
  {
    size_t  locus_count = region.loci.size();
    for(size_t l = 0; l < locus_count; ++l)
    {
      if((region.loci[l].second)->allele_count() > 2)
      {
        errors << priority(error) << "Sub_block, blocks, parameter, four_gamete_rule or ld, set to 'true', but "
               << "one or more markers in haplotyping region '"
               << region.name << "' are not snp's.  Skipping region ..." << endl;
        region.valid = false;
        break;      
      }
    }    
  }
}

void
parser::build_loci_for_pools()
{
  // - my_instructions.pool_loci contains data read from analysis_block of 
  //      parameter file.
  //   my_pool_loci contains penetrance models created from information in 
  //      my_instructions.pool_loci.
  //   my_instructions.loci contains indices of and pointers to penetrance models
  //      in my_pool_loci.
  //
  size_t  locus_count = my_instructions.pool_loci.size();
  
  // - Establish vector size here so that no reallocation will take place later.
  //   This is necessary so that pointers to the elements in the vector will 
  //   remain valid when new elements are added.
  //
  my_pool_loci.reserve(locus_count);
  
  for(size_t l = 0; l < locus_count; ++l)
  {
    MLOCUS::genotype_model  gm("");
    
    set<instructions::pool_locus::allele>::const_iterator  a_iter     = my_instructions.pool_loci[l].alleles.begin();
    set<instructions::pool_locus::allele>::const_iterator  a_end_iter = my_instructions.pool_loci[l].alleles.end();
    for(; a_iter != a_end_iter; ++a_iter)
    {
      gm.add_allele(a_iter->name);
    }
    
    MLOCUS::inheritance_model  pm(gm);
    pm.set_name(my_instructions.pool_loci[l].name);
    
    my_pool_loci.push_back(pm);
    my_instructions.loci.push_back(make_pair(my_pool_loci.size() - 1, &my_pool_loci.back())); 
  }
}


// - If there are fewer than two loci, issue error message and invalidate instructions.
//
void
parser::check_pool_locus_count()
{
  if(my_instructions.pool_loci.size() < 2)
  {
    errors << priority(critical) << "Fewer that two valid loci specified in the pools sub-block "
           << "Skipping analysis ..." << endl;
    my_instructions.valid = false;
  }
}

// - Sliding window not yet implemented for pools.
//
void
parser::check_sliding_window()
{
  if(my_instructions.sliding_window)
  {
    errors << priority(error) << "Sliding window not allowed if data sub-block parameter, "
           << "analysis_unit, equals 'pool'.  "
           << "Ignoring parameter, sliding_window ..." << endl;
    my_instructions.sliding_window = false;
  }
}

// - Four gamete rule does not apply to pools.
//
void
parser::check_four_gamete_rule()
{
  if(my_instructions.four_gamete_rule)
  {
    errors << priority(error) << "Four gamete rule not allowed if data sub-block parameter, "
           << "analysis_unit, equals 'pool'.  "
           << "Ignoring parameter, four_gamete_rule ..." << endl;
    my_instructions.four_gamete_rule = false;
  }
}

// - Four gamete rule does not apply to pools.
//
void
parser::check_maf_filter()
{
  if(my_instructions.maf_filter)
  {
    errors << priority(error) << "Minor allele frequency filtering not allowed if data sub-block parameter, "
           << "analysis_unit, equals 'pool'.  "
           << "Ignoring parameter, maf_filter ..." << endl;
    my_instructions.maf_filter = false;
  }
}

// - Four gamete rule does not apply to pools.
//
void
parser::check_ld()
{
  if(my_instructions.ld_blocks)
  {
    errors << priority(error) << "LD block determination not allowed if data sub-block parameter, "
           << "analysis_unit, equals 'pool'.  "
           << "Ignoring parameter, ld ..." << endl;
    my_instructions.ld_blocks = false;
  }
}

// - If the four gamete threshold is not at least equal to ten times the em algorithm
//   convergence criterium, make it so.
//
void
parser::check_four_gamete_threshold()
{
  if(my_instructions.four_gamete_rule)
  {
    double  fg_threshold_min = my_instructions.epsilon * 10;
    if(my_instructions.fg_threshold < fg_threshold_min)
    {
      my_instructions.fg_threshold = fg_threshold_min;
      
      errors << priority(error) << "Value of four_gamete_rule attribute, threshold, is less "
             << "than 10 times the em algorithm convergence criterion of " << my_instructions.epsilon
             << ".  setting threshold to " << fg_threshold_min << " ..." << endl;
    }
  }
}

void 
parser::print_header()
{
  my_messages << "\nParsing new analysis block ...\n" << endl;
}

void 
parser::print_footer()
{
  my_messages << "\nParsing of " << my_instructions.title << " complete.  ";

  if(my_instructions.valid)
  {
    my_messages << "Analysis valid." << endl;
  }
  else
  {
    my_messages << "Analysis is invalid.  Check your file." << endl;
  }

  my_messages << "------------------------------------------------------------------------------"
           << endl;
}

void
parser::parse_out_attr(const LSFBase* param)
{
  AttrList* a_list = param->attrs();
  
  if(a_list)
  {
    if(a_list->has_attr("OUT"))
    {
      string  out_value = a_list->StringAttr("OUT");
      if(! out_value.empty())
      {
        my_instructions.file_name_root = out_value;
      }
      else
      {
        errors << priority(warning) << "No value given for attribute, out. "
               << "Ignoring attribute ..." << endl;
      }
    }
  }
}

void
parser::parse_title(const LSFBase* param)
{
  string  temp_title = "";
  parse_string(param, temp_title, "TITLE");
  if(! temp_title.empty())
  {
    my_instructions.title = temp_title;
  }
}

void
parser::parse_region(const LSFBase* param)
{
  string  region_name = "";
  string  first_name  = "";
  string  last_name   = "";
  
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_string(param, region_name, "region");
    }
    
    if(a_list->has_attr("first"))
    {
      parse_string(param, "FIRST", first_name, "first");
    }
    else if(a_list->has_attr("start"))
    {
      parse_string(param, "START", first_name, "start");    
    }
    else if(a_list->has_attr("begin"))
    {
      parse_string(param, "BEGIN", first_name, "begin");        
    }
    
    if(a_list->has_attr("last"))
    {
      parse_string(param, "LAST", last_name, "last");
    }
    else if(a_list->has_attr("end"))
    {
      parse_string(param, "END", last_name, "end");    
    }

    if(! region_name.empty())
    {
      // - Using genome description file.
      //
      if(first_name.empty() && last_name.empty())
      {
        if(my_genome != 0)
        {
          if(my_genome->region(region_name).valid())
          {
            my_instructions.regions.push_back(instructions::region_data(region_name));
          }
          else
          {
            errors << priority(error) << "'" << region_name << "' not a valid genome region.  "
                   << "Ignoring region ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "No genome description file found.  "
                 << "Ignoring region '" << region_name << "' ..." << endl;
        }      
      }
      
      // - Region bounds described in parameter file.
      //
      else
      {
        size_t  first_index = marker_index(first_name);
        size_t  last_index  = marker_index(last_name);
      
        if(first_index != MAX_SIZE &&
           last_index  != MAX_SIZE   )
        {
          my_instructions.regions.push_back(instructions::region_data(region_name, first_index, last_index));
        }
        else
        {
          if(first_index == MAX_SIZE)
          {
            errors << priority(error)  << "'" << first_name << "' not a valid marker name.  "
                   << "Ignoring region '" << region_name << "' ..." << endl;
          }
          
          if(last_index == MAX_SIZE)
          {
            errors << priority(error)  << "'" << last_name << "' not a valid marker name.  "
                   << "Ignoring region '" << region_name << "' ..." << endl;        
          }
        }
      }
    }
    else
    {
      errors << priority(error) << "No name given for parameter, region.  Ignoring "
             << "region ..." << endl;
    }     
  }    
}

size_t
parser::marker_index(const string& marker_name) const
{
  const RPED::RefMPedInfo&  ref_mped_info = my_ref_mped->info();
  
  return  ref_mped_info.marker_find(marker_name);
}

void
parser::parse_epsilon(const LSFBase* param)
{
  double  temp_epsilon = 0.0;
  parse_real(param, temp_epsilon, "EPSILON");
  if(1.0 > temp_epsilon && temp_epsilon > 0.0)
  {
    my_instructions.epsilon = temp_epsilon;
  }        
  else
  {
    errors << priority(error) << "Invalid value given for epsilon.  "
           << "Using default value, " << EPSILON_DEFAULT << "." << endl;
  }
}

void
parser::parse_starting_points(const LSFBase* param)
{
  int  temp_starting_points = 0;
  parse_integer(param, temp_starting_points, "STARTING_POINTS");
  if(temp_starting_points > 0)
  {
    my_instructions.starting_points = temp_starting_points;
  }        
  else
  {
    errors << priority(error) << "Invalid value given for starting points.  "
           << "Using default value, " << STARTING_POINTS_DEFAULT << "." << endl;
  }
}

void  
parser::parse_dump(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.dump, "dump");
    }
    
    if(a_list->has_attr("cutoff"))
    {
      parse_real(param, "CUTOFF", my_instructions.dump_cutoff, "cutoff");
      if(! (0 <= my_instructions.dump_cutoff &&
            my_instructions.dump_cutoff <= 1   ))
      {
        errors << priority(error) << "Value of dump attribute, cutoff, "
               << "is outside of the allowed range [0, 1].  Using default value, " << DUMP_CUTOFF_DEFAULT 
               << "." << endl;
        my_instructions.dump_cutoff = DUMP_CUTOFF_DEFAULT;
      }
    }
  }
}

void
parser::parse_seed(const LSFBase* param)
{
  parse_integer(param, my_instructions.seed, "SEED");
  
  if(my_instructions.seed < 1 || MAX_INT < my_instructions.seed)
  {
    my_instructions.seed = SEED_DEFAULT;
  }
}



void
parser::parse_filters(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "MAF_FILTER")
      {
        parse_maf_filter(*iter);
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "filters sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
  }  
}

void  
parser::parse_maf_filter(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.maf_filter, "maf_filter");
    }
    
    if(a_list->has_attr("threshold"))
    {
      double  temp_maf_threshold = 1;
      parse_real(param, "THRESHOLD", temp_maf_threshold, "threshold");
      if(0 <= temp_maf_threshold && temp_maf_threshold < MAX_MAF_THRESHOLD)
      {
        my_instructions.maf_threshold = temp_maf_threshold;
      }
      else
      {
        errors << priority(error) << "Invalid value given for maf_filter "
               << "attribute, threshold, in filters sub-block.  "
               << "Using default value, " << MAF_THRESHOLD_DEFAULT << " ..." << endl;
        my_instructions.fg_threshold = MAF_THRESHOLD_DEFAULT;
      }
    }
  }  
}

void
parser::parse_blocks(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "SLIDING_WINDOW")
      {
        parse_sliding_window(*iter);
      }
      else if(param_name == "FOUR_GAMETE_RULE")
      {
        parse_four_gamete_rule(*iter);
      }
      else if(param_name == "LD")
      {
        parse_ld(*iter);
      }      
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "blocks sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
  }  
}

void  
parser::parse_sliding_window(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.sliding_window, "sliding_window");
    }
    
    if(a_list->has_attr("width"))
    {
      int  temp_window_width = 0;
      parse_integer(param, "WIDTH", temp_window_width, "width");
      if(temp_window_width > 1)
      {
        my_instructions.window_width = temp_window_width;
      }
      else
      {
        errors << priority(error) << "Invalid value given for sliding_window "
               << "attribute, width, in blocks sub-block.  "
               << "Using default value, " << WINDOW_WIDTH_DEFAULT << " ..." << endl;
        my_instructions.window_width = WINDOW_WIDTH_DEFAULT;
      }
    }
  }  
}

void  
parser::parse_four_gamete_rule(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.four_gamete_rule, "four_gamete_rule");
    }
    
    if(a_list->has_attr("threshold"))
    {
      double  temp_fg_threshold = QNAN;
      parse_real(param, "THRESHOLD", temp_fg_threshold, "threshold");
      if(! isnan(temp_fg_threshold) && 0 <= temp_fg_threshold && temp_fg_threshold < MAX_FG_THRESHOLD)
      {
        my_instructions.fg_threshold = temp_fg_threshold;
      }
      else
      {
        errors << priority(error) << "Invalid value given for four_gamete_rule "
               << "attribute, threshold, in blocks sub-block.  "
               << "Using default value, " << FG_THRESHOLD_DEFAULT << " ..." << endl;
        my_instructions.fg_threshold = FG_THRESHOLD_DEFAULT;
      }
    }
  }  
}

void  
parser::parse_ld(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.ld_blocks, "ld");
    }
    
    if(a_list->has_attr("threshold"))
    {
      double  temp_ld_threshold = QNAN;
      parse_real(param, "THRESHOLD", temp_ld_threshold, "threshold");
      if(! isnan(temp_ld_threshold) && 0 < temp_ld_threshold && temp_ld_threshold < MAX_LD_THRESHOLD)
      {
        my_instructions.ld_threshold = temp_ld_threshold;
      }
      else
      {
        errors << priority(error) << "Invalid value given ld "
               << "attribute, threshold, in blocks sub-block.  "
               << "Using default value, " << LD_THRESHOLD_DEFAULT << " ..." << endl;
        my_instructions.ld_threshold = LD_THRESHOLD_DEFAULT;
      }
    }
  }  
}

void
parser::check_for_mld_file()
{
  if(! my_data.parsed_arguments().argument_specified(APP::LOCUS_FILE))
  {
    if(my_instructions.maf_filter)
    {
      errors << priority(error) << "Parameter, maf_filter, specified in the 'filters' "
             << "sub-block without a marker locus description file.  Ignoring maf_filter ..." << endl;
      my_instructions.maf_filter = false;      
    }
      
    if(my_instructions.ld_blocks)
    {
      errors << priority(error) << "Parameter, ld, specified in the 'blocks' sub-block "
             << "without a marker locus description file.  Ignoring ld ..." << endl;
      my_instructions.ld_blocks = false;
    }
  }
}

void
parser::parse_data(const LSFBase* param)
{
  if(param->List())
  {
    LSFBase*  family_rep_ptr = 0;
    LSFBase*  pools_ptr = 0;
  
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "POOLS")
      {
        pools_ptr = *iter;
      }      
      else if(param_name == "PARTITION")
      {
        parse_partition(*iter);
      }
      else if(param_name == "FAMILY_REP")
      {
        family_rep_ptr = *iter;
      }      
      else if(param_name == "ANALYSIS_UNIT")
      {
        parse_analysis_unit(*iter);
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "data sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
    
    if(family_rep_ptr)
    {
      if(my_instructions.analysis_unit == instructions::FAMILY_REP)
      {    
        parse_family_rep(family_rep_ptr);
      }
      else
      {
        errors << priority(warning) << "Parameter, "
               << "family_rep, given with parameter, analysis_unit, not equal to 'family_rep'.  "
               << "Ignoring parameter, family_rep ..." << endl;          
      }
    }
    
    if(pools_ptr)
    {
      if(my_instructions.analysis_unit == instructions::POOL)
      {
        parse_pools(pools_ptr);
      }
      else
      {
        errors << priority(error) << "Sub-block, "
               << "pools, given with parameter, analysis_unit, not equal to 'pool'.  "
               << "Ignoring sub-block, pools ..." << endl;    
      }    
    }
  }
}

void
parser::parse_pools(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "POOL_SIZE")
      {
        parse_pool_size(*iter);
      }      
      else if(param_name == "POOL_SIZE_TRAIT")
      {
        parse_pool_size_trait(*iter);
      }            
      else if(param_name == "LOCUS")
      {
        parse_locus(*iter);
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "pools sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
  }  
}

void  
parser::parse_pool_size(const LSFBase* param)
{
  int  temp_pool_size = 0;
  parse_integer(param, temp_pool_size, "POOL_SIZE");
  if(temp_pool_size > 0)
  {
    my_instructions.pool_size = temp_pool_size;
  }        
  else
  {
    errors << priority(error) << "Invalid value given for pool_size.  "
           << "Ignoring ..." << endl;
  }
}

void
parser::parse_pool_size_trait(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list && a_list->has_attr("value"))
  {
    string  pool_size_name = "";
    parse_string(param, pool_size_name, false);
    if(! pool_size_name.empty())
    {
      partition_data  parse_results;  // This structure handy for pool_size_trait too!
      get_field_type_and_index(pool_size_name, parse_results);
      switch(parse_results.type)
      {
        case instructions::CONTINUOUS_FIELD:
        {
          map<size_t, trait_info>::iterator  t_iter     = my_trait_registry.find(parse_results.field_index);
          map<size_t, trait_info>::iterator  t_end_iter = my_trait_registry.end();
          if(t_iter != t_end_iter && t_iter->second.usage != FOR_POOL_SIZE)
          {
            errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
                   << "a value for parameter, pool_size_trait, is also given as a value for parameter, "
                   << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, pool_size_trait ..."
                   << endl;
          }
          else
          {        
            my_instructions.pool_size_trait = parse_results.field_index;
            my_trait_registry[parse_results.field_index] = trait_info(pool_size_name, FOR_POOL_SIZE);
          }
        }
            
          break;
          
        case instructions::BINARY_FIELD:
        case instructions::STRING_FIELD:        
          errors << priority(error) << "Value of parameter, pool_size_trait, not the name of a continuous trait, "
                 << "phenotype, covariate.  Ignoring parameter, pool_size_trait ..." << endl;
          break;
          
        case instructions::INVALID_FIELD:
          errors << priority(error) << "Value of parameter, pool_trait_size, not the name of a valid trait, "
                 << "phenotype, covariate or string.  Ignoring parameter, pool_trait_size ..." << endl;
          break;
          
        default:
          assert(false);
      }
    }
    else
    {
      errors << priority(error) << "Bad value given for parameter, pool_size_trait. Ignoring "
             << "parameter, pool_size_trait ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value given for parameter, pool_size_trait.  Ignoring "
           << "parameter, pool_size_trait ..." << endl;
  }      
}

void  
parser::parse_locus(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list && a_list->has_attr("value"))
  {
    string  locus_name = "";
    parse_string(param, locus_name, false);
    if(! locus_name.empty())
    {
      if(param->List())
      {
        pool_locus  pl;
        pl.name = locus_name;
        
        pool_locus::allele  last_allele("");
      
        LSFList::const_iterator  iter;
        for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
        {
          string  param_name = toUpper((*iter)->name());
          if(param_name == "ALLELE")
          {
            parse_allele(*iter, pl);
          }      
          else if(param_name == "LAST_ALLELE")
          {
            parse_last_allele(*iter, pl, last_allele);
          }            
          else
          {
            errors << priority(error) << "Parameter '" << param_name << "' in "
                   << "locus sub-block not recognized.  Ignoring parameter ..." << endl;
          }
        }
        
        // - Add last allele to locus.
        //
        if(last_allele.name.empty() || (! add_allele(pl, last_allele)))
        {
          errors << priority(error) << "Could not find valid last_allele parameter in '" << pl.name 
                 << "' locus sub-block.  Ignoring '" << pl.name << "' locus sub-block ..." << endl;
          
          return;
        }
        
        // - Add locus to instructions.
        //
        size_t  allele_count = pl.alleles.size();
        if(allele_count > 1)
        {
          vector<pool_locus>::iterator  pl_iter;
          vector<pool_locus>::iterator  pl_end_iter = my_instructions.pool_loci.end();
          pl_iter = find(my_instructions.pool_loci.begin(), pl_end_iter, pl);
          
          if(pl_iter == pl_end_iter)
          {
            my_instructions.pool_loci.push_back(pl);
          }
          else
          {
            errors << priority(error) << "Duplicate locus given.  Ignoring duplicate locus ..." 
                   << endl;
          }
        }
        else
        {
          errors << priority(error) << "Sub-block, locus, contains fewer than two valid alleles.  "
                 << "Ignoring sub-block, locus ..." << endl;
        }
      }
      else
      {
        errors << priority(error) << "No alleles given in sub-block, locus.  Ignoring "
               << "sub-block, locus ..." << endl;  
      }      
    }
    else
    {
      errors << priority(error) << "Bad value given for sub-block, locus.  Ignoring "
             << "sub-block, locus ..." << endl;    
    }
  }
  else
  {
    errors << priority(error) << "No value given for sub-block, locus.  Ignoring "
           << "sub-block, locus ..." << endl;  
  }
}      

void  
parser::parse_allele(const LSFBase* param, pool_locus& pl)
{
  AttrList*  a_list = param->attrs();
  if(a_list && a_list->has_attr("value"))
  {
    string  allele_name = "";
    parse_string(param, allele_name, false);
    if(! allele_name.empty())
    {
      pool_locus::allele  a(allele_name);
    
      if(a_list->has_attr("TRAIT"))
      {
        string  trait_name = a_list->StringAttr("TRAIT");
        if(! trait_name.empty())
        {
          partition_data  parse_results;       // This structure handy here too!
          get_field_type_and_index(trait_name, parse_results);
          
          switch(parse_results.type)
          {
            case instructions::CONTINUOUS_FIELD:
              a.index = parse_results.field_index;
              add_allele(pl, a, trait_name);
              break;
              
            case instructions::BINARY_FIELD:
            case instructions::STRING_FIELD:        
            case instructions::INVALID_FIELD:
              errors << priority(error) << "Value attribute, trait, of parameter, allele, "
                     << "not the name of a valid continuous trait, "
                     << "phenotype, covariate or string.  Ignoring parameter, allele ..." << endl;
              break;
              
            default:
              assert(false);
          }
        }
        else
        {
          errors << priority(error) << "No value given for attribute, trait, of parameter, allele. "
                 << "Ignoring parameter, allele ..." << endl;
        }
      }
      else
      {
        errors << priority(error) << "No trait attribute given for parameter, allele "
               << "Ignoring parameter, allele ..." << endl;
      }
    }
    else
    {
      errors << priority(error) << "Bad value given for parameter, allele. Ignoring "
             << "parameter, allele ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value given for parameter, allele.  Ignoring "
           << "parameter, allele ..." << endl;
  }
}

void  
parser::parse_last_allele(const LSFBase* param, const pool_locus& pl, pool_locus::allele& a)
{
  AttrList*  a_list = param->attrs();
  if(a_list && a_list->has_attr("value"))
  {
    string  allele_name = "";
    parse_string(param, allele_name, false);
    if(! allele_name.empty())
    {
      if(a_list->has_attr("TRAIT"))
      {
        errors << priority(error) << "Attribute, trait, not allowed with parameter, last allele, in "
               << "sub-block, locus.  Ignoring attribute, trait ..." << endl;
      }   
              
      a.name = allele_name;       
    }
    else
    {
      errors << priority(error) << "Bad value given for parameter, last_allele. Ignoring "
             << "parameter, last_allele ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value given for parameter, allele_trait.  Ignoring "
           << "parameter, last_allele ..." << endl;
  }
}

bool
parser::add_allele(pool_locus& pl, const pool_locus::allele& a, const string& trait_name)
{
  bool  success = false;

  map<size_t, trait_info>::iterator  t_iter     = my_trait_registry.find(a.index);
  map<size_t, trait_info>::iterator  t_end_iter = my_trait_registry.end();
  if(t_iter != t_end_iter)
  {
    if(t_iter->second.usage != FOR_ALLELE_PCT)
    {
      errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
             << "a value for parameter, allele, is also given as a value for parameter, "
             << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, allele ..."
             << endl;
    }
    else
    {
      errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", used "
             << "used more than once as a value for parameter, allele.  "
             << "Ignoring all but first instance found of parameter, allele, using this value ..."
             << endl;      
    }
  }
  else
  {       
    pair<set<pool_locus::allele>::iterator, bool>  result = pl.alleles.insert(a);
    if(result.second == true)
    {
      success = true;
    }
    else
    {
      errors << priority(error) << "allele name '" << a.name << "' repeated in locus sub-block '"
             << pl.name << "'.  Ignoring second instance of parameter with value '"
             << a.name << "' ..." << endl;
    }
    
    if(a.index != (size_t)(-1))    // last_allele has no trait.
    {
      my_trait_registry[a.index] = trait_info(trait_name, FOR_ALLELE_PCT);
    }
  }
  
  return  success;
}

void  
parser::parse_partition(const LSFBase* param)
{
  assert(my_instructions.partitions.size() == PARTITION_COUNT);
  if(! (my_instructions.partitions[0].valid() && my_instructions.partitions[1].valid()))
  {
    AttrList*  a_list = param->attrs();
    if(a_list && a_list->has_attr("value"))
    {
      string  partition_name = "";
      parse_string(param, partition_name, false);
      if(! partition_name.empty())
      {
        partition_data  parse_results;
        get_field_type_and_index(partition_name, parse_results);
        if(! partition_repeated(parse_results))
        {
          switch(parse_results.type)
          {
            case instructions::CONTINUOUS_FIELD:
              parse_sub_pops(param, parse_results);
              if(parse_results.sub_pops.empty())
              {
                fill_continuous_sub_pops(parse_results);
              }
              
              break;
              
            case instructions::BINARY_FIELD:
              parse_sub_pops(param, parse_results);
              if(parse_results.sub_pops.empty())
              {
                fill_binary_sub_pops(parse_results);
              }
                      
              break;
              
            case instructions::STRING_FIELD:
              parse_sub_pops(param, parse_results);
              if(parse_results.sub_pops.empty())
              {
                fill_string_sub_pops(parse_results);
              }          
              
              break;
              
            case instructions::INVALID_FIELD:
              errors << priority(error) << "Value of parameter, partition, not the name of a valid trait, "
                     << "phenotype, covariate or string.  Ignoring parameter, partition ..." << endl;
              break;
              
            default:
              assert(false);
          }
          
          if(parse_results.type != instructions::INVALID_FIELD)
          {
            // - Check to see that this trait not used elsewhere.  Note:  previous use in
            //   a partition parameter already checked.
            //          
            bool  field_registered = false;
            if(parse_results.type == instructions::CONTINUOUS_FIELD ||
               parse_results.type == instructions::BINARY_FIELD       )
            {
              map<size_t, trait_info>::iterator  t_iter     = my_trait_registry.find(parse_results.field_index);
              map<size_t, trait_info>::iterator  t_end_iter = my_trait_registry.end();
              if(t_iter != t_end_iter && t_iter->second.usage != FOR_PARTITIONING)
              {
                field_registered = true;
                errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
                       << "a value for parameter, partition, is also given as a value for parameter, "
                       << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, partition ..."
                       << endl;          
              }
            }
            else if(parse_results.type == instructions::STRING_FIELD)
            {
              map<size_t, trait_info>::iterator  t_iter     = my_string_registry.find(parse_results.field_index);
              map<size_t, trait_info>::iterator  t_end_iter = my_string_registry.end();
              if(t_iter != t_end_iter && t_iter->second.usage != FOR_PARTITIONING)
              {
                field_registered = true;
                errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
                       << "a value for parameter, partition, is also given as a value for parameter, "
                       << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, partition ..."
                       << endl;          
              }
            }            
            
            if(! field_registered)
            {
              // - Fill second partition first!
              //
              if(! my_instructions.partitions[1].valid())
              {
                my_instructions.partitions[1] = parse_results;
              }
              else
              {
                my_instructions.partitions[0] = parse_results;
              }
              
              if(parse_results.type == instructions::CONTINUOUS_FIELD ||
                 parse_results.type == instructions::BINARY_FIELD       )
              {
                my_trait_registry[parse_results.field_index] = trait_info(partition_name, FOR_PARTITIONING);
              }
              else if(parse_results.type == instructions::STRING_FIELD)
              {
                my_string_registry[parse_results.field_index] = trait_info(partition_name, FOR_PARTITIONING);                
              }
              else
              {
                assert(false);
              }
            }
          }
        }
        else
        {
          errors << priority(error) << "Same pedigree file field specified in two "
                 << "partition sub-blocks.  Ignoring second partition sub-block ... " << endl;
        }
      }
      else
      {
        errors << priority(error) << "Bad value given for parameter, partition. Ignoring "
               << "parameter, partition ..." << endl;
      }
    }
    else
    {
      errors << priority(error) << "No value given for parameter, partition.  Ignoring "
             << "parameter, partition ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "More than two partition sub-blocks data sub-block.  Ignoring all "
           << "but the first two ..." << endl;
  }
}

void  
parser::get_field_type_and_index(const string& partition_name, partition_data& parse_results)
{
  const RPED::RefMPedInfo&  ref_mped_info = my_ref_mped->info();
  size_t  index = ref_mped_info.string_find(partition_name);
  if(index != (size_t)(-1))
  {
    parse_results.type = instructions::STRING_FIELD;
    parse_results.field_index = index;
  }
  else
  {
    index = ref_mped_info.trait_find(partition_name);
    if(index != (size_t)(-1))
    {
      const RPED::RefTraitInfo& trait_info = ref_mped_info.trait_info(index);
      RPED::RefTraitInfo::trait_t  trait_type = trait_info.type();
      switch(trait_type)
      {
        case RPED::RefTraitInfo::continuous_trait:
          parse_results.type = instructions::CONTINUOUS_FIELD;
          parse_results.field_index = index;
          break;
          
        case RPED::RefTraitInfo::binary_trait:
          parse_results.type = instructions::BINARY_FIELD;
          parse_results.field_index = index;
          break;
          
        default:
          parse_results.type = instructions::INVALID_FIELD;          
      }  
    }
    else
    {
      parse_results.type = instructions::INVALID_FIELD;
    }    
  }
}

bool
parser::partition_repeated(const partition_data& parse_results)
{
  assert(my_instructions.partitions.size() == PARTITION_COUNT);

  return  my_instructions.partitions[1].valid() &&
          parse_results.type == my_instructions.partitions[1].type &&
          parse_results.field_index == my_instructions.partitions[1].field_index;
}

void  
parser::parse_sub_pops(const LSFBase* param, partition_data& parse_results)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      
      if(param_name == "SUB_POP")
      {
        switch(parse_results.type)
        {
          case instructions::CONTINUOUS_FIELD:
            parse_continuous_sub_pop(*iter, parse_results);
            break;
            
          case instructions::BINARY_FIELD:
            parse_binary_sub_pop(*iter, parse_results);
            break;          
            
          case instructions::STRING_FIELD:
            parse_string_sub_pop(*iter, parse_results);
            break;
            
          default:
            assert(false);                      
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "partition sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
  }  
}

void  
parser::parse_continuous_sub_pop(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    string  sub_pop_name = "";
    instructions::value  sub_pop_value;    
  
    if(a_list->has_attr("value"))
    {
      parse_string(param, sub_pop_name, true);
    }
    
    if(a_list->has_attr("sub_pop_value"))
    {
      parse_real(param, "SUB_POP_VALUE", sub_pop_value.dbl, false);      
    }

    if(! sub_pop_value.empty())
    {
      if(sub_pop_name.empty())
      {
        sub_pop_name = doub2str(sub_pop_value.dbl);      
      }
      
      pair<string, value>  sub_pop = make_pair(sub_pop_name, sub_pop_value);
      if(new_sub_pop(sub_pop, parse_results))
      {
        parse_results.sub_pops.insert(sub_pop);
      }
      else
      {
        errors << priority(error) << "Non-unique sub_pop name or sub_pop_value given in a partiion "
               << "sub-block.  Arbitrarily ignoring one of the sub_pop parameters with the duplication ..."
               << endl; 
      }
    }
    else
    {
      errors << priority(error) << "Bad or missing value for sub_pop parameter "
             << "attribute, sub_pop_value.  Ignoring sub_pop parameter ..." 
             << endl;
    }
  }
}

void  
parser::parse_binary_sub_pop(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    string  sub_pop_name = "";
    double  real_value = QNAN;
    string  string_value = "";
    instructions::value  final_value;
  
    if(a_list->has_attr("value"))
    {
      parse_string(param, sub_pop_name, true);
    }
    
    if(a_list->has_attr("sub_pop_value"))
    {
      parse_real(param, "SUB_POP_VALUE", real_value, false);
      parse_string(param, "SUB_POP_VALUE", string_value, false);      
    }
    
    const RPED::RefMPedInfo&  ref_mped_info = my_ref_mped->info();
    const RPED::RefTraitInfo& trait_info = ref_mped_info.trait_info(parse_results.field_index);
    
    // - Set value.
    //
    string  value_code = "";
    if(! SAGE::isnan(real_value))
    {
      if(real_value == trait_info.numeric_unaffected_code())
      {
        final_value.dbl = 0;
        value_code = doub2str(real_value);      
      }
      else if(real_value == trait_info.numeric_affected_code())
      {
        final_value.dbl = 1;
        value_code = doub2str(real_value);
      }
    }
    else if(! string_value.empty())
    {
      if(string_value == trait_info.string_unaffected_code())
      {
        final_value.dbl = 0;
        value_code = string_value;
      }
      else if(string_value == trait_info.string_affected_code())
      {
        final_value.dbl = 1;
        value_code = string_value;
      }
    }
    
    if(! final_value.empty())
    {
      // - Set_name if necessary.
      //
      if(sub_pop_name.empty())
      {
        sub_pop_name = value_code;
      }
    
      // - Create subpopulation.
      //
      pair<string, value>  sub_pop = make_pair(sub_pop_name, final_value);
      if(new_sub_pop(sub_pop, parse_results))
      {
        parse_results.sub_pops.insert(sub_pop);
      }
      else
      {
        errors << priority(error) << "Non-unique sub_pop name or sub_pop_value given in a partiion "
               << "sub-block.  Arbitrarily ignoring one of the sub_pop parameters with the duplication ..."
               << endl; 
      }
      
      return;
    }
  }
  
  errors << priority(error) << "Bad or missing value for sub_pop parameter "
         << "attribute, sub_pop_value.  Ignoring sub_pop parameter ..." 
         << endl;
}    

        
void
parser::parse_string_sub_pop(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    string  sub_pop_name = "";
    instructions::value  sub_pop_value;    
  
    if(a_list->has_attr("value"))
    {
      parse_string(param, sub_pop_name, true);
    }
    
    if(a_list->has_attr("sub_pop_value"))
    {
      parse_string(param, "SUB_POP_VALUE", sub_pop_value.str, false);      
    }
    
    if(! sub_pop_value.empty())
    {
      if(sub_pop_name.empty())
      {
        sub_pop_name = sub_pop_value.str;      
      }
      
      pair<string, value>  sub_pop = make_pair(sub_pop_name, sub_pop_value);
      if(new_sub_pop(sub_pop, parse_results))
      {
        parse_results.sub_pops.insert(sub_pop);
      }
      else
      {
      errors << priority(error) << "Non-unique sub_pop name or sub_pop_value given in a partiion "
             << "sub-block.  Arbitrarily ignoring one of the sub_pop parameters with the duplication ..."
             << endl; 
      }
    }
    else
    {
      errors << priority(error) << "Bad or missing value for sub_pop parameter "
             << "attribute, sub_pop_value.  Ignoring sub_pop parameter ..." 
             << endl;
    }
  }
}

// - Return 'true' if value or subpopulation name is not already contained in 
//   parse_results.sub_pops.
//
bool
parser::new_sub_pop(const pair<string, instructions::value>& sub_pop, 
                    const instructions::partition_data& parse_results)
{
  map<string, value>::const_iterator  iter;
  iter = find_if(parse_results.sub_pops.begin(), parse_results.sub_pops.end(), duplication(sub_pop));
  
  return  iter == parse_results.sub_pops.end(); 
}

void  
parser::fill_continuous_sub_pops(partition_data& parse_results)
{
  RefMultiPedigree::pedigree_const_iterator  p_iter     = my_ref_mped->pedigree_begin();
  RefMultiPedigree::pedigree_const_iterator  p_end_iter = my_ref_mped->pedigree_end();    
  for(; p_iter != p_end_iter; ++p_iter)
  {
    const RefPedInfo&  ped_info = p_iter->info();
  
    RefMultiPedigree::member_const_iterator  m_iter     = p_iter->member_begin();
    RefMultiPedigree::member_const_iterator  m_end_iter = p_iter->member_end();
    for(; m_iter != m_end_iter; ++m_iter)
    {
      string               sub_pop_name = "";
      instructions::value  sub_pop_value;  
      
      sub_pop_value.dbl = ped_info.trait(m_iter->index(), parse_results.field_index);
      sub_pop_name = doub2str(sub_pop_value.dbl);
      
      if(! sub_pop_value.empty())
      {
        parse_results.sub_pops.insert(make_pair(sub_pop_name, sub_pop_value));
      }
    }    
  }
}

// - In this case it is not necessary to consult the data as there are, by definition,
//   only two valid values.
//
void  
parser::fill_binary_sub_pops(partition_data& parse_results)
{
  const RPED::RefMPedInfo&  ref_mped_info = my_ref_mped->info();
  const RPED::RefTraitInfo& trait_info = ref_mped_info.trait_info(parse_results.field_index);  
  
  instructions::value  unaffected;
  unaffected.dbl = 0;
  
  instructions::value  affected;
  affected.dbl = 1;
  
  parse_results.sub_pops.insert(make_pair(trait_info.string_unaffected_code(), unaffected));
  parse_results.sub_pops.insert(make_pair(trait_info.string_affected_code(), affected));  
}
        
// - Create a subpopulaton for every distinct non-missing value in the partition
//   variable.
//  
void  
parser::fill_string_sub_pops(partition_data& parse_results)
{
  RefMultiPedigree::pedigree_const_iterator  p_iter     = my_ref_mped->pedigree_begin();
  RefMultiPedigree::pedigree_const_iterator  p_end_iter = my_ref_mped->pedigree_end();    
  for(; p_iter != p_end_iter; ++p_iter)
  {
    const RefPedInfo&  ped_info = p_iter->info();
  
    RefMultiPedigree::member_const_iterator  m_iter     = p_iter->member_begin();
    RefMultiPedigree::member_const_iterator  m_end_iter = p_iter->member_end();
    for(; m_iter != m_end_iter; ++m_iter)
    {
      string               sub_pop_name = "";
      instructions::value  sub_pop_value;  
      
      sub_pop_value.str = ped_info.get_string(m_iter->index(), parse_results.field_index);
      sub_pop_name = sub_pop_value.str;
      
      if(! sub_pop_value.empty())
      {
        parse_results.sub_pops.insert(make_pair(sub_pop_name, sub_pop_value));
      }
    }    
  }
}

void  
parser::parse_family_rep(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list && a_list->has_attr("value"))
  {
    string  family_rep_name = "";
    parse_string(param, family_rep_name, false);
    if(! family_rep_name.empty())
    {
      partition_data  parse_results;  // This structure handy for family_rep too!
      get_field_type_and_index(family_rep_name, parse_results);
      switch(parse_results.type)
      {
        case instructions::CONTINUOUS_FIELD:
          parse_continuous_family_rep_value(param, parse_results);
          break;
          
        case instructions::BINARY_FIELD:
          parse_binary_family_rep_value(param, parse_results);
          break;
          
        case instructions::STRING_FIELD:
          parse_string_family_rep_value(param, parse_results);
          break;
          
        case instructions::INVALID_FIELD:
          errors << priority(error) << "Value of parameter, family_rep, not the name of a valid trait, "
                 << "phenotype, covariate or string.  Ignoring parameter, family_rep ..." << endl;
          break;
          
        default:
          assert(false);
      }
      
      if(parse_results.type != instructions::INVALID_FIELD)
      {
        bool  field_registered = false;
      
        if(parse_results.type == instructions::CONTINUOUS_FIELD ||
           parse_results.type == instructions::BINARY_FIELD       )
        {
          map<size_t, trait_info>::iterator  t_iter     = my_trait_registry.find(parse_results.field_index);
          map<size_t, trait_info>::iterator  t_end_iter = my_trait_registry.end();
          if(t_iter != t_end_iter && t_iter->second.usage != FOR_FAMILY_REP)
          {
            errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
                   << "a value for parameter, family_rep, is also given as a value for parameter, "
                   << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, family_rep ..."
                   << endl;
            field_registered = true;
          }
        }
        else if(parse_results.type == instructions::STRING_FIELD)
        {
          map<size_t, trait_info>::iterator  t_iter     = my_string_registry.find(parse_results.field_index);
          map<size_t, trait_info>::iterator  t_end_iter = my_string_registry.end();
          if(t_iter != t_end_iter && t_iter->second.usage != FOR_FAMILY_REP)
          {
            errors << priority(error) << "Pedigree file field, " << t_iter->second.name << ", specified as "
                   << "a value for parameter, family_rep, is also given as a value for parameter, "
                   << trait_usage_2_string(t_iter->second.usage) << ".  Ignoring parameter, family_rep ..."
                   << endl;
            field_registered = true;
          }        
        } 
 
        if(! field_registered)
        {
          // - Set family rep values.
          //
          my_instructions.rep_field_type = parse_results.type;
          my_instructions.family_rep = parse_results.field_index;
          assert(! parse_results.sub_pops.empty());
          my_instructions.family_rep_value = parse_results.sub_pops[DUMMY_SUB_POP_NAME];
          
          if(parse_results.type == instructions::CONTINUOUS_FIELD ||
             parse_results.type == instructions::BINARY_FIELD       )
          {
            my_trait_registry[parse_results.field_index] = trait_info(family_rep_name, FOR_FAMILY_REP);
          }

          else if(parse_results.type == instructions::STRING_FIELD)
          {
            my_string_registry[parse_results.field_index] = trait_info(family_rep_name, FOR_FAMILY_REP);            
          }
          else
          {
            assert(false);
          }
        }
      }
    }
    else
    {
      errors << priority(error) << "Bad value given for parameter, family_rep. Ignoring "
             << "parameter, family_rep ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value given for parameter, family_rep.  Ignoring "
           << "parameter, family_rep ..." << endl;
  }
}

void  
parser::parse_continuous_family_rep_value(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    instructions::value  family_rep_value;    
    
    if(a_list->has_attr("family_rep_value"))
    {
      parse_real(param, "FAMILY_REP_VALUE", family_rep_value.dbl, false);      
    }
    
    if(! family_rep_value.empty())
    {
      parse_results.sub_pops.insert(make_pair(DUMMY_SUB_POP_NAME, family_rep_value));
    }
    else
    {
      errors << priority(error) << "Bad or missing value for family_rep parameter "
             << "attribute, family_rep_value.  Ignoring family_rep parameter ..." 
             << endl;
      parse_results.reset();
    }
  }
}

void  
parser::parse_binary_family_rep_value(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    double  real_value = QNAN;
    string  string_value = "";
    instructions::value  final_value;        
    
    if(a_list->has_attr("family_rep_value"))
    {
      parse_real(param, "FAMILY_REP_VALUE", real_value, false);
      parse_string(param, "FAMILY_REP_VALUE", string_value, false);            
    }

    const RPED::RefMPedInfo&  ref_mped_info = my_ref_mped->info();
    const RPED::RefTraitInfo& trait_info = ref_mped_info.trait_info(parse_results.field_index);
    
    // - Set value.
    //
    if(! SAGE::isnan(real_value))
    {
      if(real_value == trait_info.numeric_unaffected_code())
      {
        final_value.dbl = 0;
      }
      else if(real_value == trait_info.numeric_affected_code())
      {
        final_value.dbl = 1;
      }
    }
    else if(! string_value.empty())
    {
      if(string_value == trait_info.string_unaffected_code())
      {
        final_value.dbl = 0;
      }
      else if(string_value == trait_info.string_affected_code())
      {
        final_value.dbl = 1;
      }
    }    
    
    if(! final_value.empty())
    {
      parse_results.sub_pops.insert(make_pair(DUMMY_SUB_POP_NAME, final_value));
    }
    else
    {
      errors << priority(error) << "Bad or missing value for family_rep parameter "
             << "attribute, family_rep_value.  Ignoring family_rep parameter ..." 
             << endl;
      parse_results.reset();
    }
  }
}

void  
parser::parse_string_family_rep_value(const LSFBase* param, partition_data& parse_results)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    instructions::value  family_rep_value;    
    
    if(a_list->has_attr("family_rep_value"))
    {
      parse_string(param, "FAMILY_REP_VALUE", family_rep_value.str, false);      
    }
    
    if(! family_rep_value.empty())
    {
      parse_results.sub_pops.insert(make_pair(DUMMY_SUB_POP_NAME, family_rep_value));
    }
    else
    {
      errors << priority(error) << "Bad or missing value for family_rep parameter "
             << "attribute, family_rep_value.  Ignoring family_rep parameter ..." 
             << endl;
      parse_results.reset();
    }
  }
}

void
parser::parse_analysis_unit(const LSFBase* param)
{
  string  temp_unit;
  parse_string(param, temp_unit, "ANALYSIS_UNIT");
  temp_unit = toUpper(temp_unit);
  
  if(temp_unit == "EACH_INDIVIDUAL")
  {
    my_instructions.analysis_unit = instructions::EACH_INDIVIDUAL;
  }
  else if(temp_unit == "FAMILY_REP")
  {
    my_instructions.analysis_unit = instructions::FAMILY_REP;
  }
  else if(temp_unit == "FAMILY_FOUNDERS")
  {
    my_instructions.analysis_unit = instructions::FAMILY_FOUNDERS;
  }
  else if(temp_unit == "POOL")
  {
    my_instructions.analysis_unit = instructions::POOL;
  }
  else
  {
    errors << priority(error) << "Bad or missing value for analysis_unit parameter "
           << "in data sub-block.  Ignoring analysis_unit parameter ..." 
           << endl;
  }
} 

void
parser::parse_tasks(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "POP_FREQ")
      {
        parse_pop_freq(*iter);
      }
      
      /* Does not work when there is missing data.  Not part of the documentation.
      else if(param_name == "ALL_POSSIBLE_DIPLOTYPES")
      {
        parse_all_possible_diplotypes(*iter);
      }
      */
      
      else if(param_name == "ALL_POSSIBLE_DIPLOTYPES_TABLE"   ||
              param_name == "ALL_POSSIBLE_COMBINATIONS_TABLE"   )
      {
        parse_all_possible_diplotypes_table(*iter);
      }      
      else if(param_name == "MOST_LIKELY_DIPLOTYPES"   ||
              param_name == "MOST_LIKELY_COMBINATIONS"   )
      {
        parse_most_likely_diplotypes(*iter);
      }      
      else if(param_name == "ALL_POSSIBLE_HAPLOTYPES")
      {
        parse_all_possible_haplotypes(*iter);
      }            
      else if(param_name == "LIKELIHOOD_RATIO_TEST")
      {
        parse_likelihood_ratio_test(*iter);
      }            
      else if(param_name == "COMPUTE_EMPIRICAL_PVALUE")
      {
        parse_compute_empirical_pvalue(*iter);
      }                  
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "tasks sub-block not recognized.  Ignoring parameter ..." << endl;
      }
    }
  }
}

void  
parser::parse_pop_freq(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.pop_freq, "pop_freq");
    }
    
    if(a_list->has_attr("cutoff"))
    {
      parse_real(param, "CUTOFF", my_instructions.freq_cutoff, "cutoff");
      if(! (0 <= my_instructions.freq_cutoff &&
            my_instructions.freq_cutoff <= 1   ))
      {
        errors << priority(error) << "Value of pop_freq attribute, cutoff, "
               << "is outside of the allowed range [0, 1].  Using default value, " << FREQ_CUTOFF_DEFAULT 
               << "." << endl;
        my_instructions.freq_cutoff = FREQ_CUTOFF_DEFAULT;
      }
    }
  }
}

void  
parser::parse_all_possible_diplotypes(const LSFBase* param)
{
  errors << priority(error) << "all_possible_diplotypes currently disabled.  Use "
         << "all_possible_diplotypes_table for same information in an alternative format." << endl;
  //parse_boolean(param, my_instructions.all_possible_diplotypes, "all_possible_diplotypes");
}

void  
parser::parse_all_possible_diplotypes_table(const LSFBase* param)
{
  string  param_name = toUpper(param->name());
  
  if(param_name == "ALL_POSSIBLE_DIPLOTYPES_TABLE")
  {  
    parse_boolean(param, my_instructions.all_possible_diplotypes_table, "all_possible_diplotypes_table");
  }
  else
  {
    parse_boolean(param, my_instructions.all_possible_diplotypes_table, "all_possible_combinations_table");  
  }
}

void  
parser::parse_most_likely_diplotypes(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      string  param_name = toUpper(param->name());

      if(param_name == "MOST_LIKELY_DIPLOTYPES")
      {    
        parse_boolean(param, my_instructions.most_likely_diplotypes, "most_likely_diplotypes");
      }
      else
      {
        parse_boolean(param, my_instructions.most_likely_diplotypes, "most_likely_combinations");        
      }
    }
    
    if(a_list->has_attr("cutoff"))
    {
      parse_real(param, "CUTOFF", my_instructions.likely_cutoff, "cutoff");
      if(! (0 <= my_instructions.likely_cutoff &&
            my_instructions.likely_cutoff <= 1   ))
      {
        errors << priority(error) << "Value of most_likely_diplotypes attribute, cutoff, "
               << "is outside of the allowed range [0, 1].  Using default value, " 
               << LIKELY_CUTOFF_DEFAULT << "." << endl;
        my_instructions.likely_cutoff = LIKELY_CUTOFF_DEFAULT;
      }
    }
  }  
}

// - This is a hidden feature.  See Trac #1438.
//
void  
parser::parse_all_possible_haplotypes(const LSFBase* param)
{
  string  param_name = toUpper(param->name());
  parse_boolean(param, my_instructions.all_possible_haplotypes, "all_possible_haplotypes");
}

void  
parser::parse_likelihood_ratio_test(const LSFBase* param)
{
  parse_boolean(param, my_instructions.likelihood_ratio_test, "likelihood_ration_test");
}

void  
parser::parse_compute_empirical_pvalue(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.compute_empirical_pvalue, "compute_empirical_pvalue");
    }
    
    if(a_list->has_attr("permutations"))
    {
      double  permutations = 0;
      parse_real(param, "PERMUTATIONS", permutations, "permutations");
      my_instructions.permutations = static_cast<size_t>(permutations);
      if(permutations < 1)
      {
        errors << priority(error) << "Value of compute_empirical_pvalue attribute, permutations, "
               << "is not a positive integer.  Ignoring permutations attribute ... " << endl; 
        my_instructions.permutations = PERMUTATIONS_DEFAULT;
      }
    }
    
    if(a_list->has_attr("min_permutations"))
    {
      double  min_permutations = 0;
      if(my_instructions.permutations != 0)
      {
        /* Hidden feature.
        errors << priority(warning) << compute_empirical_errors attributes, min_permutations and permutations "
               << both specified.  Ignoring min_permutations ..." << endl;
        */
      }
      else
      {
        parse_real(param, "MIN_PERMUTATIONS", min_permutations, "min_permutations");
        my_instructions.min_permutations = static_cast<size_t>(min_permutations);
        if(min_permutations < 1)
        {
          /*
          errors << priority(error) << "Value of compute_empirical_pvalue attribute, min_permutations, "
                 << "not a positive value.  Using default value, " << MIN_PERMUTATIONS_DEFAULT << "." << endl;
          */
          my_instructions.min_permutations = MIN_PERMUTATIONS_DEFAULT;
        }
      }
    }            
    
    if(a_list->has_attr("max_permutations"))
    {
      double  max_permutations = 0;
      if(my_instructions.permutations != 0)
      {
        errors << priority(warning) << "compute_empirical_errors attributes, max_permutations and permutations "
               << "both specified.  Ignoring max_permutations ..." << endl;
      }
      else
      {
        parse_real(param, "MAX_PERMUTATIONS", max_permutations, "max_permutations");
        my_instructions.max_permutations = static_cast<size_t>(max_permutations);
        if(max_permutations < 1)
        {
          errors << priority(error) << "Value of compute_empirical_pvalue attribute, max_permutations, "
                 << "not a positive value.  Using default value, " << MAX_PERMUTATIONS_DEFAULT << "." << endl;
          my_instructions.max_permutations = MAX_PERMUTATIONS_DEFAULT;
        }
      }
    }        
    
    if(a_list->has_attr("width"))
    {
      if(my_instructions.permutations != 0)
      {
        errors << priority(warning) << "compute_empirical_errors attributes, width and permutations "
               << "both specified.  Ignoring width ..." << endl;
      }
      else
      {
        parse_real(param, "WIDTH", my_instructions.width, "width");
        if(my_instructions.width < 0 || my_instructions.width > 1)
        {
          errors << priority(error) << "Value of compute_empirical_pvalue attribute, width, "
                 << "is outside of allowable range of [0, 1].  Using default value, "
                 << WIDTH_DEFAULT << "." << endl;
          my_instructions.width = WIDTH_DEFAULT;
        }
      }
    }    
    
    if(a_list->has_attr("confidence"))
    {
      if(my_instructions.permutations != 0)
      {
        errors << priority(warning) << "compute_empirical_errors attributes, confidence and permutations "
               << "both specified.  Ignoring confidence ..." << endl;      
      }
      else
      {
        parse_real(param, "CONFIDENCE", my_instructions.confidence, "confidence");
        if(my_instructions.confidence < 0 || my_instructions.confidence > 1)
        {
          errors << priority(error) << "Value of compute_empirical_pvalue attribute, confidence, "
                 << "is outside of allowable range of [0, 1].  Using default value, "
                 << CONFIDENCE_DEFAULT << "." << endl;
          my_instructions.confidence = CONFIDENCE_DEFAULT;
        }
      }
    }    
  }  
}


}
}


