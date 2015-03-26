//============================================================================
// File:      parser.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   9/5/2 - created.                         djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              parser 
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "lodlink/parser.h"

using namespace std;

namespace SAGE
{

namespace LODLINK
{

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
parser::parser(const RPED::RefMultiPedigree* mp, const lodlink_data& dat, ostream& m, cerrorstream& errors)
    :  BasicParser(errors), my_ref_mped(mp), my_data(dat), my_messages(m)
{
  my_messages << endl << "Parsing LODLINK analyses ...\n" << endl;

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

// - Top level.  Iterate through parameters of lodlink_analysis
//   block.  Client code must insure that params->name() is 'lodlink'
//   or 'lodlink_analysis'.
//
void  
parser::parse_test_parameter_section(const LSFBase* params)
{
  assert(params && my_ref_mped);

  reset();
  init_parse();
  print_header();
  
  parse_out_attr(params);

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
      else if(parameter == "MODEL")
      {
        parse_model(*iter);
      }
      else if(parameter == "LINKAGE_TESTS")
      {
        parse_linkage_tests(*iter);
      }
      else if(parameter == "HOMOG_TESTS")
      {
        parse_homog_tests(*iter);
      }
      else if(parameter == "LODS")
      {
        parse_lods(*iter);
      }
      else if(parameter == "GENOTYPES")
      {
        parse_genotypes(*iter);
      }
      else
      {
        errors << priority(error) << "parameter '" << parameter << "' not recognized.  "
               << "Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "No parameters read for Lodlink analysis block." << endl;  
  }
  
  if(my_instructions.valid)
  {
    check_trait_or_marker();
  }
  
  if(my_instructions.valid)
  {
    check_sf_linkage_type_compatibility();
  }
  
  if(my_instructions.valid)
  {
    check_trait_read_in();
  }
  
  if(my_instructions.valid)
  {
    check_genotype_model_compatibility();
  }  

  print_footer();
}

void
parser::reset()
{
  my_instructions.reset();
}

void parser::init_parse()
{


  // Create the analysis identifier string.
  //
  ++my_analysis_id;
  string  id = long2str(my_analysis_id);

  // Set the title and output to defaults.
  //
  my_instructions.title          = "LODLINK Analysis " + id;
  my_instructions.file_name_root = "lodlink_analysis" + id; 
}

void parser::print_header()
{
    my_messages << "\nParsing new analysis block ...\n" << endl;
}

void parser::print_footer()
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
parser::parse_model(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  
  if(a_list)
  {
    if(a_list->has_attr("TRAIT"))
    {
      parse_trait_attr(param);
    }
    
    if(my_instructions.trait.empty() && a_list->has_attr("MARKER"))
    {
      parse_marker_attr(param);
    }
  }
}

bool  
parser::parse_marker_attr(const LSFBase* param)
{
  bool  success = false;

  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    string  temp_marker = a_list->StringAttr("MARKER");
    if(my_ref_mped->info().marker_exists(temp_marker))
    {
      my_instructions.linkage = instructions::MARKER;
      my_instructions.trait = temp_marker;
      success = true;
    }
    else
    {
      errors << priority(error) << "'" << temp_marker << "' not a valid "
             << "marker.  Ignoring ..." << endl;
    }
  }
  
  return  success;
}

bool  
parser::parse_trait_attr(const LSFBase* param)
{
  bool  success = false;

  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    string  temp_trait = a_list->StringAttr("TRAIT");
    
    // - Trait is read in and stored as a marker.
    //
    if(my_ref_mped->info().marker_exists(temp_trait))
    {
      my_instructions.linkage = instructions::TRAIT;
      my_instructions.trait = temp_trait;
      success = true;
    }
    else
    {
      errors << priority(error) << "No trait locus description read for trait, " << temp_trait 
             << ".  Ignoring ..." << endl;
    }
  }
  
  return  success;  
}

void
parser::parse_linkage_tests(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.linkage_test, "linkage_tests");
    }
    
    if(a_list->has_attr("sex_specific"))
    {
      parse_boolean(param, "SEX_SPECIFIC", my_instructions.linkage_sex_specific,
                    "sex_specific");
    }
    
    if(a_list->has_attr("homog"))
    {
      parse_boolean(param,"HOMOG", my_instructions.linkage_homog, "homog");
    }
  }
}

void
parser::parse_homog_tests(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "SMITHS_TEST")
      {
        parse_smiths_test(*iter);
      }
      else if(param_name == "MORTONS_TEST")
      {
        parse_mortons_test(*iter);
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "homog_tests sub-block not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "Sub-block, homog_tests, contains no parameters."
           << endl;
  }
}

void
parser::parse_smiths_test(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.smiths_test, "smiths_test");
    }
    
    if(a_list->has_attr("sex_specific"))
    {
      parse_boolean(param, "SEX_SPECIFIC", my_instructions.smiths_sex_specific,
                    "sex_specific");
    }
  }
}

void
parser::parse_mortons_test(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.mortons_test, "mortons_test");
    }
    
    if(a_list->has_attr("sex_specific"))
    {
      parse_boolean(param, "SEX_SPECIFIC", my_instructions.mortons_sex_specific,
                    "sex_specific");
    }
  }
  
  // - Parse group information.
  //
  bool  groups_valid = true;
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "GROUP")
      {
        if(groups_valid)
        {
          groups_valid = parse_group(*iter);
        }
        else
        {
          return;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "mortons_test sub-block not recognized.  Ignoring ..." << endl;
      }
    }
    
    if(groups_valid)
    {
      if(! groups_comprehensive())
      {
        my_instructions.groups.clear();
        errors << priority(error) << "Not all pedigrees included in groups.  "
               << "Skipping Morton's test ..." << endl;
        
        return;
      }
    }
  }
  
  // - If no group information is given, each pedigree constitutes a group.
  //
  if(my_instructions.mortons_test == true && 
     my_instructions.groups.empty() &&
     groups_valid)
  {
    build_default_groups();
  }
}

bool
parser::parse_group(const LSFBase* param)
{
  bool  groups_valid = true;
  string  group_name = "";
  parse_string(param, group_name, "GROUP");
  
  if(! group_name.empty())
  {
    group  pedigree_ids;
    parse_pedigree_ids(param, pedigree_ids);
    
    if(! pedigree_ids.empty())
    {
      verify_pedigrees(pedigree_ids);
      
      if(! pedigree_ids.empty())
      {
        groups_valid = insert_group(group_name, pedigree_ids);
      }
      else
      {
        errors << priority(error) << "No valid pedigree ids given for group '"
               << group_name << "'.  Ignoring group ..." << endl;
      }
    } 
    else
    {
      errors << priority(error) << "No pedigree ids given for group '"
             << group_name << "'.  Ignoring group ..." << endl;
    }  
  }
  else
  {
    errors << priority(error) << "Group name missing.  Ignoring group ..." << endl;
  }
  
  return  groups_valid;
}

void
parser::parse_pedigree_ids(const LSFBase* param, group& ids)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "PEDIGREE_ID")
      {
        string  temp_id = "";
        parse_string(*iter, temp_id, "PEDIGREE_ID");
        if(! temp_id.empty())
        {
          if(! ids.insert(temp_id).second)
          {
            errors << priority(error) << "Duplicate pedigree_id given in group "
                   << "sub-block.  Ignoring duplicate pedigree_id ..." << endl;
          }
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << "group sub-block not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(error) << "No parameters given in group sub-block.  "
           << "Ignoring sub-block ..." << endl;
  }
}

// - Validate the pedigree ids and add corresponding pedigree
//   pointers to the group.
//
void  
parser::verify_pedigrees(group& ids)
{
  group  temp_ids;
  set<string>::const_iterator  iter;
  for(iter = ids.begin(); iter != ids.end(); ++iter)
  {
    RPED::RefMultiPedigree::pedigree_const_pointer  ped_pointer = my_ref_mped->pedigree_find(*iter);
    if(ped_pointer)
    {
      temp_ids.insert(*iter);
    }
    else
    {
      errors << priority(error) << "Invalid pedigree_id '" << *iter << "' given "
             << "in group sub-block.  Ignoring ..." << endl;
    }
  }
  
  ids = temp_ids;  
}

bool
parser::insert_group(const string& name, const group& g)
{
  bool  groups_valid = true;
  if(! pedigrees_duplicated(g))
  {
    if(! my_instructions.groups.insert(make_pair(name, g)).second)
    {
      errors << priority(error) << "Duplicate group name '" << name << "' given in "
             << "mortons test sub-block.  Ignoring second group ..." << endl;
    }
  }
  else
  {
    groups_valid = false;
    my_instructions.groups.clear();
    errors << priority(error) << "Pedigree included in more than one group.  Skipping "
           << "Morton's test ..." << endl;
  }
  
  return  groups_valid;
}

// - Is a pedigree included in more than one group?
//
bool
parser::pedigrees_duplicated(const group& g) const
{
  bool  duplicate = false;
  group::const_iterator  p_iter;
  for(p_iter = g.begin(); p_iter != g.end(); ++p_iter)
  {
    map<string, group>::const_iterator  g_iter;
    for(g_iter = my_instructions.groups.begin(); 
        g_iter != my_instructions.groups.end();
        ++g_iter)
    {
      if(g_iter->second.find(*p_iter) != g_iter->second.end())
      {
        duplicate = true;
        break;
      }
    }
  }
  
  return  duplicate;
}

// - Are all pedigrees assigned to a group?
//
// - Assumes check for duplicated pedigrees has all ready been done.
//
bool
parser::groups_comprehensive() const
{
  size_t  pedigree_count = 0;
  
  map<string, group>::const_iterator  g_iter;
  for(g_iter = my_instructions.groups.begin(); 
      g_iter != my_instructions.groups.end();
      ++g_iter)
  {
    pedigree_count += g_iter->second.size();
  }
  
  return  pedigree_count == my_ref_mped->pedigree_count();
}

void
parser::build_default_groups()
{
  RPED::RefMultiPedigree::pedigree_const_iterator  iter;
  for(iter = my_ref_mped->pedigree_begin(); iter != my_ref_mped->pedigree_end(); ++iter)
  {
    group  temp_group;
    temp_group.insert(iter->name());
    my_instructions.groups.insert(make_pair(iter->name(), temp_group));
  }
}

void
parser::parse_lods(const LSFBase* param)
{
  lods_ptrs  pointers;
  
  init_lods_ptrs(param, pointers);
  
  string  option = "standard";
  if(pointers.option)
  {
    string  temp_option = parse_lods_option(pointers.option);
    if(! temp_option.empty())
    {
      option = temp_option;
    }
  }
  
  bool  sex_specific = false;
  if(pointers.sex_specific)
  {
    parse_boolean(pointers.sex_specific, sex_specific, "SEX_SPECIFIC");
  }
  
  option = toUpper(option);
  if(option == "NONE")
  {
    parse_lods_none(pointers);
  }
  else if(option == "STANDARD")
  {
    parse_lods_standard(pointers, sex_specific);
  }
  else if(option == "SPECIFIED")
  {
    parse_lods_specified(pointers, sex_specific);
  }
  else
  {
    assert(false);
  }
}

void
parser::init_lods_ptrs(const LSFBase* param, lods_ptrs& pointers)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        pointers.option = *iter;
      }
      else if(param_name == "SEX_SPECIFIC")
      {
        pointers.sex_specific = *iter;
      }
      else if(param_name == "MALE_FEMALE")
      {
        pointers.male_female = *iter;
      }
      else if(param_name == "AVERAGE")
      {
        pointers.average = *iter;
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in lods sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "Sub-block, lods, contains no parameters."
           << endl;
  }
}

string
parser::parse_lods_option(const LSFBase* param)
{
  string option = "";
  parse_string(param, option, "OPTION");
  option = toUpper(option);
  if(option == "NONE"      || option == "STANDARD" || 
     option == "SPECIFIED" || option == ""           )
  {
    return option;
  }
  else
  { 
    errors << priority(error) << "Value '" << option << "' for "
           << "option parameter of lods sub-block not recognized.  "
           << "Ignoring ..." << endl;
    return "";
  }
}

void
parser::parse_lods_none(const lods_ptrs& pointers)
{
  my_instructions.male_female_thetas.resize(0);
  my_instructions.average_thetas.resize(0);
  
  if(pointers.sex_specific || pointers.male_female ||
     pointers.average                                )
  {
    errors << priority(error) << "Parameter(s) specified in lods sub-block "
           << "with 'none' option.  Ignoring parameter(s) ..." << endl;
  }
}

void
parser::parse_lods_standard(const lods_ptrs& pointers, bool sex_specific)
{
  if(sex_specific)
  {
    my_instructions.male_female_thetas  = vector<theta_pair>(FIRST_THETA_PAIR,
                                                               LAST_THETA_PAIR  );
    my_instructions.average_thetas.resize(0);
  }
  else
  {
    my_instructions.male_female_thetas.resize(0);
    my_instructions.average_thetas = vector<double>(FIRST_THETA, LAST_THETA);
  }
  
  if(pointers.male_female || pointers.average)
  {
    errors << priority(error) << "male_female or average specified in lods "
           << "sub-block with 'standard' option.  Ignoring parameter ... " << endl;
  }
}

void
parser::parse_lods_specified(const lods_ptrs& pointers, bool sex_specific)
{
  my_instructions.average_thetas.resize(0);
  my_instructions.male_female_thetas.resize(0);

  if(sex_specific)
  {
    if(pointers.male_female)
    {
      parse_lods_theta_pairs(pointers.male_female);
      
      if(my_instructions.male_female_thetas.empty())
      {
        errors << priority(error) << "No valid thetas specified in lods sub-block "
               << "with option 'specified' chosen.  Ignoring lods sub-block ..." << endl;
      }
      else
      {
        sort(my_instructions.male_female_thetas.begin(), my_instructions.male_female_thetas.end());
      }
    }
    else
    {
      errors << priority(error) << "No valid male_female thetas specified "
             << "in lods sub-block with option 'specified' chosen and sex_specific "
             << "equal to 'true'.  Ignoring lods sub-block ..." << endl;
    }
  }
  else
  {
    if(pointers.average)
    {
      parse_lods_thetas(pointers.average);
      
      if(my_instructions.average_thetas.empty())
      {
        errors << priority(error) << "No valid average thetas specified in lods sub-block "
               << "with option 'specified' chosen and sex_specific equal to "
               << "false.  Ignoring lods sub-block ..." << endl;
      }
      else
      {
        sort(my_instructions.average_thetas.begin(), my_instructions.average_thetas.end());
      }
    }
    else
    {
      errors << priority(error) << "average thetas must be specified "
             << "in lods sub-block with option 'specified' chosen and sex_specific "
             << "equal to 'false'.  Ignoring lods sub-block ..." << endl; 
    }
  }
}

void
parser::parse_lods_thetas(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "THETA")
      {
        double  temp_theta = QNAN;
        parse_real(*iter, temp_theta, "THETA");
        if(theta_valid(temp_theta) && not_a_duplicate(temp_theta))
        {
          my_instructions.average_thetas.push_back(temp_theta);
        }        
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in lods sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
}

bool
parser::theta_valid(double t)
{
  if(THETA_LOWER_BOUND <= t && t <= THETA_UPPER_BOUND)
  {
    return true;
  }
  else
  {
    errors << priority(error) << t << " not a valid value for theta.  "
           << "Ignoring ..." << endl;
    return false;  
  }
}

bool
parser::not_a_duplicate(double t) const
{
  return find(my_instructions.average_thetas.begin(), 
              my_instructions.average_thetas.end(), t) 
              == my_instructions.average_thetas.end();
}

void
parser::parse_lods_theta_pairs(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "THETA")
      {
        theta_pair  temp_theta_pair;
        parse_theta_pair(*iter, temp_theta_pair);
        if(theta_pair_valid(temp_theta_pair) && 
           not_a_duplicate(temp_theta_pair))
        {
          my_instructions.male_female_thetas.push_back(temp_theta_pair);
        }        
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in lods sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
}

void  
parser::parse_theta_pair(const LSFBase* param, theta_pair& p)
{
  AttrList* a_list = param->attrs();
  
  if(a_list)
  {
    p.male_theta = a_list->FloatAttr("MALE");
    if(p.male_theta == QNAN)
    {
      errors << priority(error) << "Bad or missing value given for "
             << "attribute, male, in lods sub-block.  "
             << "Ignoring attribute ..." << endl;
    }
  
    p.female_theta = a_list->FloatAttr("FEMALE");
    if(p.female_theta == QNAN)
    {
      errors << priority(error) << "Bad or missing value given for "
             << "attribute, female, in lods sub-block.  "
             << "Ignoring attribute ..." << endl;
    }
  }
  else
  {
    errors << priority(error) << "Attributes missing for theta parameter "
           << "in lods sub-block." << endl;
  }
}

bool
parser::theta_pair_valid(theta_pair p)
{
  return theta_valid(p.male_theta) && theta_valid(p.female_theta);
}

bool
parser::not_a_duplicate(theta_pair p) const
{
  return find(my_instructions.male_female_thetas.begin(), 
              my_instructions.male_female_thetas.end(), p) 
              == my_instructions.male_female_thetas.end();
}

void
parser::parse_genotypes(const LSFBase* param)
{
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("value"))
    {
      parse_boolean(param, my_instructions.genotypes, "genotypes");
    }
    
    if(a_list->has_attr("sex_specific"))
    {
      parse_boolean(param, "SEX_SPECIFIC", my_instructions.genotypes_sex_specific,
                    "sex_specific");
    }
  }
}


// - If model is 'marker_marker', a valid marker must be supplied.  If
//   model is 'trait', a valid trait must be supplied.
//
void
parser::check_trait_or_marker()
{
  if(my_instructions.trait.empty())
  {
      errors << priority(critical) << "No valid trait or marker locus specified "
             << "for analysis.  Skipping analysis ..." << endl;
      my_instructions.valid = false;
  }
}

void  
parser::check_sf_linkage_type_compatibility()
{
  if(my_instructions.linkage == instructions::MARKER)
  {
    if(my_instructions.smiths_test)
    {
      errors << priority(error) << "Marker to marker linkage not compatible "
             << "with Smith's test.  Not performing Smith's test." << endl;
      my_instructions.smiths_test = false;
    }
    
    if(my_instructions.linkage_test && (! my_instructions.linkage_homog))
    {
      errors << priority(error) << "Marker to marker linkage not compatible "
             << "with Smith's heterogeneity model.  Not performing linkage test." << endl;      
      my_instructions.linkage_test = false;
    }
  }
}

// - Issue warning if linkage type is trait and the 'primary trait' is a marker
//   as opposed to a model produced by SEGREG and read in.
//
void
parser::check_trait_read_in()
{
  if(my_instructions.linkage == instructions::TRAIT)
  {
    size_t  trait_index = my_ref_mped->info().marker_find(my_instructions.trait);
    if(trait_index < my_data.true_marker_count())
    {
      errors << priority(warning) << "Trait, " << my_instructions.trait << ", is a "
             << " marker, not a trait locus produced by segregation analysis." << endl;
    }
  }
}

// - Genotype probabilities not allowed for marker to marker linkage.
//
void  
parser::check_genotype_model_compatibility()
{
  if(my_instructions.genotypes &&
     my_instructions.linkage != instructions::TRAIT)
  {
    errors << priority(error) << "Marker to marker linkage not compatible "
           << "with calculation of genotype probabilities.  Not calculating "            
           << "genotype probabilities ..." << endl;
    my_instructions.genotypes = false;
  }
}

}
}


