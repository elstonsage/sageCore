//============================================================================
// File:      parser.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   5/15/01 - created.                         djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              parser 
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/parser.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{

// We turn off these lint messages over the whole file.  They are turned
// back on at the end.
//lint -esym(534,*parse_string*,*parse_boolean*)
//lint -e774

// Special defines to make life somewhat easier.  Var capitalized since
// there were name conflicts.  These are undefined at the end of this file
#define mean CovariateSubmodel::ct_MEAN
#define VAR  CovariateSubmodel::ct_VAR
#define susc CovariateSubmodel::ct_SUSC
#define comp CovariateSubmodel::ct_COMP

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
parser::parser(const RPED::RefMultiPedigree* mp, ostream& m, cerrorstream& err)
    : BasicParser(err), my_mp(mp), messages(m)
{
  messages << endl << "Parsing Segreg Analyses..." << endl << endl;

  reset();
  init_genotype_map();

  idnum = 0;
}

void
parser::parse_symbols(const SymbolTable* )
{}

void 
parser::parse_parameter(const LSFBase* )
{}

// - Top level.  Iterate through parameters of segreg_analysis
//   block.  Client code must insure that params->name() is 'segreg'
//   or 'segreg_analysis'.
//   Parsing is done in two passes to control the order in which
//   sub-blocks are parsed to better deal w. dependencies among
//   them.
//
void  
parser::parse_test_parameter_section(const LSFBase* params)
{
  zero_poly_loci = false; // due to JA
  type_var_one = true;  // due to JA
  like_cutoff  = numeric_limits<double>::quiet_NaN();

  reset();
  
  assert(params && my_mp);
  assert(params->List() != 0);
  
  init_parse();

  parse_output_attribute(params);

  // - First parsing pass.  Identify and remember pointers to top-level
  //   parameters.
  //
  LSFList::const_iterator  iter;
  for(iter = params->List()->begin(); 
      iter != params->List()->end() && model_valid(); 
      ++iter)
  {
    if(*iter) 
    {
      classify_parameter(*iter);
    }
  }

  // Parse each sub model in order
  
  second_pass_parse();
  

  // Check constraints between sub models
  check_meta_constraints();

  // Print the footer message
  print_footer();
}

void parser::init_parse()
{
  messages << "Beginning new analysis block...." << endl << endl;

  // Create the analysis identifier string

  ++idnum;
  string id = long2str(idnum);

  // Set the title and output to defaults

  my_model.title          = "SEGREG Analysis " + id;
  my_model.file_name_root = "segreg_analysis" + id; 
}

void parser::print_title()
{
  messages << endl
           << "Parsing Analysis: '" << my_model.title << "'...." 
           << endl << endl;
}

void parser::print_footer()
{
  messages << endl
           << "Analysis parsing complete.  ";

  if(model_valid())
    messages << "Analysis valid." << endl;
  else
    messages << "Analysis is invalid.  Check your file." << endl;

  messages << "------------------------------------------------------------------------------"
           << endl << endl;
}

void  
parser::parameter_duplicated(const string& param_name)
{
  errors << priority(critical) << "Parameter, " << param_name << ", specified "
         << "more than once.  Skipping analysis ..." << endl;
  my_model.m_class = model_INVALID; 
}

void  
parser::parse_test_parameter(const LSFBase* )
{
  // This is a virtual function in the base class.
}



// - Identify and remember pointers to top level
//   parameters.
//
void
parser::classify_parameter(const LSFBase* param)
{
  string param_name = toUpper(param->name());

  if(param_name == "TITLE")
  {
    set_parameter(param_name, &ptrs.title, param);
  }
  else if(param_name == "OUTPUT")
  {
    set_parameter(param_name, &ptrs.output, param);
  }
  else if(param_name == "OUTPUT_OPTIONS")
  {
    set_parameter(param_name, &ptrs.output_options, param);
  }
  else if(param_name == "TRAIT")
  {
    set_parameter(param_name, &ptrs.trait, param);
  }
  else if(param_name == "TYPE_MEAN")
  {
    set_parameter(param_name, &ptrs.type_mean, param);
  }
  else if(param_name == "TYPE_SUSCEPT")
  {
    set_parameter(param_name, &ptrs.type_suscept, param);
  }
  else if(param_name == "TYPE_VAR")
  {
    set_parameter(param_name, &ptrs.type_var, param);
  }
  else if(param_name == "GENO_FREQ")
  {
    set_parameter(param_name, &ptrs.geno_freq, param);
  }
  else if(param_name == "RESID")
  {
    set_parameter(param_name, &ptrs.resid, param);
  }
  else if(param_name == "TRANSMISSION")
  {
    my_model.trans_missing = false;
  
    set_parameter(param_name, &ptrs.transmission, param);
  }
  else if(param_name == "TRANSFORMATION")
  {
    set_parameter(param_name, &ptrs.transformation, param);
  }
  else if(param_name == "MEAN_COV")
  {
    set_parameter(param_name, &ptrs.mean_cov, param);
  }
  else if(param_name == "VAR_COV")
  {
    set_parameter(param_name, &ptrs.var_cov, param);
  }
  else if(param_name == "SUSCEPT_COV")
  {
    set_parameter(param_name, &ptrs.suscept_cov, param);
  }
  else if(param_name == "COMPOSITE_TRAIT")
  {
    set_parameter(param_name, &ptrs.composite_trait, param);
  }
  else if(param_name == "FPMM")
  {
    set_parameter(param_name, &ptrs.fpmm, param);
  }
  else if(param_name == "ASCERTAINMENT")
  {
    set_parameter(param_name, &ptrs.ascertainment, param);
  }
  else if(param_name == "PREV_CONSTRAINT")
  {
    set_parameter(param_name, &ptrs.prev_constraint, param);
  }
  else if(param_name == "PREV_CONSTRAINTS")
  {
    set_parameter(param_name, &ptrs.prev_constraints, param);
  }
  else if(param_name == "PREV_ESTIMATE")
  {
    set_parameter(param_name, &ptrs.prev_estimate, param);
  }
  else if(param_name == "CLASS")
  {
    set_parameter(param_name, &ptrs.class_, param);
  }
  else if(param_name == "MAXFUN")
  {
    set_parameter(param_name, &ptrs.maxfun_options, param);
  }
  else
  {
    errors << priority(error) << "Parameter '" << param_name
           << "' not recognized.  Ignoring ..." << endl;
  }
}

// - Perform actual parsing in an order that takes into account
//   dependencies between sub-models and checks inter-sub-model
//   constraints.
//
void
parser::second_pass_parse()
{
  // TITLE
  if(process_block(ptrs.title))
  {
    parse_title(ptrs.title);
  }
  
  // OUTPUT PARAMETER (for backward compatibility)
  if(process_block(ptrs.output))
  {
    parse_output_parameter(ptrs.output);
  }
  
  
  print_title();
  
  // PRIMARY TRAIT
  if(model_valid())
  {
    if(ptrs.trait)
    {
      parse_primary_trait(ptrs.trait);
    }
    else
    {
      no_primary_trait();
    }
  }

  // OUTPUT OPTIONS -- Requires primary trait type for some output options,
  //   so appears here rather than before the primary trait.
  if(process_block(ptrs.output_options))
  {
    parse_output_options(ptrs.output_options);
  }
  
  // Now we know output, so we can configure the maxfun options.  This
  // is basically independent of model options.
  
  if(process_block(ptrs.maxfun_options))
  {
    string fname = my_model.get_file_name_root() + ".max";
    
    MAXFUN::parseDebugParams(my_model.my_maxfun_debug, ptrs.maxfun_options, errors);
    
    my_model.my_maxfun_debug.setDebugOutput(fname, true);
  }

  // TYPE MEAN
  if(process_block(ptrs.type_mean))
  {
    if(my_model.get_primary_trait_type() != pt_BINARY)
    {
      my_model.mean_missing = false;

      parse_mean_sub_model(ptrs.type_mean, "mean");
    }
    else
    {
      not_relevant_error_message("type_mean sub-block", pt_BINARY);
    }
  }

  // TYPE SUSCEPT
  if(process_block(ptrs.type_suscept))
  {
    if(my_model.get_primary_trait_type() != pt_CONTINUOUS)
    {
        my_model.susc_missing = false;

        parse_type_suscept_sub_model(ptrs.type_suscept);
    }
    else
    {
      not_relevant_error_message("type_suscept sub-block", pt_CONTINUOUS);
    }
  }

  // COMPOSITE TRAIT
  if(process_block(ptrs.composite_trait))
  {
    if(my_model.primary_trait_type == pt_BINARY ||
       my_model.primary_trait_type == pt_ONSET )
    {
      not_relevant_error_message("composite_trait sub-block",
                                 my_model.primary_trait_type);
    }
    else
    {
      parse_cov_sub_model(ptrs.composite_trait, comp);
    }
  }
  
  // MEAN COVARIATE
  if(process_block(ptrs.mean_cov))
  {
    if(my_model.primary_trait_type == pt_BINARY)
    {
      not_relevant_error_message("mean_cov sub-block", pt_BINARY);
    }
    else
    {
      parse_cov_sub_model(ptrs.mean_cov, mean);
    }
  }
  
  // VARIANCE COVARIATE
  if(process_block(ptrs.var_cov))
  {
    if(my_model.primary_trait_type == pt_BINARY)
    {
      not_relevant_error_message("var_cov sub-block", pt_BINARY);
    }
    else
    {
      parse_cov_sub_model(ptrs.var_cov, VAR);
    }
  }
  
  // SUSCEPTIBILITY COVARIATE
  if(process_block(ptrs.suscept_cov))
  {
    if(my_model.primary_trait_type == pt_CONTINUOUS)
    {
      not_relevant_error_message("suscept_cov sub-block", pt_CONTINUOUS);
    }
    else
    {
      parse_cov_sub_model(ptrs.suscept_cov, susc);
    }
  }
  
  // CLASS
  if(model_valid())
  {
    // Set the default model class.  This may be altered in
    // parse_model_class()

    if(my_model.primary_trait_type == pt_BINARY)
    {
      my_model.m_class = model_MLM;
    }

    if(ptrs.class_)
    {
      parse_model_class(ptrs.class_);
    }
  }
  
  // FPMM
  if(model_valid())
  {
    if(ptrs.fpmm)
    {
      parse_fpmm_sub_model(ptrs.fpmm);
      fpmm_parsed = true;
    }
    else
    {
      if(my_model.primary_trait_type == pt_ONSET)
      {
        errors << priority(critical) << "No fpmm sub-block when onset model specified."
               << " Analysis invalid. Skipping analysis." << endl;
        my_model.m_class = model_INVALID;
      }
    }
  }

  // - One type requires an alternative default for the genotype frequency
  //   model, but we don't know which type is primary (mean or susc) until
  //   we have parsed the onset sub block (done in FPMM).
  //
  if(!my_model.get_type_missing() && 
     my_model.type_dependent_sub_model().option() == genotype_specific_mean_sub_model::one)
  {
    my_model.fix_freq_not_none_one_type();
  }
  
  // RESID
  if(model_valid())
  {
    if(ptrs.resid)
    {
      parse_resid_sub_model(ptrs.resid);
      resid_parsed = true;
    }
    else
    {
      //lint -e{534}
      my_model.resid_sub_model.set_as_default(my_model.m_class);
    }
  }
  
  // TRANSFORMATION
  
  // By default, there is no transformation for binary traits.
  if(my_model.primary_trait_type == pt_BINARY)
     my_model.transf_sub_model.set_option(MAXFUN::TransformationSubmodel::no_trans);

  if(process_block(ptrs.transformation))
  {
    parse_transformation_sub_model(ptrs.transformation);
    transformation_parsed = true;
  }
  
  // GENO FREQUENCY
  if(process_block(ptrs.geno_freq))
  {
    parse_freq_sub_model(ptrs.geno_freq);
  }
  
  // TRANSMISSION
  if(process_block(ptrs.transmission))
  {
    parse_transmission_sub_model(ptrs.transmission);
    transmission_parsed = true;
  }

  // VARIANCE
  if(process_block(ptrs.type_var))
  {
    if(my_model.primary_trait_type == pt_BINARY)
    {
      not_relevant_error_message("type_var sub-block", pt_BINARY);
    }
    else
    {
      parse_var_sub_model(ptrs.type_var);
    }
  }
  
  // ASCERTAINMENT
  if(process_block(ptrs.ascertainment))
  {
    parse_ascertainment_sub_model(ptrs.ascertainment);
  }

  // PREVALENCES

  // We can now set the prevalence fpmm and onset options.

  my_model.prev_sub_model.set_fpmm_option(my_model.m_class == model_FPMM);

  my_model.prev_sub_model.set_onset_option(my_model.primary_trait_type == pt_ONSET);
  
  // PREV_CONSTRAINT
  if(process_block(ptrs.prev_constraint))
  {
    if(!ptrs.prev_constraints)
    {
      if(my_model.primary_trait_type == pt_CONTINUOUS)
      {
        not_relevant_error_message("prev_constraint sub-block", pt_CONTINUOUS);
      }
      else
      {
        //lint -e{534}
        parse_prev_constraint(ptrs.prev_constraint);
      }
    }
    else
    {
      errors << priority(error) << "Both prev_constraint and prev_constraints"
             << " sub-blocks detected.  Will ignore prev_constraint sub-block."
             << endl;
    }
  }

  // PREV_CONSTRAINT
  if(process_block(ptrs.prev_constraints))
  {
    if(my_model.primary_trait_type == pt_CONTINUOUS)
    {
      not_relevant_error_message("prev_constraints sub-block", pt_CONTINUOUS);
    }
    else
    {
      parse_prev_constraints(ptrs.prev_constraints);
    }
  }
  
  // PREV_ESTIMATE
  if(process_block(ptrs.prev_estimate))
  {
    if(my_model.primary_trait_type == pt_CONTINUOUS)
    {
      not_relevant_error_message("prev_estimate sub-block", pt_CONTINUOUS);
    }
    else
    {
      parse_prev_estimate_sub_model(ptrs.prev_estimate);
    }
  }
}

void
parser::not_relevant_error_message(const string& sub_model_name, primary_type type)
{
  errors << priority(warning) << sub_model_name <<" not relevant for "
         << primary_type_2_string(type) << " trait.  "
         << "Ignoring ..." << endl;
}

void
parser::reset()
{
  my_model.reset(errors);
  transmission_parsed = false;
  fpmm_parsed = false;
  resid_parsed = false;
  transformation_parsed = false;
  file_name_root_parsed = false;

  ptrs = base_ptrs();
}

void
parser::init_genotype_map()
{
  //lint --e{534}

  my_genotype_map.clear();
  my_genotype_map.insert(pair<string, unsigned>("**", AA | AB | BB));
  my_genotype_map.insert(pair<string, unsigned>("A*", AA | AB ));
  my_genotype_map.insert(pair<string, unsigned>("B*", AB | BB));
  my_genotype_map.insert(pair<string, unsigned>("AA", AA));
  my_genotype_map.insert(pair<string, unsigned>("AB", AB));
  my_genotype_map.insert(pair<string, unsigned>("BB", BB));
}

// - Get 'value' and 'fixed' attributes.  Populate model_input.  Last 
//   argument determines whether a missing attribute is flagged as an
//   error or a warning.  def (default is a C++ keyword) is used for value 
//   if fixed is specified, but no value is given.  Value of QNAN for def
//   means that if fixed is 'true' a value must be supplied for 'val'.
//
// - value_open argument of 'true' means value or val attributes can be 
//   used for real number argument associated w. this parameter.  'false'
//   means only val can be used.
//
// - Return value of false means stop parsing and invalidate analysis.
//
bool
parser::get_attributes(model_input& mi, double def, const LSFBase* param,
                       const string& name_phrase, bool value_open) const
{
  // - 'val' attribute.
  //
  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
 
    // - Check for non-meaningful attributes.
    //
    for(iter = a_list->begin(); iter != a_list->end(); ++iter)
    {
      string  attr_name = AttrNameMgr.query(iter->first);
      if(! (toUpper(attr_name) == "VAL"         || 
            toUpper(attr_name) == "FIXED"       ||
            toUpper(attr_name) == "VALUE"       ||
            toUpper(attr_name) == "INTERACTION" ||
            toUpper(attr_name) == "LOWER_BOUND" ||
            toUpper(attr_name) == "UPPER_BOUND"    ))
      {
        errors << priority(error) << "Attribute '" << attr_name << "' of " << name_phrase 
               << " not understood.  Ignoring ..." << endl;
      }
    }
  
    if(value_open)
    {
      iter = a_list->find("VALUE");
      if(iter != a_list->end())
      {
        AttrVal a  = iter->second;
        if(a.has_value())
        {
          if(finite(a.Real()))
          {
            mi.value = a.Real();
          }
          else
          {
            errors << priority(error) << "Value of " << name_phrase
                   << " not understood.  Ignoring ..." << endl;
          }
        }
      }
    }
    
    iter = a_list->find("VAL");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          mi.value = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for 'val' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
    
    // - 'fixed' attribute.
    //
    iter = a_list->find("FIXED");
    if(iter != a_list->end())
    {
      string value; 
      value = toUpper(iter->second.String());
      if(value == "TRUE" || value == "YES" || value == "ON")
      {
        mi.fixed = true;
        
        // - User specified value for fixed, but not for val.  Use default
        //   for val.  If left as QNAN, sub-model will ignore entire model_input.
        // 
        if(SAGE::isnan(mi.value))
        {
        
          // - Fixed default is available.
          //
          if(! SAGE::isnan(def))
          {
            mi.value = def;
          }
          else
          {
            errors << priority(critical) << "No value given for " << name_phrase
                   << " with attribute, fixed, equal to 'true'.  Skipping analysis ..." << endl;
            return false;
          }
          
        }
      }
      else if(value == "FALSE" || value == "NO" || value == "OFF")
      {
        mi.fixed = false;    
        
        // - See above.
        //
        if(SAGE::isnan(mi.value))
        {
          mi.value = def;
        }
      }
      
      // - User tried unsuccessfully to give a value for fixed.  Don't make any
      //   assumptions about his intent.
      //
      else
      {
        errors << priority(error) << "Value of 'fixed' attribute of " 
               << name_phrase << " not understood.  Ignoring ..." << endl;
      }
    }
  }
  
  // - No attributes given for parameter.
  //
  else
  {
    errors << priority(warning) << "No 'val' or 'fixed' attribute(s) given for " << name_phrase 
           << "." << endl;
  }
  
  return true;
}

RPED::RefTraitInfo::trait_t
parser::primary_trait_type() const
{
  return my_mp->info().get_trait_type(my_model.primary_trait);
}

// - Parse 'type' sub-parameter (mean, var and suscept) and set 
//   model_inputs accordingly.
//
bool
parser::get_types(model_input types[], const LSFBase* param,
                  const string& keyword)
{
  assert(keyword == "mean"   ||
         keyword == "var"    ||
         keyword == "suscept"  ); 

  string sub_parameter = "";
  parse_string(param, sub_parameter);
  sub_parameter = toUpper(sub_parameter);
  
  std::map<string, unsigned>::const_iterator  map_iter;
  map_iter = my_genotype_map.find(sub_parameter);
  if(map_iter != my_genotype_map.end())
  {
    model_input  attributes;
    string name_phrase = keyword + " " + sub_parameter + " in type_" + keyword + " sub-block";
    
    if(! get_attributes(attributes, QNAN, param, name_phrase, false))
    {
      return false;
    }
    
    for(genotype_index i = index_AA; i < index_BB + 1; 
                i = static_cast<genotype_index>(i + 1))
    { 
      if(index_2_info(i) & map_iter->second)
      {
        types[i] = attributes;
      }
    }
  }
  else
  {
    errors << priority(error) << "Value '" << sub_parameter << "'  "
           << "for " << keyword <<" parameter in type_" << keyword << " sub-block not recognized.  "
           << "Ignoring ..." << endl;
  }

  
  return true;
}

// --- GENOTYPE SPECIFIC MEAN SUB-BLOCK ---
//
// - Parser for mean and suscept sub-models.
//
void  
parser::parse_mean_sub_model(const LSFBase* param, const string& keyword)
{
  assert(keyword == "mean" || keyword == "suscept");

  genotype_specific_mean_sub_model::sm_option  option;
  option = genotype_specific_mean_sub_model::one;
  
  // - model_input default values are: value == qNaN and fixed == false.
  // 
  model_input  means[NUM_OF_TYPES];     
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          string option_string = toUpper(a.String());
          if(option_string == "ONE")
          {
            option = genotype_specific_mean_sub_model::one;
          }
          else if(option_string == "TWO")
          {
            option = genotype_specific_mean_sub_model::two;
          }
          else if(option_string == "THREE")
          {
            option = genotype_specific_mean_sub_model::three;
          }
          else if(option_string == "TWO_DOM")
          {
            option = genotype_specific_mean_sub_model::two_dom;
          }
          else if(option_string == "TWO_REC")
          {
            option = genotype_specific_mean_sub_model::two_rec;
          }
          else if(option_string == "THREE_ADD")
          {
            option = genotype_specific_mean_sub_model::three_add;
          }
          else if(option_string == "THREE_DEC")
          {
            option = genotype_specific_mean_sub_model::three_dec;
          }
          else if(option_string == "THREE_INC")
          {
            option = genotype_specific_mean_sub_model::three_inc;
          }
          else
          { 
            errors << priority(error) << "Value '" << option_string << "' for"
                   << " option parameter of type_" << keyword
                   << " sub-block not recognized or not allowed. "
                   << " Ignoring ..." << endl;
          }
        }
      }
      else if(param_name == toUpper(keyword))
      {
        if(! get_types(means, *iter, keyword))
        {
          my_model.m_class = model_INVALID;
          return;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in type_" << keyword << " sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "Sub-block, type_" << keyword << " contains no parameters."
           << endl;
  }
  
  if(keyword == "mean")
  {
    if(! my_model.mean_sub_model.set(option, means[index_AA], 
                                     means[index_AB], means[index_BB],
                                     my_model.get_primary_trait_type()))
    {
      my_model.m_class = model_INVALID;
    }
  }
  else
  {
    if(! my_model.susc_sub_model.set(option, means[index_AA], 
                                     means[index_AB], means[index_BB],
                                     my_model.get_primary_trait_type()))
    {
      my_model.m_class = model_INVALID;
    }
  }
} 
  
// --- GENOTYPE SPECIFIC VARIANCE SUB-BLOCK ---
// 
void  
parser::parse_var_sub_model(const LSFBase* param)
{
  genotype_specific_variance_sub_model::sm_option  option 
        = genotype_specific_variance_sub_model::one;
  //bool  option_specified = false;
  
  // - model_input default values are: value == qNaN and fixed == false.
  // 
  model_input  vars[NUM_OF_TYPES];    
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          //option_specified = true;
        
          string option_string = toUpper(a.String());
          
          if(option_string == "ONE")
          {
            option = genotype_specific_variance_sub_model::one;
          }
          else if(option_string == "TWO")
          {
            option = genotype_specific_variance_sub_model::two;
            type_var_one = false; // due to JA
          }
          else if(option_string == "THREE")
          {
            option = genotype_specific_variance_sub_model::three;
            type_var_one = false; // due to JA
          }
          else if(option_string == "TWO_DOM")
          {
            option = genotype_specific_variance_sub_model::two_dom;
            type_var_one = false; // due to JA
          }
          else if(option_string == "TWO_REC")
          {
            option = genotype_specific_variance_sub_model::two_rec;
            type_var_one = false; // due to JA
          }
          else if(option_string == "THREE_ADD")
          {
            option = genotype_specific_variance_sub_model::three_add;
            type_var_one = false; // due to JA
          }
          else
// The _prop_mean options are being removed.
/*
          //lint -e{774}
          if(! PROP_MEAN_OPTIONS_AVAILABLE)
          {
            errors << priority(critical) << "var_prop_mean and std_prop_mean options "
                   << "in type_var sub-block not available in "
                   << "this release.  Skipping analysis ..." << endl;
            my_model.m_class = model_INVALID;
            return;
          }
          
          else if(option_string == "VAR_PROP_MEAN")
          {
            option = genotype_specific_variance_sub_model::var_prop_mean;
          }
          else if(option_string == "STD_PROP_MEAN")
          {
            option = genotype_specific_variance_sub_model::std_prop_mean;
          }
          else
*/
          { 
            //option_specified = false;
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of type_var sub-block not recognized.  "
                   << "Ignoring ..." << endl;
          }
        }
      }
      else if(param_name == "VAR")
      {
        if(! get_types(vars, *iter, "var"))
        {
          my_model.m_class = model_INVALID;
          return;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in type_var sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "type_var sub-block contains no parameters."
           << endl;
  }
  
  if(my_model.primary_trait_type == pt_ONSET){
      if (! type_var_one) {
          errors << priority(critical) << " Cannot specify more than one variance with age of " << 
          " onset model ... Ignoring Analysis " << endl;
       my_model.m_class = model_INVALID;
   } 
 }  

  if(! my_model.var_sub_model.set(option, vars[index_AA], vars[index_AB], vars[index_BB],
                                  my_model.mean_missing, my_model.trans_missing,
                                  my_model.get_primary_trait_type()))
  {
    my_model.m_class = model_INVALID;
  }
} 


// --- GENOTYPE FREQUENCY SUB-BLOCK ---
//
void
parser::get_probs_fixed(bool& probs_fixed, bool& probs_fixed_specified, 
                        const LSFBase* param)
{
  string  value = "";
  parse_string(param, value);
  
  if(value.length())
  {
    value = toUpper(value);
    if(value == "TRUE" || value == "YES" || value == "ON")
    {
      probs_fixed = true;
      probs_fixed_specified = true;
    }
    else if(value == "FALSE" || value == "NO" || value == "OFF")
    {
      probs_fixed = false;
      probs_fixed_specified = true;    
    }
    else
    {
      errors << priority(error) << "Value of geno_freq sub-block parameter, " 
             << "probs_fixed, not understood.  Ignoring ..." << endl;
    }
  }
}

void
parser::get_probs(double probs[], const LSFBase* param)
{
  string  name_phrase = "prob in geno_freq sub-block";

  string prob_name = "";
  parse_string(param, prob_name);
  prob_name = toUpper(prob_name);
  
  if(prob_name == "AA" || prob_name == "AB" || prob_name == "BB")
  {
    // - Get 'val' attribute.
    //
    AttrList* a_list = param->attrs();
    AttrList::const_iterator iter;
    if(a_list)
    {
      // - Check for non-meaningful attributes.
      //
      for(iter = a_list->begin(); iter != a_list->end(); ++iter)
      {
        string  attr_name = AttrNameMgr.query(iter->first);
        if(! (toUpper(attr_name) == "VAL"   ||
              toUpper(attr_name) == "VALUE"    ))
        {
          if(toUpper(attr_name) == "FIXED")
          {
            errors << priority(error) << "Attribute 'fixed' not used in parameter, "
                   << name_phrase << ".  Ignoring ...  Use 'probs_fixed'." << endl;
          }
          else
          
          {
            errors << priority(error) << "Attribute '" << attr_name << "' of " << name_phrase
                   << " not understood.  Ignoring ..." << endl;
          }
        }
      }
    
      iter = a_list->find("VAL");
      if(iter != a_list->end())
      {
        AttrVal a  = iter->second;
        if(a.has_value())
        {
          if(finite(a.Real()))
          {
            double  prob_value = a.Real();
            if(prob_name == "AA")
            {
              probs[index_AA] = prob_value;
            }
            else if(prob_name == "AB")
            {
              probs[index_AB] = prob_value;
            }
            else if(prob_name == "BB")
            {
              probs[index_BB] = prob_value;
            }
          }
          else
          {
            errors << priority(error) << "Value for 'val' attribute of " << name_phrase
                   << " not understood.  Ignoring ..." << endl;
          }
        }
      }
    }
    else
    {
      errors << priority(warning) << "No attribute(s) given for " << name_phrase
             << "." << endl;
    }
  }
  else 
  {
    errors << priority(error) << "Value given for " << name_phrase << "not understood.  "
           << "Ignoring ..." << endl;
  }
}

void
parser::get_freq_A(double& freq_A, const LSFBase* param)
{
  string  name_phrase = "freq_A in geno_freq sub-block";

  // - Get 'val' attribute.
  //
  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    // - Check for non-meaningful attributes.
    //
    for(iter = a_list->begin(); iter != a_list->end(); ++iter)
    {
      string  attr_name = AttrNameMgr.query(iter->first);
      if(! (toUpper(attr_name) == "VAL"   || 
            toUpper(attr_name) == "VALUE"    ))
      { 
        if(toUpper(attr_name) == "FIXED")
        {
          errors << priority(error) << "Attribute 'fixed' not used in parameter, "
                 << name_phrase << ".  Use 'probs_fixed'." << endl;
        }
        else
        {
          errors << priority(error) << "Attribute '" << attr_name << "' of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  
    iter = a_list->find("VALUE");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          freq_A = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for 'value' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
 
    iter = a_list->find("VAL");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          freq_A = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for 'val' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  }
  else
  {
    errors << priority(warning) << "No attribute(s) given for " << name_phrase
           << "." << endl;
  }
}

void  
parser::parse_freq_sub_model(const LSFBase* param)
{
  genotype_frequency_sub_model::sm_option  option 
        = genotype_frequency_sub_model::hwe;
  //bool         option_specified = false;
  
  bool         probs_fixed = PROBS_FIXED_DEFAULT;
  bool         probs_fixed_specified = false;
  double       probs[] = { QNAN, QNAN, QNAN };
  double       freq_A = QNAN;
  model_input  corr(QNAN, CORR_DEFAULT_FIXED); 
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          string option_string = toUpper(a.String());
          if(option_string == "HWE")
          {
            option = genotype_frequency_sub_model::hwe;
            //option_specified = true;
          }
          else if(option_string == "NHWE")
          {
            option = genotype_frequency_sub_model::nhwe;
            //option_specified = true;
          }
          else
          {
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of geno_freq sub-block not recognized.  "
                   << "Using ..." << endl;
          }
        }
      }
      else if(param_name == "PROBS_FIXED")
      {
        get_probs_fixed(probs_fixed, probs_fixed_specified, *iter);  
      }
      else if(param_name == "PROB")
      {
        get_probs(probs, *iter);
      }
      else if(param_name == "FREQ_A")
      {
        get_freq_A(freq_A, *iter);        
      }
      else if(param_name == "CORR")
      {
      
        string  name_phrase = "corr in geno_freq sub-block";
        if(! get_attributes(corr, QNAN, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in geno_freq sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "geno_freq sub-block contains no parameters."
           << endl;
  }

  genotype_specific_mean_susc_sub_model& type_sub_model = 
                                       my_model.type_dependent_sub_model();

  if(! my_model.freq_sub_model.set(option, freq_A, probs[index_AA], probs[index_AB], 
                                   probs[index_BB], corr, type_sub_model, 
                                   my_model.get_type_missing(), my_model.trans_missing, probs_fixed))
  {
    my_model.m_class = model_INVALID;
  } 
}
  

// --- RESIDUAL CORRELATION SUB-BLOCK ---
//    
void  
parser::parse_resid_sub_model(const LSFBase* param)
{
  if(param->List())
  {
    residual_correlation_sub_model::sm_option  option
          = residual_correlation_sub_model::equal_po_ss;
    
    model_input  fm_mi(QNAN, true);
    model_input  mo_mi(QNAN, RESID_DEFAULT_FIXED);
    model_input  fo_mi(QNAN, RESID_DEFAULT_FIXED);
    model_input  ss_mi(QNAN, RESID_DEFAULT_FIXED);
  
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          string option_string = toUpper(a.String());
          if(option_string == "EQUAL_PO_SS")
          {
            option = residual_correlation_sub_model::equal_po_ss;
          }
          else if(option_string == "EQUAL_PO")
          {
            option = residual_correlation_sub_model::equal_po;
          }
          else if(option_string == "ARB")
          {
            option = residual_correlation_sub_model::arb;
          }
          else
          {
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of resid sub-block not recognized.  "
                   << "Using ..." << endl;
          }
        }
      }
      else if(param_name == "FM")
      {
        string  name_phrase = "fm in resid sub-block";
        
        // - A default value IS available for this parameter.
        //
        if(! get_attributes(fm_mi, RESID_SP_DEFAULT_VALUE, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
      }
      else if(param_name == "MO")
      {
        string  name_phrase = "mo in resid sub-block";
        
        if(! get_attributes(mo_mi, QNAN, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
      }
      else if(param_name == "FO")
      {
        string  name_phrase = "fo in resid sub-block";
        
        if(! get_attributes(fo_mi, QNAN, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
      }
      else if(param_name == "SS")
      {
        string  name_phrase = "ss in resid sub-block";
        
        if(! get_attributes(ss_mi, QNAN, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in resid sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }

    if(! my_model.resid_sub_model.set(option, fm_mi, mo_mi, fo_mi, ss_mi,
                                      my_model.m_class))
    {
      my_model.m_class = model_INVALID;
    }  

    my_model.resid_sub_model.calculate_alpha_and_delta();
  }
  else
  {
    my_model.resid_sub_model.set_as_default(my_model.m_class);

    my_model.resid_sub_model.calculate_alpha_and_delta();

    errors << priority(warning) << "resid sub-block contains no parameters."
           << endl;
  }
}
  

// --- TRANSMISSION SUB-BLOCK ---
//    
// - Parse 'transmission' sub-parameter and set model_inputs accordingly.
//
bool
parser::get_taus(model_input taus[], const LSFBase* param)
{
  string sub_parameter = "";
  parse_string(param, sub_parameter);
  sub_parameter = toUpper(sub_parameter);
  
  std::map<string, unsigned>::const_iterator  map_iter;
  map_iter = my_genotype_map.find(sub_parameter);
  if(map_iter != my_genotype_map.end())
  {
    model_input  attributes;
    string name_phrase = "tau " + sub_parameter + " in transmission sub-block";
    
    if(! get_attributes(attributes, QNAN, param, name_phrase, false))
    {
      return false;
    }
    
    for(genotype_index i = index_AA; i < index_BB + 1; 
                i = static_cast<genotype_index>(i + 1))
    { 
      if(index_2_info(i) & map_iter->second)
      {
        taus[i] = attributes;
      }
    }
  }
  else
  {
    errors << priority(error) << "Value '" << sub_parameter << "'  "
           << "for tau parameter in transmission sub-block not recognized.  "
           << "Ignoring ..." << endl;
  }
  
  return true;
}

void
parser::get_no_bounds(bool& no_bounds, const LSFBase* param)
{
  string  value = "";
  parse_string(param, value);
  
  if(value.length())
  {
    value = toUpper(value);
    if(value == "TRUE" || value == "YES" || value == "ON")
    {
      no_bounds = true;
    }
    else if(value == "FALSE" || value == "NO" || value == "OFF")
    {
      no_bounds = false;
    }
    else
    {
      errors << priority(error) << "Value of transmission sub-block parameter, " 
             << "no_bounds, not understood.  Ignoring ..." << endl;
    }
  }
}

void  
parser::parse_transmission_sub_model(const LSFBase* param)
{
  transmission_sub_model::sm_option  option 
        = transmission_sub_model::homog_no_trans;
  //bool  option_specified = false;
  
  // - model_input default values are: value == qNaN and fixed == false.
  // 
  model_input  taus[3];

  bool  no_bounds = false;     
  bool  no_bounds_specified = false;
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
     string  param_name = toUpper((*iter)->name());
      if(param_name == "OPTION")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          //option_specified = true;
        
          string option_string = toUpper(a.String());
          
          if(option_string == "NO_TRANS")
          {
            option = transmission_sub_model::no_trans;
          }
          else if(option_string == "HOMOG_NO_TRANS")
          {
            option = transmission_sub_model::homog_no_trans;
          }
          else if(option_string == "HOMOG_MENDELIAN")
          {
            option = transmission_sub_model::homog_mendelian;
          }
          else if(option_string == "HOMOG_GENERAL")
          {
            option = transmission_sub_model::homog_general;
          }
          else if(option_string == "GENERAL")
          {
            option = transmission_sub_model::general;
          }
          else if(option_string == "TAU_AB_FREE")
          {
            option = transmission_sub_model::tau_ab_free;
          }
          else if(option_string == "MITOCHONDRIAL")
          {
            option = transmission_sub_model::mitochondrial;
          }
          else
          { 
            //option_specified = false;
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of transmission sub-block not recognized.  "
                   << "Using ..." << endl;
          }
        }
      }
      else if(param_name == "TAU")
      {
        if(! get_taus(taus, *iter))
        {
          my_model.m_class = model_INVALID;
          return;
        }
      }
      else if(param_name == "NO_BOUNDS")
      {
        get_no_bounds(no_bounds, *iter);
        no_bounds_specified = true;
      }  
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in transmission sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  
  else
  {
    errors << priority(warning) << "Transmission sub-block contains no parameters."
           << endl;
  }
  
  // - no_bounds not relevant for all options.
  //
  if(no_bounds_specified && (option == transmission_sub_model::homog_no_trans  ||
                             option == transmission_sub_model::homog_mendelian ||
                             option == transmission_sub_model::no_trans        ||
                             option == transmission_sub_model::mitochondrial      ))
  {
    errors << priority(error) << "Specifed transmission sub-block parameter, no_bounds, "
           << "not relevant for option, " << transmission_sub_model::option_2_parameter(option) 
           << ".  Ignoring ..." << endl;
  }

  if(! my_model.transm_sub_model.set(option, taus[index_AA], 
                                   taus[index_AB], taus[index_BB], no_bounds,
                                   my_model.get_type_missing()))
  {
    my_model.m_class = model_INVALID;
  }
}  

// --- TRANSFORMATION SUB-BLOCK ---
//

void  
parser::parse_transformation_sub_model(const LSFBase* param)
{
  MAXFUN::TransformationSubmodel::sm_option  option;

  if(my_model.primary_trait_type == pt_BINARY)
      option = MAXFUN::TransformationSubmodel::no_trans;
  else
      option = MAXFUN::TransformationSubmodel::box_cox;       

  if(!MAXFUN::parseTransformationSubmodel(my_model.transf_sub_model, param, option, errors))
  {
    my_model.m_class = model_INVALID;
  }
}


// --- COVARIATE ---
//
void  
parser::get_interaction(bool& interaction, const LSFBase* param, 
                        covariate_type type, const string& name_phrase) // const
{
  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    iter = a_list->find("INTERACTION");
    if(iter != a_list->end())
    {
      // - Interactions not allowed for composite trait
      //   sub-model.
      //
      if(type == comp)
      {
        errors << priority(error) << "Interaction specified for "
               << name_phrase << ".  Ignoring ..." << endl;
        return;
      }
    
      string value; 
      value = toUpper(iter->second.String());
      if(value == "TRUE" || value == "YES" || value == "ON")
      {
        interaction= true;
      }
      else if(value == "FALSE" || value == "NO" || value == "OFF")
      {
        interaction = false;    
      }
      else
      {
        errors << priority(error) << "Value of 'interaction' attribute of " 
               << name_phrase << " not understood.  Ignoring ..." << endl;
      }
    }
  }
}

bool
parser::get_covariate(model_input& coeff, bool& interaction, 
                      covariate_type type, const LSFBase* param) //  const
{
  string  name_phrase = "covariate in " + CovariateSubmodel::get_cov_type_subblock_name(type) + " sub-block";
  
  if(! get_attributes(coeff, QNAN, param, name_phrase, false))
  {
    return false;
  }
  
  get_interaction(interaction, param, type, name_phrase);
  
  // - Model has been invalidated if ! INTERACTIONS_AVAILABLE.
  //
  if(!model_valid())
  {
    return false;
  }
  
  if(type == comp)
  {
    interaction = false;
  }
  
  return true;
}

void  
parser::reset_covariate_sub_model(covariate_type type)
{
  switch(type)
  {
    case mean:
      my_model.mean_cov_sub_model = 
          MeanCovariateSubmodel(&(my_model.mean_sub_model), errors);
      break;
    case VAR:
      my_model.var_cov_sub_model =
          VarianceCovariateSubmodel(&(my_model.mean_sub_model), errors);
      break;
    case susc:
      my_model.susc_cov_sub_model =
          SusceptibilityCovariateSubmodel(&(my_model.susc_sub_model), errors);
      break;
    case comp:
      my_model.comp_trait_sub_model =
          CompositeTraitSubmodel(&(my_model.mean_sub_model), errors);
      break;
    default:
      assert(false);
  }
}

void  
parser::parse_cov_sub_model(const LSFBase* param, covariate_type type)
{
  reset_covariate_sub_model(type);
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "COVARIATE") 
      {
        model_input  coeff = model_input(QNAN, COV_DEFAULT_FIXED);
        bool         interaction = COV_DEFAULT_INTERACTION; 
        string       trait_name = "";
        
        parse_string(*iter, trait_name);
        
        if(! get_covariate(coeff, interaction, type, *iter))
        {
          my_model.m_class = model_INVALID;
          return;
        }
        
        bool  set_valid = false;
        switch(type)
        {
          case mean:
            set_valid = my_model.mean_cov_sub_model.add_covariate
                (my_mp,
                 trait_name,
                 coeff.value,
                 interaction,
                 coeff.fixed,
                 my_model.mean_missing);
                 
            break;
          case VAR:
            set_valid = my_model.var_cov_sub_model.add_covariate
                (my_mp,
                 trait_name,
                 coeff.value,
                 interaction,
                 coeff.fixed,
                 my_model.mean_missing);
            break;
          case susc:
            set_valid = my_model.susc_cov_sub_model.add_covariate
                (my_mp,
                 trait_name,
                 coeff.value,
                 interaction,
                 coeff.fixed,
                 my_model.mean_missing);
            break;
          case comp:
            set_valid = my_model.comp_trait_sub_model.add_covariate
                (my_mp,
                 trait_name,
                 coeff.value,
                 coeff.fixed,
                 my_model.mean_missing);
            break;
          default:
            assert(false);
        }
        
        if(! set_valid)
        {
          my_model.m_class = model_INVALID;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in "
               << CovariateSubmodel::get_cov_type_subblock_name(type) << " sub-block not recognized.  "
               << "Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << CovariateSubmodel::get_cov_type_subblock_name(type) << " contains no parameters." 
           << endl;
  }
}


// --- FINITE POLYGENIC MIXED MODEL SUB-BLOCK ---
//    
void  
parser::parse_fpmm_sub_model(const LSFBase* param)
{
  model_input  var(QNAN, FPMM_DEFAULT_FIXED);
  double       freq = QNAN;
  size_t       loci = (size_t)(-1);

  bool         found_onset = false;

  bool         found_freq  = false; // due to JA
  bool         found_var   = false; // due to JA
  
  if(! FPMM_AVAILABLE)
  {
    errors << priority(critical) << "FPMM sub-block not available in "
           << "this release.  Skipping analysis ..." << endl;
    my_model.m_class = model_INVALID;
    return;
  }
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "VAR")
      {
        found_var = true; // due to JA
        string  name_phrase = "var in fpmm sub-block";
        if(! get_attributes(var, QNAN, *iter, name_phrase, true))
        {
          my_model.m_class = model_INVALID;
          return; 
        }
        
      }
      else if(param_name == "FREQ")
      {
        AttrVal  a = attr_value(*iter, 0);
        if(a.has_value())
        {
          if(finite(a.Real()))
          {
            freq = a.Real();
            found_freq = true; // due to JA
          }
          else
          {
            errors << priority(error) << "Value for 'freq' parameter in fpmm sub-block"
                   << " not understood.  Ignoring ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "Value for 'freq' parameter in fpmm sub-block"
                 << " missing.  Ignoring ..." << endl;
        }
      }
      else if(param_name == "LOCI")
      {
        AttrVal  a = attr_value(*iter, 0);
        if(a.has_value())
        {
          if(finite(a.Real()) && a.Int() >= 0)
          {
                  loci = (uint) a.Int();
                  if (loci == 0)  // due to JA
                  { zero_poly_loci = true;} 
                  else { zero_poly_loci = false; }
          }
          else
          {
            errors << priority(error) << "Value for 'loci' parameter in fpmm sub-block"
                   << " not understood.  Ignoring ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "Value for 'loci' parameter in fpmm sub-block"
                 << " missing.  Ignoring ..." << endl;
        }
      }
      else if(param_name == "ONSET")
      {
        // Verify that onset is specified in the primary trait

        if(my_model.primary_trait_type == pt_ONSET)
        {
          found_onset = true;

          parse_onset_sub_model(*iter);
        }
        else
        {
          errors << priority(critical) << "Onset sub-block detected, but "
                 << "primary trait is not listed as age_onset.  Skipping "
                 << "analysis ..." << endl;
          my_model.m_class = model_INVALID;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in fpmm sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
     if (zero_poly_loci && found_var) 
      {
      errors << priority(critical) << " Cannot specify polygenic locus variance with " << 
                        "zero polygenic loci ... Skipping analysis" << endl;
                        my_model.m_class = model_INVALID;
      } // due to JA
     if (zero_poly_loci && found_freq) 
      {
      errors << priority(critical) << " Cannot specify polygenic locus frequency with " << 
                        "zero polygenic loci ... Skipping analysis " << endl;
                        my_model.m_class = model_INVALID;
      } // due to JA
  }
  else
  {
    errors << priority(warning) << "fpmm sub-block contains no parameters."
           << endl;
  }

  if(my_model.primary_trait_type == pt_ONSET && !found_onset)
  {
    errors << priority(critical) << "No onset sub-block found within fpmm sub-block "
           << "when onset model specified. Analysis invalid. Skipping analysis." << endl;
    my_model.m_class = model_INVALID;
  }

  // - Check model validity since setting of onset is nested w/i this sub-model
  //   and may have invalidated the model.
  //
  
  if(model_valid())
  {
    if(zero_poly_loci) loci =0; // just a check
    if(! my_model.fpmm_sub_model.set(var, freq, loci))
    {
      my_model.m_class = model_INVALID;
    }
  }  
}

// --- ONSET SUB-BLOCK ---
//    
void  
parser::parse_onset_sub_model(const LSFBase* param)
{
  onset_sub_model::type_option   t_option = onset_sub_model::t_A;
  onset_sub_model::multi_option  m_option = onset_sub_model::m_N;

  string  onset;
  string  exam;

  
  bool type_dependent_found = false; // added by JA to tighten up parser
  
//  bool age_onset_found = false;
//  bool age_exam_found = false;

  if(! ONSET_AVAILABLE)
  {
    errors << priority(critical) << "Onset sub-block not available in "
           << "this release.  Skipping analysis ..." << endl;
    my_model.m_class = model_INVALID;
    return;
  }

  assert(param != 0);
  assert(param->List() != 0);
  
  LSFList::const_iterator  iter;
  for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
  {
    string  param_name = toUpper((*iter)->name());
    if(param_name == "TYPE_DEPENDENT")
    {
      AttrVal a = attr_value(*iter, 0);
      if(a.has_value())
      {
        string t_option_string = toUpper(a.String());
        if(t_option_string == "A")
        {
          t_option = onset_sub_model::t_A;
          type_dependent_found = true; // new addition by JA
        }
        else if(t_option_string == "S")
        {
          t_option = onset_sub_model::t_S;
          type_dependent_found = true; // new addition by JA
        }
        else
        {
          errors << priority(error) << "Value '" << t_option_string << "' for "
                 << "type_dependent parameter of onset sub-block not recognized.  "
                 << "Ignoring ... " << endl;
        }
      }        
    }
    else if(param_name == "MULTI_DEPENDENT")
    {
      if (zero_poly_loci)  // addition by JA
      {
        errors << priority(critical) << " Cannot specify if polygenic loci affect age of onset or " <<
         "susceptibility with zero polygenic loci ... Ignoring" << endl;
          my_model.m_class = model_INVALID;
          return;
      } 
      AttrVal a = attr_value(*iter, 0);
      if(a.has_value())
      {
        string m_option_string = toUpper(a.String());
        if(m_option_string == "N")
        {
          m_option = onset_sub_model::m_N;
        }
        else if(m_option_string == "A")
        {
          m_option = onset_sub_model::m_A;
        }
        else if(m_option_string == "S")
        {
          m_option = onset_sub_model::m_S;
        }
        else
        {
          errors << priority(error) << "Value '" << m_option_string << "' for "
                 << "multi_dependent parameter of onset sub-block not recognized.  "
                 << "Ignoring ... " << endl;
        }
      }        
    }
    else if(param_name == "AGE_ONSET")
    {
      string  value;
      parse_string(*iter, value);

      // Don't have to check trait here because onset sub model checks it

      if(value.size())
      {
        onset = value;
      }
      else
      {
        errors << priority(error) << "No value given for 'age_onset' parameter "
               << "in onset sub-block." << endl;
      }
//      age_onset_found = true;
    }
    else if(param_name == "AGE_EXAM")
    {
      string  value;
      parse_string(*iter, value);

      // Don't have to check trait here because onset sub model checks it

      if(value.size())
      {
        exam = value;
      }
      else
      {
        errors << priority(error) << "No value given for 'age_exam' parameter "
               << "in onset sub-block." << endl;
      }
//      age_exam_found = true;
    }
    else
    {
      errors << priority(error) << "Parameter '" << param_name << "' in onset sub-block "
             << "not recognized.  Ignoring ..." << endl;
    }
  }
/*
  if (!(age_onset_found)) {
      errors << priority(error) << " Must specify age_onset in onset sub-block when trait is declared binary with variable " <<
      " age of onset  ... Ignoring " << endl;
      my_model.m_class = model_INVALID;
      return;
      }

  if (!(age_exam_found)) {
      errors << priority(error) << " Must specify age of exam in onset sub-block when trait is declared binary with variable " <<
      " age of onset  ... Ignoring " << endl;
      my_model.m_class = model_INVALID;
      return;
      }
*/
  if (zero_poly_loci) {
     if (m_option != onset_sub_model::m_N) {
      errors << priority(error) << " Cannot allow polygenic loci to affect either age of onset or susceptibility " <<
      "with no polygenic loci ... Ignoring" << endl; 
       my_model.m_class = model_INVALID;
       return; }
  }   

  // Set the model.  This returns an error state boolean

  if(! my_model.ons_sub_model.set(t_option, m_option, my_mp,
                                  my_model.get_primary_trait(), onset, exam))
  {
    my_model.m_class = model_INVALID;

    return;
  }

  // Checks for other potential errors

  // Check that age onset and age exam are not in the covariates

  if(my_model.mean_cov_sub_model.has_covariate(onset) ||
     my_model.mean_cov_sub_model.has_covariate(exam)  ||
     my_model.var_cov_sub_model .has_covariate(onset) ||
     my_model.var_cov_sub_model .has_covariate(exam)  ||
     my_model.susc_cov_sub_model.has_covariate(onset) ||
     my_model.susc_cov_sub_model.has_covariate(exam) )
  {
    errors << priority(critical) << "Age of Onset and Age at Exam cannot be "
           << "in covariate lists.  Skipping analysis ..." << endl;
    my_model.m_class = model_INVALID;

    return;
  }

  // Check that the !type_dependent parameter (mean/susc) is set to a single
  // mean.  Also, determine if a commingling analysis is required.

  if(t_option == onset_sub_model::t_A)
  {
    if(my_model.susc_sub_model.option() != genotype_specific_sub_model::one)
    {
      errors << priority(critical) << "Susceptibilities have more than one "
             << "type, but age is type dependent in age of onset model. "
             << "Analysis invalid.  Skipping ... " << endl;
      my_model.m_class = model_INVALID;

      return;

    }
  }
  else // suscept is type dependent
  {
    if(my_model.mean_sub_model.option() != genotype_specific_sub_model::one)
    {
      errors << priority(critical) << "Means have more than one "
             << "type, but susceptibility is type dependent in age of onset model. "
             << "Analysis invalid.  Skipping ... " << endl;
      my_model.m_class = model_INVALID;

      return;

    }
  }
}

// --- ASCERTAINMENT SUB-BLOCK ---
//    
void  
parser::parse_ascertainment_sub_model(const LSFBase* param)
{
  ascertainment_sub_model::s_sm_option  s_option = ascertainment_sub_model::none;
  ascertainment_sub_model::v_sm_option  v_option = ascertainment_sub_model::not_specified;
  string          psf_indic;                  
  string          thresh_indicator;
  vector<double>  psf_includes; 
  double          thresh = QNAN;
  double          thresh_high = QNAN;
  double          thresh_low = QNAN;
  double          indic_thresh = QNAN;

  bool  s_option_specified = false;
  
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "COND_SET")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value())
        {
          string option_string = toUpper(a.String());
          if(option_string == "NONE")
          {
            s_option = ascertainment_sub_model::none;
            s_option_specified = true;
          }
          else if(option_string == "FOUNDERS")
          {
            s_option = ascertainment_sub_model::founders;
            s_option_specified = true;
          }
          else if(option_string == "PSF")
          {
            s_option = ascertainment_sub_model::psf;
            s_option_specified = true;
          }
          else if(option_string == "FOUNDERS_PLUS_PSF")
          {
            s_option = ascertainment_sub_model::founders_plus_psf;
            s_option_specified = true;
          }
          else
          {
            errors << priority(error) << "Value '" << option_string << "' for "
                   << " parameter, cond_set, of ascertainment sub-block not recognized.  "
                   << "Ignoring ... " << endl;
          }
        }        
      }
      else if(param_name == "PSF_INDIC")
      {
        string  value;
        parse_string(*iter, value);
        if(value.size())
        {
          psf_indic = value;
        }
        else
        {
          errors << priority(error) << "No value given for 'psf_indic' parameter "
                 << "in ascertainment sub-block." << endl;
        }
      }
      else if(param_name == "PSF_INDIC_INCLUDE")
      {
        AttrVal  a = attr_value(*iter, 0);
        if(a.has_value())
        {
          if(finite(a.Real()))
          {
            psf_includes.push_back(a.Real());
          }
          else
          {
            errors << priority(error) << "Value for parameter, psf_indic_include, in ascertainment sub-block"
                   << " not understood.  Ignoring ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "Value for parameter, psf_indic_include, in ascertainment sub-block"
                 << " missing.  Ignoring ..." << endl;
        }
      }
      else if(param_name == "COND_VAL")
      {
        if(my_model.primary_trait_type == pt_CONTINUOUS)
        {
          AttrVal a = attr_value(*iter, 0);
          if(a.has_value())
          {
            string option_string = toUpper(a.String());
            if(option_string == "ACTUAL")
            {
              v_option = ascertainment_sub_model::actual;
            }
            else if(option_string == "GTE_THRESH")
            {
              v_option = ascertainment_sub_model::gte_thresh;
            }
            else if(option_string == "LTE_THRESH")
            {
              v_option = ascertainment_sub_model::lte_thresh;
            }
            else if(option_string == "THRESH_INDIC")
            {
              v_option = ascertainment_sub_model::thresh_indic;
            }
            else
            {
              errors << priority(error) << "Value '" << option_string << "' for "
                     << "cond_val parameter of ascertainment sub-block not recognized.  "
                     << "Ignoring ..." << endl;
            }
          }
          
          // - Get cond_val attributes.
          //
          if(v_option != ascertainment_sub_model::actual)
          {
            LSFBase::AList*  a_list = (*iter)->attrs();
            AttrList::const_iterator  l_iter = a_list->find("THRESH");
            if(l_iter != a_list->end())
            {
              AttrVal a2  = l_iter->second;
              if(a2.has_value())
              {
                if(finite(a2.Real()))
                {
                  thresh = a2.Real();
                }
                else
                {
                  errors << priority(error) << "Value for attribute, thresh, of " 
                         << "parameter, cond_val, in ascertainment sub-block "
                         << " not understood.  Ignoring ..." << endl;
                }
              }
            }
            
            l_iter = a_list->find("THRESH_INDIC_HIGH");
            if(l_iter != a_list->end())
            {
              AttrVal a2  = l_iter->second;
              if(a2.has_value())
              {
                if(finite(a2.Real()))
                {
                  thresh_high = a2.Real();
                }
                else
                {
                  errors << priority(error) << "Value for attribute, thresh_indic_high, of " 
                         << "parameter, cond_val, in ascertainment sub-block "
                         << " not understood.  Ignoring ..." << endl;
                }
              }
            }
            
            l_iter = a_list->find("THRESH_INDIC_LOW");
            if(l_iter != a_list->end())
            {
              AttrVal a2  = l_iter->second;
              if(a2.has_value())
              {
                if(finite(a2.Real()))
                {
                  thresh_low = a2.Real();
                }
                else
                {
                  errors << priority(error) << "Value for attribute, thresh_indic_low, of " 
                         << "parameter, cond_val, in ascertainment sub-block "
                         << " not understood.  Ignoring ..." << endl;
                }
              }
            }
          }
        }
        else  // primary trait is continuous or age of onset.
        {
          // We only set the v_option if it hasn't already been set

          if(v_option == ascertainment_sub_model::not_specified)
            v_option = ascertainment_sub_model::actual;

          AttrVal a = attr_value(*iter, 0);
          if(a.has_value())
          {

            string option_string = toUpper(a.String());
            if(option_string != "ACTUAL")
            {
              errors << priority(error) << "Parameter, cond_val, in ascertainment sub-block "
                   << "not relevant for binary or age of onset traits.  Ignoring ..." << endl;
            }
          }
        }
      }
      
      else if(param_name == "THRESH_INDIC")
      {
        string  value;
        parse_string(*iter, value);

        if(value.size())
        {
          thresh_indicator = value;
          
          // - Get thresh attribute.
          //
          LSFBase::AList*  a_list = (*iter)->attrs();
          AttrList::const_iterator  l_iter = a_list->find("THRESH");
          if(l_iter != a_list->end())
          {
            AttrVal a  = l_iter->second;
            if(a.has_value())
            {
              if(finite(a.Real()))
              {
                indic_thresh = a.Real();
              }
              else
              {
                errors << priority(error) << "Value for attribute, thresh, of " 
                       << "parameter, thresh_indic, in ascertainment sub-block "
                       << " not understood.  Ignoring ..." << endl;
              }
            }
          }
        }
      }
      else if(param_name == "ONSET_OPTION")
      {
        // If we're not doing age of onset, this should not be parsed
        if(my_model.primary_trait_type == pt_ONSET)
        {
          AttrVal a = attr_value(*iter, 0);
          if(a.has_value())
          {
            string option_string = toUpper(a.String());
            if(option_string == "ACTUAL")
            {
              v_option = ascertainment_sub_model::actual;
            }
            else if(option_string == "BY_ONSET")
            {
              v_option = ascertainment_sub_model::onset;
            }
            else  // value not recognized
            {
              errors << priority(error) << "Value '" << option_string << "' for "
                     << "onset_option parameter of ascertainment sub-block not recognized.  "
                     << "Ignoring ..." << endl;
            }
          }
          else // No value specified
          {
            errors << priority(error) << "No onset_option "
                   << "value specified in ascertainment sub-block.  "
                   << "Ignoring ..." << endl;
          }
        }
        else // Not an onset option
        {
          errors << priority(warning) << "Parameter, onset_option, in "
                 << "ascertainment sub-block only required for age of onset "
                 << "traits.  Ignoring ..." << endl;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in ascertainment sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "ascertainment sub-block contains no parameters."
           << endl;
  }
  
  // - cond_set option default is 'psf' if psf indicator trait is specified.
  //
  if(! s_option_specified && psf_indic.size())
  {
    s_option = ascertainment_sub_model::psf;
  }
  
  // - cond_val option default is 'actual', but it is not relevant when cond_set 
  //   option is 'none'.
  //
  if(s_option != ascertainment_sub_model::none && 
     v_option == ascertainment_sub_model::not_specified) 
  {
    v_option = ascertainment_sub_model::actual;
  }
  
  if(! my_model.ascer_sub_model.set(s_option, v_option, my_mp, psf_indic, psf_includes, 
                                    thresh_indicator, thresh, thresh_high, 
                                    thresh_low, indic_thresh))
  {
    my_model.m_class = model_INVALID;
  }
}

// --- PREVALENCE CONSTRAINT SUB-BLOCK ---
//    
size_t
parser::get_int_attribute(const LSFBase* param, const std::string& attribute,
                                                const std::string& parameter,
                                                const std::string& sub_block )
{
  size_t  return_value = (size_t)(-1);
  string  name_phrase = parameter + " in " + sub_block + " sub-block";

  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    // - Check for non-meaningful attributes.
    //
    for(iter = a_list->begin(); iter != a_list->end(); ++iter)
    {
      string  attr_name = AttrNameMgr.query(iter->first);
      if(! (toUpper(attr_name) == "R"     ||
            toUpper(attr_name) == "N"     ||
            toUpper(attr_name) == "VAL"   ||
            toUpper(attr_name) == "VALUE"    ))
      {
        if(toUpper(attr_name) == "FIXED")
        {
          errors << priority(error) << "Attribute, fixed, not used in parameter, "
                 << name_phrase << ".  Ignoring ..." << endl;
        }
        else
        {
          errors << priority(error) << "Attribute, " << attr_name << ", of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  
    iter = a_list->find(toUpper(attribute));
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        int  value = a.Int();
        if(finite(value))
        {
          if(value >= 0)
          {
            return_value = (uint) value;
          }
        }
        else
        {
          errors << priority(error) << "Value for '" << attribute << "'  attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  }
  else
  {
    errors << priority(warning) << "No attribute given for " << name_phrase
           << "." << endl;
  }
  
  return return_value;
}

bool
parser::prev_attribute_invalid(double value, const std::string& name)
{
  if(!finite(value))
  {
    errors << priority(critical) << "Attribute, " << name << ", in prev_constraint "
           << "sub-block missing or invalid.  Skipping analysis ..." << endl;
    return true;
  }
  else
  {
    return false;
  }
}

void
parser::parse_prev_constraint(const LSFBase* param)
{
  if(param->List())
  {
    psm_builder builder(&my_model.mean_cov_sub_model,
                        &my_model.susc_cov_sub_model,
                        &my_model.var_cov_sub_model);

    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "COVARIATE")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value() && a.String().size() > 0)
        {
          // Get the covariate string

          string covariate = a.String();
          double val       = get_real_attribute(*iter, "val", "covariate", "prev_constraint");

          psm_builder::add_cov_ret_type rt = builder.add_covariate(covariate,val);

          switch(rt)
          {
            case psm_builder::valid     :
            case psm_builder::duplicate : break;

            case psm_builder::unknown   :
            {
              errors << priority(critical) << "covariate '" << covariate 
                     << "' specified in prev_constraint sub-block must also be a specifed in "
                     << "mean_cov sub-block, suscept_cov sub-block or var_cov sub-block. "
                     << "Skipping analysis ..." << endl;
              my_model.m_class = model_INVALID;
              break;
            }

            case psm_builder::value_err :
            {
              // Value errors may occur, but if so, have already been
              // reported by the get_real_attribute function, so we can
              // safely break here.  There is no need to report them.

              break;
            }
          }
        }
        else // Covariate not specified
        {
          errors << priority(warning) << "Covariate name " 
                 << "not specified in prev_constraint sub-block.  "
                 << "Covariate will be ignored." << endl;
        }
      }

      else if(param_name == "AGE")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value() && a.String().size() > 0)
        {
          // Get the covariate string

          string value = toUpper(a.String());
          
          double val = QNAN;

          if(value == "INFINITY")
            val = POSITIVE_INF;
          else
            val = get_real_attribute(*iter, "value", "age", "prev_constraint");

          size_t rt = builder.set_age(val);

          if(rt)
          {
              errors << priority(critical)
                     << "Age specified in prev_constraint sub-block must be a non-negative"
                     << " number or 'INFINITY'.  "
                     << "Skipping analysis ..." << endl;
              my_model.m_class = model_INVALID;
              break;
          }
        }
        else // Covariate not specified
        {
          errors << priority(warning) << "Age value " 
                 << "not specified in prev_constraint sub-block.  "
                 << "Age will be ignored." << endl;
        }
      }

      else if(param_name == "R")
      {

        AttrVal  a = attr_value(*iter, 0);
        if(a.has_value())
        {
          if(finite(a.Real()) && a.Real() >= 0)
          {
            size_t err = builder.set_number_affected(a.Real());

            if(err)
              errors << priority(error) << "Value for 'R' parameter must be a finite value "
                     << "greater than 0.  Ignoring ..." << endl; 
          }
          else
          {
            errors << priority(error) << "Value for 'R' parameter in prev_constraint sub-block"
                   << " not understood.  Ignoring ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "Value for 'R' parameter in prev_constraint sub-block"
                 << " missing.  Ignoring ..." << endl;
        }

      }
      else if(param_name == "N")
      {

        AttrVal  a = attr_value(*iter, 0);
        if(a.has_value())
        {
          if(finite(a.Real()) && a.Real() >= 0)
          {
            size_t err = builder.set_sample_size(a.Real());

            if(err)
              errors << priority(error) << "Value for 'N' parameter must be a finite value "
                     << "greater than 0.  Ignoring ..." << endl; 
          }
          else
          {
            errors << priority(error) << "Value for 'N' parameter in prev_constraint sub-block"
                   << " not understood.  Ignoring ..." << endl;
          }
        }
        else
        {
          errors << priority(error) << "Value for 'N' parameter in prev_constraint sub-block"
                 << " missing.  Ignoring ..." << endl;
        }

      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in prev_constraint sub-block "
              << "not recognized.  Ignoring ..." << endl;        
      }
    }

    prevalence_sub_model::add_elt_ret_type rt =
        my_model.prev_sub_model.add_constraint(builder);

    switch(rt)
    {
      case prevalence_sub_model::valid : break;

      case prevalence_sub_model::na_not_finite :
      {
        errors << priority(critical) << "Attribute, R, in prev_constraint "
               << "sub-block missing or invalid.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::ss_not_finite :
      {
        errors << priority(critical) << "Attribute, N, in prev_constraint "
               << "sub-block missing or invalid.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::na_gt_ss :
      {
        errors << priority(critical) << "Number of affected persons, R, "
               << "may not be greater than sample size, N, in prev_constraint "
               << "sub-block.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::age_required :
      {
        errors << priority(critical) << "Age is required for prevalence constraints "
               << "when age of onset is used.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::age_present :
      {
        errors << priority(critical) << "Age may not be present in prevalence constraints "
               << "except when age of onset is used.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::duplicate :
      {
        errors << priority(warning) << "duplicate constraint detected in "
               << "prev_constraints sub-block.  Second copy will be ignored." << endl;
        break;
      }
    }
  }
  else
  {
    errors << priority(warning) << "prev_constraint sub-block contains no parameters."
         << endl;
  }
}

void  
parser::parse_prev_constraints(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "CONSTRAINT")
      {
        parse_prev_constraint(*iter);

        if(my_model.m_class == model_INVALID) break;
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in prev_constraints sub-block "
              << "not recognized.  Ignoring ..." << endl;        
      }
    }
  }
  else
  {
    errors << priority(warning) << "prev_constraints sub-block contains no parameters."
         << endl;
  }
}

// --- PREVALENCE ESTIMATE SUB-BLOCK ---
//
double  
parser::get_real_attribute(const LSFBase* param, const std::string& attribute,
                                                 const std::string& parameter, 
                                                 const std::string& sub_block )
{
  double  return_value = QNAN;
  string  name_phrase = parameter + " in " + sub_block + " sub-block";

  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    // - Check for non-meaningful attributes.
    //
    for(iter = a_list->begin(); iter != a_list->end(); ++iter)
    {
      string  attr_name = AttrNameMgr.query(iter->first);
      if(! (toUpper(attr_name) == "R"     ||
            toUpper(attr_name) == "N"     ||
            toUpper(attr_name) == "VAL"   ||
            toUpper(attr_name) == "VALUE"    ))
      {
        if(toUpper(attr_name) == "FIXED")
        {
          errors << priority(error) << "Attribute, fixed, not used in parameter, "
                 << name_phrase << ".  Ignoring ..." << endl;
        }
        else
        
        {
          errors << priority(error) << "Attribute, " << attr_name << ", of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  
    iter = a_list->find(toUpper(attribute));
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          return_value = a.Real();
        }
        else
        {
          errors << priority(error) << "Value for '" << attribute << "'  attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }
    }
  }
  else
  {
    errors << priority(warning) << "No attribute given for " << name_phrase
           << "." << endl;
  }
  
  return return_value;
}
    
void  
parser::parse_prev_estimate_sub_model(const LSFBase* param)
{
  if(param->List())
  {
    psm_builder builder(&my_model.mean_cov_sub_model,
                        &my_model.susc_cov_sub_model,
                        &my_model.var_cov_sub_model);

    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "COVARIATE")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value() && a.String().size() > 0)
        {
          // Get the covariate string

          string covariate = a.String();
          double val       = get_real_attribute(*iter, "val", "covariate", "prev_estimate");

          psm_builder::add_cov_ret_type rt = builder.add_covariate(covariate,val);

          switch(rt)
          {
            case psm_builder::valid     :
            case psm_builder::duplicate : break;

            case psm_builder::unknown   :
            {
              errors << priority(critical) << "covariate '" << covariate 
                     << "' specified in prev_estimate sub-block must also be a specifed in "
                     << "mean_cov sub-block, suscept_cov sub-block or var_cov sub-block. "
                     << "Skipping analysis ..." << endl;
              my_model.m_class = model_INVALID;
              break;
            }

            case psm_builder::value_err :
            {
              // Value errors may occur, but if so, have already been
              // reported by the get_real_attribute function, so we can
              // safely break here.  There is no need to report them.

              break;
            }
          }
        }
        else // Covariate not specified
        {
          errors << priority(warning) << "Covariate name " 
                 << "not specified in prev_estimate sub-block.  "
                 << "Covariate will be ignored." << endl;
        }
      }

      else if(param_name == "AGE")
      {
        AttrVal a = attr_value(*iter, 0);
        if(a.has_value() && a.String().size() > 0)
        {
          // Get the covariate string

          string value = toUpper(a.String());
          
          double val = QNAN;

          if(value == "INFINITY")
            val = POSITIVE_INF;
          else
            val = get_real_attribute(*iter, "value", "age", "prev_estimate");

          size_t rt = builder.set_age(val);

          if(rt)
          {
              errors << priority(critical)
                     << "Age specified in prev_estimate sub-block must be a non-negative"
                     << " number or 'INFINITY'.  "
                     << "Skipping analysis ..." << endl;
              my_model.m_class = model_INVALID;
              break;
          }
        }
        else // Covariate not specified
        {
          errors << priority(warning) << "Age value " 
                 << "not specified in prev_estimate sub-block.  "
                 << "Age will be ignored." << endl;
        }
      }

      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in prev_estimate sub-block "
              << "not recognized.  Ignoring ..." << endl;        
      }
    }

    prevalence_sub_model::add_elt_ret_type rt =
        my_model.prev_sub_model.add_estimate(builder);

    switch(rt)
    {
      case prevalence_sub_model::valid : break;

      case prevalence_sub_model::duplicate :
      {
        errors << priority(warning) << "duplicate estimate detected in "
               << "prev_estimates sub-block.  Second copy will be ignored." << endl;

      }

      case prevalence_sub_model::age_required :
      {
        errors << priority(critical) << "Age is required for prevalence estimates "
               << "when age of onset is used.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }

      case prevalence_sub_model::age_present :
      {
        errors << priority(critical) << "Age may not be present in prevalence estimates "
               << "except when age of onset is used.  Skipping analysis ..." << endl;
        my_model.m_class = model_INVALID;
        break;
      }
      
      case prevalence_sub_model::na_not_finite :
      case prevalence_sub_model::ss_not_finite :
      case prevalence_sub_model::na_gt_ss      :
        // These error cases are not applicable to the prevalence estimates
        break;

    }
  }
  else
  {
    errors << priority(warning) << "prev_estimate sub-block contains no parameters."
         << endl;
  }
}

// --- TYPE SUSCEPTIBILITY SUB-BLOCK ---
//    
void  
parser::parse_type_suscept_sub_model(const LSFBase* param)
{
  if(! TYPE_SUSCEPT_AVAILABLE)
  {
    errors << priority(critical) << "type_suscept sub-block not available in "
           << "this release.  Skipping analysis ..." << endl;
    my_model.m_class = model_INVALID;
    return;
  }
  else
  {
    parse_mean_sub_model(param, "suscept");
  }
}


// --- NOT IN SUB-BLOCKS ---
//
void
parser::parse_primary_trait(const LSFBase* param)
{
  string  value = "";
  parse_string(param, value);
    
  if(my_mp->info().trait_exists(value))   
  {
    my_model.primary_trait = value;
    
    // - Check trait type.
    //
    RPED::RefTraitInfo::trait_t  trait_t = primary_trait_type();
    if(! (trait_t == RPED::RefTraitInfo::binary_trait     ||
          trait_t == RPED::RefTraitInfo::continuous_trait   ))
    {
      my_model.m_class = model_INVALID;
      errors << priority(critical) << "Trait '" << value
             << "' must be continuous or binary.  Skipping analysis ..." << endl;
      return;    
    } 
    
    if(trait_t == RPED::RefTraitInfo::continuous_trait)
    {
      my_model.primary_trait_type = pt_CONTINUOUS;
    }
    else
    {
      my_model.primary_trait_type = pt_BINARY;
    }  
  }
  else
  {
    errors << priority(critical) << "Trait '" << value 
           << "' does not exist.  Skipping analysis ..." << endl;
    my_model.m_class = model_INVALID;
    return;
  }
  
  // - Parse type attribute.  
  //
  AttrList* a_list = param->attrs();
  AttrList::const_iterator  iter = a_list->find("TYPE");
  if(iter != a_list->end())
  {
    string  ptt = iter->second.String();
    if(toUpper(ptt) == "CONTINUOUS")
    {
      my_model.primary_trait_type = pt_CONTINUOUS;
    }
    else if(toUpper(ptt) == "BINARY")
    {
      my_model.primary_trait_type = pt_BINARY;
    }
    else if(toUpper(ptt) == "AGE_ONSET")
    {
      if(ONSET_AVAILABLE)
      {
        // Verify that the trait is binary (affection status)

        RPED::RefTraitInfo::trait_t  trait_t = primary_trait_type();
        if(trait_t == RPED::RefTraitInfo::binary_trait)
        {
          my_model.primary_trait_type = pt_ONSET;
          my_model.m_class = model_FPMM;
        }
        else
        {
          errors << priority(critical) << "Trait '" << value << "' must "
                 << "be binary for age of onset.  Skipping analysis ..."
                 << endl;
          my_model.m_class = model_INVALID;
        }
      }
      else
      {
        errors << priority(critical) << "Age of onset not available in this "
               << "version of SAGE.  Skipping analysis ..." << endl;
      }
    }
    else
    {
      errors << priority(error) << "Value of trait type attribute not "
             << "understood.  Ignoring ..." << endl;
    }
  }
}

void
parser::parse_title(const LSFBase* param)
{
  string  value = "";
  parse_string(param, value);
  if(value.size())
  {
    my_model.title = value;
  }
  else
  {
    errors << priority(warning) << "No value given for parameter, title.  "
           << "Using '" << my_model.title << "'." << endl; 
  }
}

// - New implemation.  For consistency among SAGE programs.
//
void
parser::parse_output_attribute(const LSFBase* param)
{
  AttrList* a_list = param->attrs();
  if(a_list)
  {
    string  out_value = a_list->StringAttr("OUT");
    string  output_value = a_list->StringAttr("OUTPUT");
    if(out_value.size())
    {
      my_model.file_name_root = out_value;
      file_name_root_parsed = true;
    }
    else if(output_value.size())
    {
      my_model.file_name_root = output_value;
      file_name_root_parsed = true;
    }
  }
}

// - As originally specified.
//
void
parser::parse_output_parameter(const LSFBase* param)
{
  if(! file_name_root_parsed)
  {
    string  value = "";
    parse_string(param, value);
    if(value.size())
    {
      my_model.file_name_root = value;
    }
  }
}



void
parser::type_warning(const string& type)
{
  errors << priority(warning) << "class '" << type << "' is incompatible "
         << "with a binary primary trait.  Ignoring ..." << endl;
}

void  
parser::parse_model_class(const LSFBase* param)
{
  string  value = "";
  parse_string(param, value);

  // Determine the model class as given by user

  if(toUpper(value) == "A")         my_model.m_class = model_A;
  else if(toUpper(value) == "D")    my_model.m_class = model_D;
  else if(toUpper(value) == "FPMM") my_model.m_class = model_FPMM;
  else if(toUpper(value) == "MLM") my_model.m_class = model_MLM; // for debugging purposes
  else if(value.size() == 0)        return;
  else
  {
    errors << priority(error) << "Invalid value for parameter, class, '"
           << value << "'.  Using model class, '" 
           << model_class_2_string(my_model.m_class) 
           << "'." << endl;

    return;
  }

  // Check for incompatibilities - This is only possible if the user gave a
  // class.  Note that MLM, since it is not set by the user (being
  // predicated on a binary trait) cannot be incompatible.

  switch(my_model.m_class)
  {
    case model_A :
    case model_D :
        switch(my_model.primary_trait_type)
        {
          case pt_BINARY :
            my_model.m_class = model_MLM;
            type_warning(toUpper(value));
            break;
          case pt_ONSET  :
            my_model.m_class = model_FPMM;
            errors << priority(warning) << "Age of onset requires FPMM model class."
                   << " Ignoring class '" << value << "' given." << endl;
            break;

          case pt_NONE       :

            SAGE_internal_error();
            break;

          case pt_CONTINUOUS :
          default            :
            break;
        }
        
        break;

    case model_FPMM :
        if(! FPMM_AVAILABLE)
        {
          errors << priority(critical) << "FPMM model class not available in "
                 << "this release.  Skipping analysis ..." << endl;
          my_model.m_class = model_INVALID;
        }
        break;

    case model_MLM     :
    case model_INVALID :
    default            :
        break;
  }
}
    
void  
parser::parse_output_options(const LSFBase* param)
{
  if(param->List())
  {
    LSFList::const_iterator  iter;
    for(iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "EACH_PEDIGREE")
      {
        if(! EACH_PEDIGREE_AVAILABLE)
        {
          errors << priority(error) << "each_pedigree option "
                 << "in output_options sub-block not available in "
                 << "this release.  Ignoring ..." << endl;
        }
        else
        {
          parse_boolean(*iter, my_model.each_pedigree);
        }
      }
      else if(param_name == "PEN_FUNC_OUT")
      {
        if(! PEN_FUNC_OUT_AVAILABLE)
        {
          errors << priority(error) << "pen_func_out option "
                 << "in output_options sub-block not available in "
                 << "this release.  Ignoring ..." << endl;
        }
        else
        {
          if(my_model.get_primary_trait_type() == pt_CONTINUOUS)
          {
            parse_boolean(*iter, my_model.pen_func_output);
          }
          else
          {
            errors << priority(error) << "pen_func_out option "
                   << "is only valid for continous traits.  Ignoring ..."
                   << endl;
          }
        }
      }
      else if(param_name == "TYPE_PROB")
      {
        if(! TYPE_PROBABILITIES_AVAILABLE)
        {
          errors << priority(error) << "type_prob option "
                 << "in output_options sub-block not available in "
                 << "this release.  Ignoring ..." << endl;
        }
        else
        {
          parse_boolean(*iter, my_model.type_prob);
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name
               << "' not recognized.  Ignoring ..." << endl; 
      }
    }
  }
}

// --- META-CONSTRAINTS ---
//
// - Check sub-model/sub-model and sub-model/top-level parameter
//   restrictions.
//
void
parser::check_meta_constraints()
{
  bitset<i_num_i> incon = my_model.check_consistency();
  
  // Print all discovered errors with this model
  for(int i = 0; i < i_num_i; ++i)
  {
    //lint --e{732}
    if(incon.test(i))
    {
      switch(i)
      {
        case i_transm_one_type:
          errors << priority(critical) << "Specified transmission sub-block option not allowed "
                 << "with only one type specified.  Skipping analysis ..." << endl;
          my_model.m_class = model_INVALID;
          break;
        
        case i_transm_nhwe:
           errors << priority(critical) << "transmission sub-block option is "
                  << "homog_no_trans, homog_mendelian or homog_general "
                  << "with option, nhwe, in geno_freq sub-block.  "
                  << "Skipping analysis ..." << endl;
          my_model.m_class = model_INVALID;
          break;
          
        case i_prim_cov:
          errors << priority(critical) << "Primary trait also used in mean_cov, var_cov, "
                 << "suscept_cov or composite trait sub-block.  Skipping analysis ..." << endl;
          my_model.m_class = model_INVALID;
          break;
         
        case i_transm_pen_out:
          errors << priority(error) << "pen_func_out in output_options only relevant "
                 << "when transmission sub-block option is homog_mendelian.  Ignoring ..." << endl;      
          my_model.pen_func_output = false;
          break;
         
        case i_type_prob_one_type:
          errors << priority(error) << "type_prob in output_options only relevant "
                 << "when type_mean or type_suscept sub-block option is not 'one'.  Ignoring ..." << endl;
          break;
          
        case i_freq_not_none_one_type:
          assert(false);      // Parser should insure that this doesn't happen.
          break;
          
        case i_mv_types:
          errors << priority(critical) << "Options for mean and variance sub-models "
                 << "incompatible in terms of number of types.  Skipping analysis ..." << endl;      
          my_model.m_class = model_INVALID;
          break;
         
        case i_freq_nhwe_two_types:
          assert(false);      // Parser should insure that this doesn't happen.
          break;
          
        case i_transm_no_trans_two_types:
          assert(false);      // Parser should insure that this doesn't happen.
          break;

        case i_binary_nt_no_residuals:
          errors << priority(critical) << "The maximum likelihood for a MLM model "
                 << "with no residuals and no transmission is identical to that "
                 << "for a one susceptibility model, and there is an infinity of "
                 << "sets of maximum likelihood estimates that yield the same "
                 << "Maximum likelihood.  Analysis will be ignored."
                 << endl;
          my_model.m_class = model_INVALID;
          break;
        case i_mc_cov_non_exclusive:
          print_covariate_non_exclusive_error
              (my_model.mean_cov_sub_model,
               my_model.comp_trait_sub_model);
          my_model.m_class = model_INVALID;
          break;
        case i_ms_cov_non_exclusive:
          print_covariate_non_exclusive_error
              (my_model.mean_cov_sub_model,
               my_model.susc_cov_sub_model);
          my_model.m_class = model_INVALID;
          break;
        case i_mv_cov_non_exclusive:
          print_covariate_non_exclusive_error
              (my_model.mean_cov_sub_model,
               my_model.var_cov_sub_model);
          my_model.m_class = model_INVALID;
          break;
        case i_sc_cov_non_exclusive:
          print_covariate_non_exclusive_error
              (my_model.susc_cov_sub_model,
               my_model.comp_trait_sub_model);
          my_model.m_class = model_INVALID;
          break;
        default:
          assert(false);
      }
    }
  }

  if(model_valid())
  {
    check_ignore_3();
  }
  
  if(model_valid())
  {
    check_ignore_4();
  }
  
  if(model_valid())
  {
    check_ignore_6();
  }
}

void
parser::print_covariate_non_exclusive_error
    (const CovariateSubmodel& sm1, const CovariateSubmodel& sm2) const
{
  list<string> shared_traits = get_shared_covariates(sm1, sm2);
  
  errors << priority(critical) 
         << sm1.get_subblock_name() << " sub-block and "
         << sm2.get_subblock_name() << " sub-block share trait(s) (";
         
  for(list<string>::const_iterator i = shared_traits.begin(); i != shared_traits.end(); ++i)
  {
    if(i != shared_traits.begin())
      errors << ", ";
    
    errors << *i;
  }

  errors << ").  These two sub-blocks must not share traits for this "
         << "model.  Analysis will be ignored."
         << endl;
}

void
parser::no_primary_trait()
{
  errors << priority(critical) << "No primary trait specified.  Skipping analysis ..."
         << endl;
  my_model.m_class = model_INVALID;
}

// - fpmm sub-model relevant only for model class fpmm.
//
void
parser::check_ignore_3()
{
  if(fpmm_parsed && my_model.m_class != model_FPMM)
  {
    errors << priority(error) << "fpmm sub-block only relevant when "
           << "value of class parameter is 'fpmm'.  Ignoring ..." << endl;
  }
} 

// - residual correlation sub-model not relevant for model class fpmm.
//
void
parser::check_ignore_4()
{
  if(resid_parsed && my_model.m_class == model_FPMM)
  {
    errors << priority(error) << "resid sub-block not relevant when "
           << "value of class parameter is 'fpmm'.  Ignoring ..." << endl;
  }
} 

// - transformation sub-model not relevant for discrete traits.
//
void
parser::check_ignore_6()
{
  if(transformation_parsed && my_model.primary_trait_type == pt_BINARY)
  {
    errors << priority(error) << "transformation sub-block not relevant for "
           << "binary primary trait.  Ignoring ..." << endl;
    my_model.transf_sub_model.set_option(MAXFUN::TransformationSubmodel::no_trans);
  }
}

//lint +esym(534,*parse_string*)
//lint +e774

}
}

#undef mean
#undef VAR
#undef susc
#undef comp
