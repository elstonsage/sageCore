//============================================================================
//
// File:      Parser.cpp                      
//                                                                          
// Author:    Stephen Gross
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//
//============================================================================

#include "util/StringUtils.h"
#include "util/RegexUtils.h"
#include "ageon/Parser.h"

namespace SAGE {
namespace AO {

//============================================================================
//
// Parser(...) CONSTRUCTOR
//
//============================================================================
Parser::Parser(
	const RPED::RefMultiPedigree & RMP, 
	      ostream                & m, 
	      cerrorstream           & errors) 
	:
	BasicParser (errors), 
	my_RMP      (RMP), 
	my_messages (m), 
	my_errors   (errors), 
	idnum       (0)
{
  my_messages << endl 
              << "Parsing Ageonset Analyses..." 
              << endl 
              << endl;

  reset();
}

//============================================================================
//
//  parse_test_parameter_section(...)
//
//============================================================================
void  
Parser::parse_test_parameter_section(const LSFBase* params)
{
  reset              ();
  init_parse         ();
  process_parameters (params);
  print_footer       ();
}

//============================================================================
//
//  init_parse()
//
//============================================================================
void Parser::init_parse()
{
  my_messages << "Beginning new analysis block...." 
              << endl << endl;

  my_model.title() = AO_DEFAULT_ANALYSIS_TITLE + long2str(++idnum);

  print_title();
}

//============================================================================
//
//  print_title()
//
//============================================================================
void Parser::print_title()
{
  my_messages << "Parsing Analysis: '" 
              << my_model.get_title() 
              << "'...."
              << endl << endl;
}

//============================================================================
//
//  print_footer()
//
//============================================================================
void Parser::print_footer()
{
  my_messages << endl
              << "Analysis parsing complete.  ";

  if( my_model.is_valid() )
  { 
    my_messages << "  Analysis valid." << endl;
  }
  else
    my_messages << "  Analysis is invalid. Check your file." << endl;

  my_messages << "=============================================================="
              << endl 
              << endl;
}

//============================================================================
//
//  reset()
//
//============================================================================
void
Parser::reset()
{
  my_model.reset();
}

//============================================================================
//
//  process_parameters(...)
//
//============================================================================
void
Parser::process_parameters(const LSFBase* params)
{
  // 1. Loop through parameters and process them:

  for(LSFList::const_iterator iter  = params->List()->begin (); 
                              iter != params->List()->end   (); ++iter)
  {
    if(*iter)
    {
      string param_name = toUpper((*iter)->name());

           if(param_name == "TITLE")              parse_title               (*iter);
      else if(param_name == "DEBUG")              my_model.debug() = true;
      else if(param_name == "MAXFUN")             MAXFUN::parseDebugParams(my_model.getDebugCfg(), *iter, my_errors);
      else if(param_name == "AFFECTEDNESS")       parse_affectedness_trait  (*iter);
      else if(param_name == "AGE_OF_ONSET" ||
              param_name == "AGE_ONSET" ) parse_age_of_onset_trait  (*iter);
      else if(param_name == "AGE_OF_EXAM" ||
              param_name == "AGE_EXAM" )   parse_age_of_exam_trait   (*iter);
      else if(param_name == "SUSCEPT_COV")        parse_covariate_sub_block (*iter, type_suscept);
      else if(param_name == "MEAN_COV")           parse_covariate_sub_block (*iter, type_mean);
      else if(param_name == "VAR_COV")            parse_covariate_sub_block (*iter, type_var);
      else if(param_name == "TRANSFORMATION")     parse_transform_sub_block (*iter);
      else if(param_name == "POOL")               parse_pool                (*iter);
      else if(param_name == "ALLOW_AVERAGING" ||
              param_name == "AA")                 parse_allow_averaging     (*iter); 

      else errors << priority(error) 
                  << "Parameter '" 
                  << param_name
                  << "' not recognized.  Skipping ..." 
                  << endl;
    }
  }

  if(    my_model.get_affectedness_trait() != ""
      && my_model.get_age_of_onset_trait() != ""
      && my_model.get_age_of_exam_trait() != "" )
  {
    my_model.ofilename() = "ageon_analysis" + long2str(idnum);

    if(params->attrs())
    {
      for(AttrList::const_iterator attr_itr  = params->attrs()->begin ();
                                   attr_itr != params->attrs()->end   (); ++attr_itr)
      {
        if(toUpper(AttrNameMgr.query(attr_itr->first)) == "OUT" && (attr_itr->second.String()) != "")
          my_model.ofilename() = attr_itr->second.String();
      }
    }

    my_model.verify();
  }
  else
  {
    errors << priority(error) 
           << "In analysis '" 
           << my_model.get_title()
           << "': Parameters 'affectedness', 'age_of_onset', and 'age_of_exam' are required.  "
           << "Can't continue analysis.  Skipping..."
           << endl; 
  }
}

//============================================================================
//
// parse_title(...)
//
//============================================================================
void
Parser::parse_title(const LSFBase* param)
{
  // 0. Set up local variables:

  string value = "";

  // 1. Fetch string value:

  parse_string(param, value);

  // 2. Check it and assign it to the model:

  if(value.size())
    my_model.title() = value;
  else
    errors << priority(warning) 
           << "In analysis '" 
           << my_model.get_title()
           << "': No value given for parameter, title.  "
           << "Using '" 
           << my_model.get_title() 
           << "'." 
           << endl; 
}

//============================================================================
//
// parse_affectedness_trait(...)
//
//============================================================================
void
Parser::parse_affectedness_trait(const LSFBase* param)
{
  string value = parse_trait(param, "affectedness");

  if( value.size() )
    my_model.affectedness_trait() = value;
  else
    errors << priority(error) 
           << "In analysis '" 
           << my_model.get_title()
           << "': Invalid value for parameter 'affectedness'.  "
           << endl; 
}

//============================================================================
//
// parse_age_of_onset_trait(...)
//
//============================================================================
void
Parser::parse_age_of_onset_trait(const LSFBase* param)
{
  string value = parse_trait(param, "age_of_onset");

  if( value.size() )
    my_model.age_of_onset_trait() = value;
  else
    errors << priority(error) 
           << "In analysis '" 
           << my_model.get_title()
           << "': Invalid value for parameter 'age_of_onset'.  "
           << endl; 
}

//============================================================================
//
// parse_age_of_exam_trait(...)
//
//============================================================================
void
Parser::parse_age_of_exam_trait(const LSFBase* param)
{
  string value = parse_trait(param, "age_of_exam");

  if( value.size() )
    my_model.age_of_exam_trait() = value;
  else
    errors << priority(error) 
           << "In analysis '" 
           << my_model.get_title()
           << "': Invalid value for parameter 'age_of_exam'.  "
           << endl; 
}

//============================================================================
//
// parse_trait(...)
//
//============================================================================
string
Parser::parse_trait(const LSFBase * param, string a_name)
{
  AttrVal a = attr_value(param,0);

  if( a.has_value() )
  {
    string tname = strip_ws(a.String());
    size_t t = my_RMP.info().trait_find(tname);

    if( t < my_RMP.info().trait_count() )
    {
      return tname;
    }
    else
    {
      errors << priority(error)
             << "In analysis '" 
             << my_model.get_title()
             << "': " << a_name << " parameter value '" << tname << "' not found." << endl;
    }

  }
  else
  {
    errors << priority(error)
           << "In analysis '" 
           << my_model.get_title()
           << "': " << a_name << " parameter is missing name." << endl;
  }

  return "";
}
//============================================================================
//
// parse_covariate_sub_block(...)
//
//============================================================================
void
Parser::parse_covariate_sub_block(const LSFBase* param, trait_type t)
{
  // 0. Set up local variables:

  string name = "";

  // 1. Loop through parameters:

  for(LSFList::const_iterator iter  = param->List()->begin(); 
                              iter != param->List()->end  (); ++iter)
  {
    if(*iter) 
    {
      name = toUpper((*iter)->name());

      if( name == "COVARIATE" ) 
        parse_covariate(*iter, t);
      else
      {
        errors << priority(critical)
               << "In analysis '" 
               << my_model.get_title()
               << "': Unrecognized value parameter '"
               << name
               << "' exists in parameter file." 
               << endl;

        exit(0);
      }
    }
  }
}

//============================================================================
//
// parse_covariate(...)
//
//============================================================================
void 
Parser::parse_covariate(const LSFBase * param, trait_type t)
{
  // 0. Set up local variables:

  AttrVal attr_val;
  string  attr_name       = "",
          name            = param->name();
  bool    new_initial_est = false,
          new_fixed       = false,
          fixed           = false;
  double  initial_est     = 1.0;

  // 1. Process user-specified options:

  for(AttrList::const_iterator attr_itr  = param->attrs()->begin ();
                               attr_itr != param->attrs()->end   (); ++attr_itr)
  {
    attr_name = toUpper(AttrNameMgr.query(attr_itr->first));
    attr_val  = attr_itr->second;

         if(attr_name == "VALUE") { name            = attr_val.String(); }
    else if(attr_name == "VAL")   { new_initial_est = true; 
                                    initial_est     = attr_val.Real(); }
    else if(attr_name == "FIXED") { new_fixed       = true; 
                                    fixed           = toUpper(attr_val.String()) == "YES" || 
                                                      toUpper(attr_val.String()) == "TRUE" ? true : false; }
    else                          
    {
      errors << priority(critical) 
             << "In analysis '" 
             << my_model.get_title()
             << "': Unrecognized value parameter '"
             << attr_name 
             << "' exists in parameter file." 
             << endl; 
      exit(0); 
    }
  }

  // 2. Override the name variable if we are parsing a transformation parameter (lambda1/lambda2):

  string comp_name = toUpper(name);

  if(comp_name == toUpper(my_model.GetParameterMgr().getParameter("Transformation", "Lambda1").getName()))
    name = my_model.GetParameterMgr().getParameter("Transformation", "Lambda1").getName();

  if(comp_name == toUpper(my_model.GetParameterMgr().getParameter("Transformation", "Lambda2").getName()))
    name = my_model.GetParameterMgr().getParameter("Transformation", "Lambda2").getName();

  // 3. Add the new parameter (or update the existing parameter) to the model:

  my_model.add_trait(t, name, new_initial_est, initial_est, new_fixed, fixed);
}

//============================================================================
//
// parse_transform_sub_block(...)
//
//============================================================================
void
Parser::parse_transform_sub_block(const LSFBase* param)
{
  // 0. Set up local variables:

  string name;

  // 1. Loop through parameters:

  for(LSFList::const_iterator iter  = param->List()->begin (); 
                              iter != param->List()->end   (); ++iter)
  {
    if(*iter)
    {
      name = toUpper((*iter)->name());

      if(name == "LAMBDA1" || name == "LAMBDA2")
        parse_covariate(*iter, type_var);
      else
      {
        errors << priority(critical)
               << "In analysis '" 
               << my_model.get_title()
               << "': Unrecognized value parameter '"
               << name
               << "' exists in parameter file."
               << endl;

        exit(0);
      }
    }
  }
}

//============================================================================
//
// parse_allow_averaging(...)
//
//============================================================================
void
Parser::parse_allow_averaging(const LSFBase* param)
{
  // 0. Verify parameter:

  if(!*param)
    return;
                  
  // 1. Set up local variables:

  bool aa = false;
                  
  // 2. Fetch string value, convert to upper case:

  parse_boolean(param, aa);

  // 3. Check AA value:

  if( aa ) my_model.allow_averaging() = true;
  else     my_model.allow_averaging() = false;
}

//============================================================================
//
// parse_pool(...)
//
//============================================================================
void
Parser::parse_pool(const LSFBase* param)
{
  if(param->attrs() && param->attrs()->has_attr("value"))
  {
    // Set up variables;
    std::string              pool_str    = param->attrs()->find("value")->second.String();
    std::vector<std::string> assignments;
    
    // Split along the comma into assignment statements:
    UTIL::StringUtils::splitString(pool_str, ",", assignments);

    // Loop across the assignment statements:
    for(size_t i = 0; i < assignments.size(); ++i)
    {
      // Split along whitespace and '=" into tokens:
      std::vector<std::string> tokens;
      
      UTIL::StringUtils::splitMultiDelimitedString(assignments[i], " \t=", tokens);
      
      // We need at least TWO tokens:
      if(tokens.size() < 2)
      {
        errors << priority(warning)
               << "In analysis '" 
               << my_model.get_title()
               << "': Pool statements must consist of at least two tokens (eg: \"?U=AU\"). Ignoring badly formed portion of pool statement ('"
               << assignments[i]
               << "')..."
               << std::endl;

        continue;
      }
      
      // Fetch the target code:
      size_t to_code = PoolingCfg::getClass(tokens[tokens.size() - 1]);
      
      // Make sure it's valid:
      if(to_code == (size_t)-1)
      {
        errors << priority(warning)
               << "In analysis '" 
               << my_model.get_title()
               << "': Ignoring badly formed reassignment code ('"
               << tokens[tokens.size() - 1]
               << "') in portion of pool statement ('"
               << assignments[i]
               << "')..."
               << std::endl;

        continue;
      }
      
      // Try to reassign stuff:
      for(size_t j = 0; j < tokens.size() - 1; ++j)
      {
        // Fetch from source code:
        size_t from_code = PoolingCfg::getClass(tokens[j]);
        
        // Make sure it's valid:
        if(from_code == (size_t)-1)
        {
          errors << priority(warning)
                 << "In analysis '" 
                 << my_model.get_title()
                 << "': Ignoring badly formed reassignment code in portion ('"
                 << tokens[j]
                 << "') of pool statement ('"
                 << assignments[i]
                 << "')..."
                 << std::endl;

          continue;
        }
        
        // Try to carry out the reasssignment:
        if( my_model.get_pooling_cfg().setPool(from_code, to_code) )
        {
          my_model.pool_class() = true;

          my_messages << "Reassigning individuals from class '" 
                      << PoolingCfg::getCode(from_code) 
                      << "' to '"
                      << PoolingCfg::getCode(to_code)
                      << "'."
                      << std::endl;
        }
        else // Invalid conversion
        {
          errors << priority(warning)
                 << "In analysis '" 
                 << my_model.get_title()
                 << "': Ignoring illogical reassignment from class '"
                 << PoolingCfg::getCode(from_code) 
                 << "' to '"
                 << PoolingCfg::getCode(to_code)
                 << "'."
                 << std::endl;

          continue;
        }
      }
    }
  }
  else // Pool doesn't have a value attribute!
  {
    errors << priority(warning)
           << "In analysis '" 
           << my_model.get_title()
           << "': Pool statement requires a value (ie: 'pool=\"?U=AU\"'). Ignoring..."
           << std::endl;
  }
}

} // End namespace AO
} // End namespace SAGE
