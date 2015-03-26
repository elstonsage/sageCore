//=============================================================================
//
// File:      Parser.cpp                      
//                                                                          
// Author:    Stephen Gross
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//=============================================================================

#include "assoc/Parser.h"
#include "func/Expression.h"

namespace SAGE  {
namespace ASSOC {

Configuration 
Parser::parseAssocBlock(size_t analysis_num, const LSFBase* params, 
                        const RPED::MultiPedigree& mp,     
                        std::ostream& info, cerrorstream& errors)
{
  Configuration config;
  
  bool  polygenic_eff = true;
  bool  family_eff = false;
  bool  marital_eff = true;
  bool  sibling_eff = true;
  
  info << "Beginning new analysis block....\n" << endl;

  if(! params->List())
  {
    errors << priority(critical) << "User has specified an analysis block with no content. Exiting..." << std::endl;
  
    throw std::exception();
  }

  bool  batch_mode = false;      // Added 10/5/7.  -djb
  vector<LSFList::const_iterator>  covariates;
  LSFBase*  trans_param = 0;
  LSFBase*  resid_param = 0;
  
  for(LSFList::const_iterator iter = params->List()->begin(); iter != params->List()->end(); ++iter)
  {
    if(*iter)
    {
      string param_name = toUpper((*iter)->name());

      if(param_name == "TITLE")              
      {
        parseTitle(config, *iter, errors);
      }
      else if(param_name == "BATCH")
      {
        batch_mode = true;       // Must be parsed after covariates.
      }
      /*
      else if(param_name == "SCORE_TEST_ONLY")
      {
        config.setScoreTestOnly(true);
      }
      */
      else if(param_name == "MAXFUN")             
      {
        config.setMaxfunDebug(true);
      }
      else if(param_name == "REVERSE_SORT")             
      {
        config.reverse_sort = true;
      }
      else if(param_name == "PRIMARY_TRAIT" || param_name == "TRAIT")              
      {
        parsePrimaryTrait(config, mp, *iter, errors);
      }
      else if(param_name == "COVARIATE" || param_name == "COV")
      {
        covariates.push_back(iter);  // Parse covariates once primary trait is known.
      }
      else if(param_name == "FAMILY_EFFECT" || param_name == "FE")                 
      {
        parseFixedEffect(*iter, config.getFixedEffect(Configuration::FAMILY), errors, family_eff);
      }
      else if(param_name == "MARITAL_EFFECT" || param_name == "ME")
      {
        parseFixedEffect(*iter, config.getFixedEffect(Configuration::MARITAL), errors, marital_eff);
      }
      else if(param_name == "SIBLING_EFFECT" || param_name == "SIBSHIP_EFFECT" || param_name == "SE")
      {
        parseFixedEffect(*iter, config.getFixedEffect(Configuration::SIBLING), errors, sibling_eff);
      }
      else if(param_name == "POLYGENIC_EFFECT" || param_name == "PE")
      {
        parseFixedEffect(*iter, config.getFixedEffect(Configuration::POLYGENIC), errors, polygenic_eff);
      }
      else if(param_name == "CLASS_EFF")
      {
        parseUserEffect(*iter, config, mp.info(), errors);
      }      
      else if(param_name == "ALLOW_AVERAGING" || param_name == "AA")
      {
        parseAllowAveraging(config, *iter, info, errors);
      }
      else if(param_name == "TRANSFORMATION" || param_name == "TRANSFORM" || param_name == "TRANS")
      {
        trans_param = *iter;
      }
      else if(param_name == "OMIT_COMPLETE_SUMMARY")
      {
        config.setOmitCompleteSummary(true);
      }      
      else if(param_name == "SUMMARY_DISPLAY")
      {
        parseSummaryDisplay(config, *iter, errors);
      }
      else if(param_name == "RESIDUALS")
      {
        resid_param = *iter;
      }      
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' not recognized.  Skipping ..."  << endl;
      }
    }
  }
  
  if(trans_param)
  {
    parseTransSubModel(config, trans_param, info, errors);
  }
  
  size_t  covariate_count = covariates.size();
  for(size_t c = 0; c < covariate_count; ++c)
  {
    parseCovariate(config, mp, *(covariates[c]), info, errors);
  }
  
  if(batch_mode)
  {
    parseBatch(config, mp, errors);
  }
  
  if(config.getModelList().size() > 1)
  {
    config.removeBaselineModel();
  }

  // Make sure the title is non-empty, otherwise set it to something:
  
  if(config.getTitle() == "")
  {
    std::ostringstream t; t << "assoc"; if(analysis_num) { t << analysis_num; }
    config.setTitle(t.str());
  }

  // 3. Set the ofilename to be equal to the model name:

  config.setOfilename(config.getTitle());

  // 4. Now check to see if the user specified an output filename option to override the model name:

  if(params->attrs())
  {
    for(AttrList::const_iterator attr_itr = params->attrs()->begin(); attr_itr != params->attrs()->end(); ++attr_itr)
    {
      if(toUpper(AttrNameMgr.query(attr_itr->first)) == "OUT" && (attr_itr->second.String()) != "")
      {
        config.setOfilename(attr_itr->second.String());
      }
    }
  }
  
  //prohibitBinaryTransformation(config, errors);           // Added 6-6-7. djb

  config.initializeModels(polygenic_eff, sibling_eff, marital_eff, family_eff);
  
  if(resid_param)
  {
    parseResiduals(config, resid_param, errors);
  }
      
//#define DEBUG_PARSER 
#ifdef DEBUG_PARSER
  cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << config.dump() << endl;
  Configuration::const_effect_iterator  ue_iter = config.userEffectBegin();
  Configuration::const_effect_iterator  ue_end_iter = config.userEffectEnd();
  for(; ue_iter != ue_end_iter; ++ue_iter)
  {
    cout << ue_iter->dump() << endl;
  }
  
  config.dumpCovariates();
  config.dumpModels();
  cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  exit(0);
  
#endif

  /*
  cout << config.getTransConfig().dump();
  exit(0);
  */
  
  return config;
}


void
Parser::parseTitle(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  string value = "";

  APP::ParsingFunctions::parse_string(param, value, false);

  if(value.size()) 
  {
    config.setTitle(value);
  }
  else
  {
    errors << priority(warning) << "No value given for parameter, title.  Using '" << config.getTitle() << "'." << std::endl; 
  }
}


void
Parser::parsePrimaryTrait(Configuration& config, const RPED::MultiPedigree& mp, const LSFBase* param, cerrorstream& errors)
{
  // 2. Grab string value:

  string value = "";

  APP::ParsingFunctions::parse_string(param, value, false);

  // 3. Check for a valid value, and assign it:

  if(!value.size())
  {
    errors << priority(warning) << "No value given for parameter, primary_trait.  Using '" << config.getPrimaryTraitName() << "'."  << endl; 
    return;
  }

  config.setPrimaryTraitName(value);

  // Make sure the trait exists!
  if(!mp.info().trait_exists(config.getPrimaryTraitName()))
  {
    errors << priority(error) << "Primary trait '" << config.getPrimaryTraitName() << "' not in the pedigree data; Skipping this analysis..." << endl;
    
    throw std::exception();
  }

  // 4. Check for binary trait use:

  if(mp.info().get_trait_type(config.getPrimaryTraitName()) == RPED::RefTraitInfo::binary_trait)
  {
    config.setDependentTraitType(BINARY);
  }
}


void
Parser::parseBatch(Configuration& config, const RPED::MultiPedigree& mp, cerrorstream& errors)
{
  // Loop across all COVARIATES (per RCE 11/8) and add each one as a unique covariate to a unique model
  // in the configuration if it's not already in the configuration as a covariate
  // in all models.
  for(size_t trait_id = 0; trait_id < mp.info().trait_count(); ++trait_id)
  {
    const string& cov_name = mp.info().trait_info(trait_id).name();
    
    if(! config.hasNullCovariate(cov_name)                                           && 
       toUpper(config.getPrimaryTraitName()) != toUpper(cov_name)                    &&
       ! modelNameInUse(cov_name, config, errors)                                    &&
       mp.info().trait_info(trait_id).usage() == RPED::RefTraitInfo::trait_covariate    )
       
    {
      Configuration::CovariateInfo cinfo;
      Configuration::ModelList  model_list(1, cov_name);
      cinfo.cfg.param_name = cov_name;
      config.addCovariate(cinfo, model_list);
    }
  }
}

bool
Parser::modelNameInUse(const string& name, Configuration& config, cerrorstream& errors)
{
  if(config.hasModel(name))
  {
    errors << priority(warning) << "Not creating model '" << name << "' in batch mode.  "
           << "A model by this name has already been specified." << endl;
  
    return  true;
  }
  else
  {
    return  false;
  }
}


void
Parser::parseFixedEffect(const LSFBase* lsf_param, Parameter& param,  
                                 cerrorstream& errors, bool& use_effect)
{
  AttrVal  attr_val;
  string  value;
  string  value_str;
  string  attr_name;

  // - Include in model?
  //
  APP::ParsingFunctions::parse_string(lsf_param, value, false);
  if(value != "")
  {
    bool  in_use = use_effect;
    
    //APP::ParsingFunctions::parse_boolean(lsf_param, in_use, false);
    
    APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_boolean(lsf_param, in_use, false);
    if(status == APP::LSFConvert::GOOD) 
    {   
        use_effect = in_use;
    }
    else
    {
      errors << priority(error) << "Value of parameter '" << lsf_param->name() << "' not understood.  "
             << "Ignoring this parameter ..." << endl;
             
      return;
    }
  }

  // - Initial value and "fixity"
  //
  assert(isnan(param.initial_estimate));
  assert(param.initial_type == MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL);
  if(lsf_param->attrs())
  {
    for(AttrList::const_iterator attr_itr = lsf_param->attrs()->begin(); attr_itr != lsf_param->attrs()->end(); ++attr_itr)
    {
      attr_name = toUpper(AttrNameMgr.query(attr_itr->first));
      attr_val  = attr_itr->second;

      if(attr_name == "VAL")   
      { 
        param.initial_estimate = attr_val.Real();
      }
      else if(attr_name == "FIXED") 
      { 
        param.initial_type = toUpper(attr_val.String()) == "TRUE" ? MAXFUN::Parameter::FIXED : 
                                                                           MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      }
    }
  }
  
  // - Specifing a fixed variance component without a value is an error.
  //
  if(param.initial_type == MAXFUN::Parameter::FIXED && isnan(param.initial_estimate))
  {
    errors << priority(error) << "Variance component '" << param.param_name << "' designated as FIXED, "
           << "but no value is given.  This variance component will not be included in the model." << endl;
    use_effect = false;
    
    return;
  }

  // - Remove variance component from model if fixed at 0.
  //
  if(param.initial_estimate == 0.0                  && 
     param.initial_type == MAXFUN::Parameter::FIXED && 
     use_effect == true)
  {
    errors << priority(information) << "User has set '" << param.param_name << "' to 0, and has also set fixed to TRUE.\n"
           << "Please note that this is equivalent to specifying '" << param.param_name << "=FALSE'.\n" << endl;

    use_effect = false;
  }
  
  // - Ensure that initial estimate is w/i bounds.
  //
  if(use_effect == true              && 
     ! isnan(param.initial_estimate) && 
     param.initial_estimate < param.lower_bound)
  {
    if(param.initial_type == MAXFUN::Parameter::FIXED)
    {
      errors << priority(critical) << "Fixed, non-zero value less than " << param.lower_bound
             << " not allowed for variance component '" << param.param_name << "'." << endl;
             
      exit(0);    
    }
    else
    {
      errors << priority(error) << "Initial value less than " << param.lower_bound
             << " not allowed for variance component '" << param.param_name << "'.  "
             << "Setting initial value to " << 2 * param.lower_bound << " ..." << endl;
             
      param.initial_estimate = 2 * param.lower_bound;
    }
  }
}


void
Parser::parseUserEffect(const LSFBase* lsf_param, Configuration& config, const RPED::RefMPedInfo& rmp, cerrorstream& errors)
{
  AttrVal  attr_val;
  string  value;
  string  value_str;
  string  attr_name;
  
  Parameter  param("Variance components", "", 
                   MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, QNAN, VARIANCE_COMPONENT_LB, MAXFUN::MF_INFINITY);

  // - Initial value and "fixity"
  //
  if(lsf_param->attrs())
  {
    for(AttrList::const_iterator attr_itr = lsf_param->attrs()->begin(); attr_itr != lsf_param->attrs()->end(); ++attr_itr)
    {
      attr_name = toUpper(AttrNameMgr.query(attr_itr->first));
      attr_val  = attr_itr->second;

      if(attr_name == "VALUE")
      {
        param.param_name = attr_val.String();
        if(fixed_effect(param.param_name))
        {
          errors << priority(error) << "Categorical trait name '" << param.param_name << "' is a reserved word. "
                 << "Please rename this trait.  Ignoring parameter 'class_eff' ..." << endl;
                 
          return;        
        }
        
        if(! isCategoricalTrait(param.param_name, rmp))
        {
          errors << priority(error) << "Categorical trait '" << param.param_name << "' does not exist. "
                 << "Ignoring parameter 'class_eff' ..." << endl;
                 
          return;
        }
      }
      else if(attr_name == "VAL")   
      { 
        param.initial_estimate = attr_val.Real();
      }
      else if(attr_name == "FIXED") 
      { 
        param.initial_type = toUpper(attr_val.String()) == "TRUE" ? MAXFUN::Parameter::FIXED : 
                                                                           MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      }
    }
  }
  
  // - Specifing a fixed variance component without a value is an error.
  //
  if(param.initial_type == MAXFUN::Parameter::FIXED && isnan(param.initial_estimate))
  {
    errors << priority(error) << "Variance component '" << param.param_name << "' designated as FIXED, "
           << "but no value is given.  This variance component will not be included in the model." << endl;
    
    return;
  }

  // - Remove variance component from model if fixed at 0.
  //
  if(param.initial_estimate == 0.0                  && 
     param.initial_type == MAXFUN::Parameter::FIXED   )
  {
    errors << priority(information) << "User has set '" << param.param_name << "' to 0, and has also set fixed to TRUE.\n"
           << "Ignoring 'class_eff' parameter ...\n" << endl;

    return;
  }
  
  // - Ensure that initial estimate is w/i bounds.
  //
  if(! isnan(param.initial_estimate) && 
     param.initial_estimate < param.lower_bound)
  {
    if(param.initial_type == MAXFUN::Parameter::FIXED)
    {
      errors << priority(critical) << "Fixed, non-zero value less than " << param.lower_bound
             << " not allowed for variance component '" << param.param_name << "'." << endl;
             
      exit(0);    
    }
    else
    {
      errors << priority(error) << "Initial value less than " << param.lower_bound
             << " not allowed for variance component '" << param.param_name << "'.  "
             << "Setting initial value to " << 2 * param.lower_bound << " ..." << endl;
             
      param.initial_estimate = 2 * param.lower_bound;
    }
  }
  
  config.addUserEffect(param);
}


bool
Parser::isCategoricalTrait(const string& trait_name, const RPED::RefMPedInfo& rmp)
{
  return  rmp.trait_exists(trait_name) && rmp.get_trait_type(trait_name) == RPED::RefTraitInfo::categorical_trait;
}


void 
Parser::parseCovariate(Configuration& config, const RPED::MultiPedigree& mp, 
                       const LSFBase* lsf_param, std::ostream& info, cerrorstream& errors)
{
  string  attr_name = "";
  AttrVal  attr_val;
  Configuration::CovariateInfo  cinfo;
  
  // Set default options:
  Configuration::ModelList  model_list;
  
  cinfo.cfg.initial_estimate = SAGE::QNAN;
  assert(cinfo.cfg.initial_type == MAXFUN::Parameter::INDEPENDENT);
  
  if(lsf_param->attrs())
  {
    for(AttrList::const_iterator attr_itr = lsf_param->attrs()->begin(); attr_itr != lsf_param->attrs()->end(); ++attr_itr)
    {
      attr_name = toUpper(AttrNameMgr.query(attr_itr->first));
      attr_val  = attr_itr->second; 

      if(attr_name == "MODELS") 
      {
        UTIL::StringUtils::splitMultiDelimitedString(attr_val.String(), " \t\n,", model_list);
      }
      else if(attr_name == "TEST")   
      { 
        errors << priority(error) << "The 'covariate' parameter attribute 'TEST' has been replaced "
               << "by 'MODELS'. Please see the description of the 'MODELS' attribute in the current "
               << "User Document.  Ignoring attribute 'TEST' ..." << endl;
      }
      else if(attr_name == "VALUE")
      { 
        cinfo.cfg.param_name = attr_val.String();
        
        // - Check for duplicate covariate.
        //
        if(config.hasCovariate(cinfo.cfg.param_name))
        {
          errors << priority(warning) << "Covariate '" << cinfo.cfg.param_name << "' specified "
                 << "more than once; ignoring all but first instance ..." << endl;
          return;
        }
        
        // - Check that the covariate is not the primary trait.
        //
        if(toUpper(config.getPrimaryTraitName()) == toUpper(cinfo.cfg.param_name))
        {
          errors << priority(error) << "Covariate '" << cinfo.cfg.param_name << "' is the "
                 << "primary trait; ignoring this covariate ..." << endl;
          return;
        }
      }
      else if(attr_name == "VAL")    
      { 
        cinfo.cfg.initial_estimate = attr_val.Real();
      }
      else if(attr_name == "FIXED")  
      { 
        cinfo.cfg.initial_type = toUpper(attr_val.String()) == "TRUE" ? MAXFUN::Parameter::FIXED : 
                                                                               MAXFUN::Parameter::INDEPENDENT;
      }
      else
      { 
        errors << priority(critical) << "Unrecognized value parameter '" << attr_name << "' exists in parameter file." << endl; 

        exit(0); 
      }
    }
  }

  // Post-process a few things:
  
  // Check for error conditions:

  // Make sure the covariate is in the multipedigree:
  if(cinfo.cfg.param_name == "")
  {
    errors << priority(error) << "No name given for one of the covariates in the analysis; "
           << "ignoring covariate ..." << endl;
  }
  
  // Not in the multipedigree:
  else if(! mp.info().trait_exists(cinfo.cfg.param_name))
  {
    errors << priority(error) << "Covariate '" << cinfo.cfg.param_name 
           << "' not in the pedigree data; ignoring covariate ..." << endl;
  }
  
  // - Specifing a fixed variance covariate without specifying a value is an error.
  //
  else if(cinfo.cfg.initial_type == MAXFUN::Parameter::FIXED && isnan(cinfo.cfg.initial_estimate))
  {
    errors << priority(error) << "Covariate '" << cinfo.cfg.param_name << "' designated as FIXED, "
           << "but no value is given; ignoring covariate ..." << endl;
                               
    return;
  }  
  
  // - Test covariate may not be fixed.
  //
  else if(cinfo.cfg.initial_type == MAXFUN::Parameter::FIXED && 
          (! model_list.empty())                             )
  {
    errors << priority(error) << "Test Covariate '" << cinfo.cfg.param_name << "' designated as FIXED.  "
                                 "Ignoring covariate ..." << endl;
    return;
  }
  else
  {
    // - Baseline is a 'reserved word.'
    //
    vector<string>::iterator  found = find(model_list.begin(), model_list.end(), "Baseline");
    if(found != model_list.end())
    {
      model_list.erase(found);
      errors << priority(error) << "Name 'Baseline' given in a covariate model list.  This "
             << "is a reserved word.  Ignoring ..." << endl;
    }
    
    config.addCovariate(cinfo, model_list);
  }
}

void
Parser::parseTransSubModel(Configuration& config, 
                            const LSFBase* param, ostream& info, cerrorstream& errors)
{
  // - Default as set at Configuration construction: ge; transform diff; lambda1 = 1, not fixed;
  //   lambda2 = 0, fixed.
  //
  MFSUBMODELS::Transformation::Configuration&  trans_config = config.getTransConfig();
  
  if(param->List())
  {
    const LSFBase*  lambda1 = 0;
    const LSFBase*  lambda2 = 0;
    
    for(LSFList::const_iterator iter = param->List()->begin(); iter != param->List()->end(); ++iter)
    {
      string param_name = toUpper((*iter)->name());
            
      if(param_name == "OPTION")  
      {
        MFSUBMODELS::Transformation::Parser::parse_option(trans_config, *iter, info, errors);
        if(trans_config.get_type() == MFSUBMODELS::Transformation::BOX_COX)
        {
          errors << priority(critical) << "Box-Cox transformation not allowed in ASSOC.  Ending program ..." << endl;
          exit(1);
        }
      }
      else if(param_name == "BOTH_SIDES")
      {
        config.setTransformBothSides(true);
      }
      else if(param_name == "LAMBDA1") 
      {
        lambda1 = *iter;
      }
      else if(param_name == "LAMBDA2") 
      {
        lambda2 = *iter;
      }
    }
    
    if(config.getDependentTraitType() == BINARY && config.getTransformBothSides())       
    {      
      errors << priority(error) << "Transformation parameter 'both sides' may not "      
             << "be specified with a binary dependent trait.  Skipping this analysis ..." << endl;       
             throw  std::exception();       
    }  
    
    if(lambda1)
    {
      parseLambda1(trans_config, lambda1, info, errors);
    }
    
    if(lambda2)
    {
      parseLambda2(config, trans_config, lambda2, info, errors);         
    }    
  }
}


void
Parser::parseLambda1(MFSUBMODELS::Transformation::Configuration& trans_config, 
                     const LSFBase* param, ostream& info, cerrorstream& errors)
{
  if(! param->attrs())
  {
    errors << priority(warning) << "No attributes specified for parameter 'lambda1' of the transformation sub-block." << endl;
    return;
  }
                 
  AttrList::const_iterator attr_iter;
                       
  // - Fixity
  //  
  attr_iter = param->attrs()->find("FIXED");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(toUpper(attr_iter->second.String()) == "TRUE")
    {
      trans_config.set_lambda1_fixed(true);
    }
    else if(toUpper(attr_iter->second.String()) == "FALSE")
    {
      trans_config.set_lambda1_fixed(false);
    }
    else
    {   
      errors << priority(error) << "Value for 'fixed' attribute of 'lambda1' not understood.  Ignoring ..." << endl;
    }
  }
                                                                    
  // - Lower bound
  //
  attr_iter = param->attrs()->find("LOWER_BOUND");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    double  users_lower_bound = attr_iter->second.Real();
    if(finite(users_lower_bound))
    {
      // Default lower bound changed to -infinity by rce on 10-25-11.
      //
      if(users_lower_bound > -numeric_limits<double>::infinity() && users_lower_bound <= 1.0)
      {
        trans_config.set_lambda1_lower_bound(users_lower_bound);
      }
      else
      {
        errors << priority(error) << "Value for 'lower_bound' attribute of 'lambda1' outside of allowed range.  Setting to 1.0 ..." << endl;
        trans_config.set_lambda1_lower_bound(1.0);    
      }
    }
    else
    {
      errors << priority(error) << "Value for 'lower_bound' attribute of 'lambda1' not understood.  Ignoring ..." << endl;
    }
  }
  
  // - Upper bound
  //
  attr_iter = param->attrs()->find("UPPER_BOUND");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    double  users_upper_bound = attr_iter->second.Real();
    if(finite(users_upper_bound))
    {
      if(users_upper_bound >= 1)
      {
        trans_config.set_lambda1_upper_bound(users_upper_bound);
      }
      else
      {
        errors << priority(error) << "Value for 'upper_bound' of attribute of 'lambda1' outside of allowed range.  Setting to 1.0 ..." << endl;
        trans_config.set_lambda1_upper_bound(1.0);    
      }
    }
    else
    {
      errors << priority(error) << "Value for 'upper_bound' attribute of 'lambda1' not understood.  Ignoring ..." << endl;
    }
  }
  
  // - Check bounds against each other.
  //
  if(trans_config.get_lambda1().lower_bound >= trans_config.get_lambda1().upper_bound)
  {
    errors << priority(error) << "The value of 'lambda1' 'lower_bound' is greater than or equal to "
           << "the value of 'lambda1' 'upper_bound'.  Setting 'lower_bound' to -1 and 'upper_bound' "
           << "to infinity ..." << endl;
    trans_config.set_lambda1_lower_bound(-1.0);
    trans_config.set_lambda1_upper_bound(numeric_limits<double>::infinity());
  }
  
  // - Value
  //                           
  attr_iter = param->attrs()->find("VAL");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(trans_config.get_lambda1().initial_type != MAXFUN::Parameter::FIXED)
    {
      errors << priority(error) << "Initial value may not be specified for 'lambda1' if it is to be estimated."  
             << "  Ignoring ..." << endl;
      return;    
    }
    
    double  users_init_est = attr_iter->second.Real();
    if(finite(users_init_est))
    {
      if(trans_config.get_lambda1().lower_bound <= users_init_est &&
         users_init_est <= trans_config.get_lambda1().upper_bound)
      {
        trans_config.set_lambda1_init_est(attr_iter->second.Real());
      }
      else
      {
        errors << priority(error) << "Value specified for 'val' attribute of 'lambda1' is out of bounds.  Ignoring ..." << endl;
      }
    }
    else
    {   
      errors << priority(error) << "Value for 'val' attribute of 'lambda1' not understood.  Ignoring ..." << endl;
    }
  }
}


void
Parser::parseLambda2(Configuration& config, MFSUBMODELS::Transformation::Configuration& trans_config, 
                     const LSFBase* param, ostream& info, cerrorstream& errors)
{
  if(! param->attrs())
  {
    errors << priority(warning) << "No attributes specified for parameter 'lambda2' of the transformation sub-block." << endl;
    return;
  }
                 
  AttrList::const_iterator attr_iter;
                       
  // - Fixity
  //  
  attr_iter = param->attrs()->find("FIXED");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(toUpper(attr_iter->second.String()) == "TRUE")
    {
      trans_config.set_lambda2_fixed(true);
    }
    else if(toUpper(attr_iter->second.String()) == "FALSE")
    {
      if(config.getTransformBothSides())
      {
        trans_config.set_lambda2_fixed(false);
      }
      else
      {
        errors << priority(error) << "lambda2 may not be estimated when the 'difference' " 
               << "is transformed.  Fixing lambda2 at 0 ..." << endl;
      }
    }
    else
    {   
      errors << priority(error) << "Value for 'fixed' attribute of 'lambda2' not understood.  Ignoring ..." << endl;
    }
  }

  // - Value
  //                           
  attr_iter = param->attrs()->find("VAL");
  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(trans_config.get_lambda2().initial_type != MAXFUN::Parameter::FIXED)
    {
      errors << priority(error) << "Initial value may not be specified for 'lambda2' if it is to be estimated."  
             << "  Ignoring ..." << endl;    
    }
    else if(finite(attr_iter->second.Real()))
    {
      trans_config.set_lambda2_init_est(attr_iter->second.Real());
    }
    else
    {   
      errors << priority(error) << "Value for 'val' attribute of 'lambda2' not understood.  Ignoring ..." << endl;
    }
  }
                                                                    
  // - Lower bound
  //
  attr_iter = param->attrs()->find("LOWER_BOUND");
  if(attr_iter != param->attrs()->end())
  {
    errors << priority(error) << "A lower bound may not be specified for 'lambda2'.  Ignoring ..." << endl;
  }
  
  // - Upper bound
  //
  attr_iter = param->attrs()->find("UPPER_BOUND");
  if(attr_iter != param->attrs()->end())
  {
    errors << priority(error) << "An upper bound may not be specified for 'lambda2'.  Ignoring ..." << endl;
  }
}


void
Parser::parseAllowAveraging(Configuration& config, const LSFBase* param, ostream& info, cerrorstream& errors)
{
  // 0. Set up local variables:

  string aa = "";
    
  // 1. Fetch & convert value:

  APP::ParsingFunctions::parse_string(param, aa);

  aa = toUpper(aa);

  // 2. Check aa option:

       if(aa == "MEAN" || aa == "") config.setAllowAveraging(true);
  else if(aa == "NONE")             config.setAllowAveraging(false);
  else                              info << "Error: User has specified '"
                                         << aa
                                         << "' for the allow-averaging options. Acceptable values are '', 'mean', or 'none'."
                                         << std::endl;
}


// - Added 6/6/7.  -djb
//
void
Parser::prohibitBinaryTransformation(Configuration& config, cerrorstream& errors)
{
  if(config.getDependentTraitType() == BINARY &&
     config.getTransConfig().get_type() != MFSUBMODELS::Transformation::NONE)
  {
    errors << priority(error) << "Transformation not allowed for binary dependent trait. "
                              << "Omitting transformation ..." << endl;
    config.setTransformationOption(MFSUBMODELS::Transformation::NONE);
  }
}

// - Trac #1764.  Summary display.
//
void 
Parser::parseSummaryDisplay(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  LSFList::const_iterator  p_iter     = param->List()->begin();
  LSFList::const_iterator  p_end_iter = param->List()->end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    if(*p_iter)
    {
      string param_name = toUpper((*p_iter)->name());

      if(param_name == "ORDER")              
      {
        parseOrder(config, *p_iter, errors);
      }
      else if(param_name == "FILTERS")
      {
        parseFilters(config, *p_iter, errors);
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' not recognized.  Skipping ..."  << endl;
      }
    }
  }
}

// - Trac #1764.  Summary display.
//
void  
Parser::parseOrder(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  string  order;
  
  APP::ParsingFunctions::parse_string(param, order);
  order = toUpper(order);
  
  if(order == "AS_INPUT")
  {
    config.setDisplayOrder(Configuration::AS_INPUT);
  }
  else if(order == "LRT")
  {
    config.setDisplayOrder(Configuration::BY_LRT);
  }
  else if(order == "WALD")
  {
    config.setDisplayOrder(Configuration::BY_WALD);
  }
  else if(order == "LARGER_PVALUE")
  {
    config.setDisplayOrder(Configuration::BY_LARGER_PVALUE);
  }
  else if(order == "PVALUE_RATIO")
  {
    config.setDisplayOrder(Configuration::BY_PVALUE_RATIO);
  }  
  else
  {
    errors << priority(error) << "Bad or missing value for 'order' parameter in 'summary_display' "
           << "sub-block.  Ignoring 'order' parameter ..." << endl;
  }
}
            
// - Trac #1764.  Summary display.
//
void  
Parser::parseFilters(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  // - Use two passes because "all" option overides the others.
  //
  LSFList::const_iterator p_iter     = param->List()->begin();
  LSFList::const_iterator p_end_iter = param->List()->end();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    if(*p_iter)
    {
      string  param_name = toUpper((*p_iter)->name());

      if(param_name == "ALL")              
      {
        parseAll(config, *p_iter, errors);
      }
      else if(! (param_name == "WALD" || param_name == "LRT" || param_name == "LIMIT_NUMBER"))
      {
        errors << priority(error) << "Parameter '" << param_name << "' in sub-block '" 
               << "filters' not recognized.  Skipping parameter '" << param_name 
               << "' ..."  << endl;      
      }
    }
  }
  
  p_iter = param->List()->begin();
  for(; p_iter != p_end_iter; ++p_iter)
  {
    if(*p_iter)
    {
      string  param_name = toUpper((*p_iter)->name());
      
      if(param_name == "WALD")
      {
        parseWald(config, *p_iter, errors);
      }
      else if(param_name == "LRT")
      {
        parseLRT(config, *p_iter, errors);
      }
      else if(param_name == "LIMIT_NUMBER")
      {
        parseLimitNumber(config, *p_iter, errors);
      }
    }  
  }
}

void
Parser::parseAll(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  bool  display_all = config.getDisplayAll();

  APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_boolean(param, display_all, false);
  if(status == APP::LSFConvert::GOOD)
  {
    config.setDisplayAll(display_all);
    if(display_all)                       // 'all' overrides all other filters.
    {
      config.setFilterByWald(false);
      config.setFilterByLRT(false);
      config.setLimitNumberDisplayed(false);
    }
  }
  else
  {
    errors << priority(error) << "Value for parameter 'all' in 'summary_display' sub-block "
           << "not understood.  Ignoring value for parameter 'all' ..." << endl;  
  }
}

void
Parser::parseWald(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  if(config.getDisplayAll() == false)
  {
    bool  filter_by_wald = config.getFilterByWald();

    APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_boolean(param, filter_by_wald, false);
    if(status == APP::LSFConvert::GOOD)
    {
      config.setFilterByWald(filter_by_wald);
    }
    else
    {
      errors << priority(error) << "Value for parameter 'wald' in 'summary_display' sub-block "
             << "not understood.  Ignoring value for parameter 'wald' ..." << endl;  
    }
  }
  else
  {
    errors << priority(error) << "Ignoring 'summary_display' sub-block parameter 'wald' because "
           << "parameter 'all' equals 'true' ..." << endl;
  }
  
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("max"))
    {
      double  maximum_pvalue = 0;
      APP::ParsingFunctions::parse_real(param, "MAX", maximum_pvalue, false);
      if(! (0.0 < maximum_pvalue && maximum_pvalue <= 1.0))
      {
        errors << priority(error) << "Value of attribute 'max' for parameter 'wald' in "
               << "'filters' sub-block is invalid.  Ignoring this value ..." << endl;
      }
      else
      {
        config.setWaldFilterThreshold(maximum_pvalue);
      }
    }
  }
}

void
Parser::parseLRT(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  if(config.getDisplayAll() == false)
  {
    bool  filter_by_LRT = config.getFilterByLRT();

    APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_boolean(param, filter_by_LRT, false);
    if(status == APP::LSFConvert::GOOD)
    {
      config.setFilterByLRT(filter_by_LRT);
    }
    else
    {
      errors << priority(error) << "Value for parameter 'lrt' in 'summary_display' sub-block "
             << "not understood.  Ignoring value for parameter 'lrt' ..." << endl;  
    }
  }
  else
  {
    errors << priority(error) << "Ignoring 'summary_display' sub-block parameter 'lrt' because "
           << "parameter 'all' equals 'true' ..." << endl;
  }

  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("max"))
    {
      double  maximum_pvalue = 0;
      APP::ParsingFunctions::parse_real(param, "MAX", maximum_pvalue, false);
      if(! (0.0 < maximum_pvalue && maximum_pvalue <= 1.0))
      {
        errors << priority(error) << "Value of attribute 'max' for parameter 'lrt' in "
               << "'filters' sub-block is invalid.  Ignoring this value ..." << endl;
      }
      else
      {
        config.setLRTFilterThreshold(maximum_pvalue);
      }
    }
  }  
}

void
Parser::parseLimitNumber(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  if(config.getDisplayAll() == false)
  {
    bool  limit_number_displayed = config.getLimitNumberDisplayed();

    APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_boolean(param, limit_number_displayed, false);
    if(status == APP::LSFConvert::GOOD)
    {
      config.setLimitNumberDisplayed(limit_number_displayed);
    }
    else
    {
      errors << priority(error) << "Value for parameter 'limit_number' in 'summary_display' sub-block "
             << "not understood.  Ignoring value for parameter 'limit_number' ..." << endl;  
    }
  }
  else
  {
    errors << priority(error) << "Ignoring 'summary_display' sub-block parameter 'limit_number' because "
           << "parameter 'all' equals 'true' ..." << endl;
  }
  
  AttrList*  a_list = param->attrs();
  if(a_list)
  {
    if(a_list->has_attr("number"))
    {
      int  number = config.getDisplayNumber();
      APP::LSFConvert::error_t  status = APP::ParsingFunctions::parse_integer(param, "NUMBER", number, false);
      if(status == APP::LSFConvert::GOOD)
      {
        if(number < 0)
        {
          errors << priority(error) << "Value of attribute 'number' for parameter 'limit_number' in "
                 << "'filters' sub-block is invalid.  Ignoring this value ..." << endl;
        }
        else
        {
          config.setDisplayNumber(number);
        }
      }
      else
      {
        errors << priority(error) << "Value for attribute 'number' of parameter 'limit_number' "
               << "in sub-block 'filters' not understood.  Ignoring this attribute ..." << endl;
      }
    }
  }
}


void
Parser::parseResiduals(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  LSFList::const_iterator      iter = param->List()->begin();
  LSFList::const_iterator  end_iter = param->List()->end();
  for(; iter != end_iter; ++iter)
  {
    string  param_name = toUpper((*iter)->name());
    if(param_name == "MODEL")
    {
      if((*iter)->attrs())
      {
        parseModel(config, *iter, errors);
      }
      else
      {
        errors << priority(error) << "No attributes given for parameter, model, in "
               << "residuals sub-block.  Ignoring parameter ..." << endl;
      }
    }
    else   
    {      
      errors << priority(error) << "Parameter '" << param_name 
             << "' in residuals sub-block not recognized.  Skipping ..."  << endl;
    }
  }              
}

void
Parser::parseModel(Configuration& config, const LSFBase* param, cerrorstream& errors)
{
  string  model_name = "";
  bool  null = false;
  bool  test = true;

  AttrList::const_iterator      attr_itr = param->attrs()->begin();
  AttrList::const_iterator  end_attr_itr = param->attrs()->end();
  for(; attr_itr != end_attr_itr; ++attr_itr)
  {
    string  attr_name = toUpper(AttrNameMgr.query(attr_itr->first));
    if(attr_name == "VALUE")
    {
      model_name = (attr_itr->second).String();
      if(toUpper(model_name) == "BASELINE")
      {
        model_name = "Baseline";
        null = true;
        test = false;
      }
      
      if(! config.hasModel(model_name))
      {
        errors << priority(error) << "Model name '" << model_name << "' given in residuals sub-block "
               << "not valid.  Skipping parameter, models ..." << endl;
        return;
      }
    }
    else if(attr_name == "NULL")
    {
      string  null_value = toUpper((attr_itr->second).String());
      if(null_value == "TRUE")
      {
        null = true;
      }
      else if(null_value == "FALSE")
      {
        null = false;
      }
      else
      {
        errors << priority(error) << "Attribute value '" << null_value << "' in residuals " 
               << "sub-block not understood.  Ignoring ..." << endl;
      }
    }
    else if(attr_name == "TEST")
    {
      string  test_value = toUpper((attr_itr->second).String());
      if(test_value == "TRUE")
      {
        test = true;
      }
      else if(test_value == "FALSE")
      {
        test = false;
      }
      else
      {
        errors << priority(error) << "Attribute value '" << test_value << "' in residuals "
               << "sub-block not understood.  Ignoring ..." << endl;
      }
    }
    else
    {
      errors << priority(error) << "Attribute '" << attr_name 
             << "' of parameter, model, in residuals sub-block not recognized.  Skipping ..."  << endl;
    }
  }
  
  Configuration::Model&  model = config.getModel(model_name);
  model.null_residuals = null;
  model.test_residuals = test;
}
            

} 
} 
