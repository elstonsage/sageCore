#include "func/evalfunc.h"
#include "util/RegexUtils.h"
#include "func/ParentOfOrigin.h"

namespace SAGE {
namespace FUNC {

//================================================================
//
//  eval_func(...)
//
//================================================================
void
FunctionEvaluator::evaluateFunction(RPED::RefMultiPedigree & mp, cerrorstream errors, const LSFBase* params)
{
  // Loop across all parameters, looking for "FUNCTION":

  for(LSFList::const_iterator iter = params->List()->begin(); iter != params->List()->end(); ++iter)
  {
    // Skip this parameter if it is invalid, or not a "FUNCTION" subblock:
    if(!*iter || (toUpper((*iter)->name()) != "FUNCTION"))
      continue;

    // Set up list of variable to which this function will apply:
    std::vector<std::string> names(0);

    // Check to see if this function has the LIST attribute: (see trac ticket #1492, columbus sailed the ocean blue):
    if((*iter)->attrs() && (*iter)->attrs()->has_attr("LIST") && attr_value(*iter, "LIST").has_value())
    {
      std::string val = toUpper(attr_value(*iter, "LIST").String());
    
      if(val == "MARKERS")
      {
        for(size_t i = 0; i < mp.info().marker_count(); ++i)
          names.push_back(mp.info().marker_info(i).name());
      }
      else if(val == "TRAITS")
      {
        for(size_t i = 0; i < mp.info().trait_count(); ++i)
          names.push_back(mp.info().trait_info(i).name());
      }
      else
      {
        errors << priority(error) << "Unknown value '" << val << "' for attribute LIST. Ignoring function block..." << std::endl;
        return;
      }
    }

    // Parse the function into a Functionparser::TraitData:
    FunctionParser::TraitData trait_data = FunctionParser::create_trait_data(*iter, errors);

    if(names.size() == 0) // No list attribute found; evaluate normally:
    {
      evaluateSingleFunction(mp, errors, trait_data);
    }
    else // We've got a list of substitution names to deal with!
    {
      // Loop across all names; substitute then and invoke the function evaluator:
      for(std::vector<std::string>::const_iterator name_itr = names.begin(); name_itr != names.end(); ++name_itr)
      {
        // Copy over a temporary TraitData:
        FunctionParser::TraitData temp = trait_data;

        // Append the current name to the trait_name:
        temp.trait_name += "_" + *name_itr;        

        // Replace $name$ with the actual name:
        temp.expr = UTIL::RegexUtils::searchAndReplace(temp.expr, "\\$name\\$", *name_itr);

        // Replace $name$ with the acutal name in any constant expressions:
        for(size_t i = 0; i < temp.constants.size(); ++i)
          temp.constants[i].second = UTIL::RegexUtils::searchAndReplace(temp.constants[i].second, "\\$name\\$", *name_itr);

        // Evaluate the function:
        evaluateSingleFunction(mp, errors, temp);
      }
    }

  } // End parameter loop
}

//==========================================
//
//  evaluateSingleFunction(...)
//
//==========================================
void 
FunctionEvaluator::evaluateSingleFunction(
  RPED::MultiPedigree & mp,
  cerrorstream errors,
  const FunctionParser::TraitData & trait_data)
{
  // Skip this function if the parser says so:
  if(trait_data.skip)
    return;

  // Otherwise, does it match a pre-determined thing?
  if(UTIL::RegexUtils::doesPatternMatch("^mean_adj.*", trait_data.expr) ||
     UTIL::RegexUtils::doesPatternMatch("^var_adj.*",  trait_data.expr) ||
     UTIL::RegexUtils::doesPatternMatch("^z_score.*",  trait_data.expr) || 
     UTIL::RegexUtils::doesPatternMatch("^merge.*",    trait_data.expr))
  {
    AdjustedTraitCreator::createAdjustedTrait(mp, errors, trait_data);
  }

  // How about the 'TAI' thingy?
  else if(UTIL::RegexUtils::doesPatternMatch("^u?tai.*", trait_data.expr))
  {
    transmitted_allele_indicator(errors).createTaiStatusTrait(mp, trait_data);
  }

  // Or 'trim' or 'winsor' ?
  else if(UTIL::RegexUtils::doesPatternMatch("^trim.*",   trait_data.expr) ||
          UTIL::RegexUtils::doesPatternMatch("^winsor.*", trait_data.expr))
  {
    trim_winsor_process(errors).create_twp(mp, trait_data);
  }

  else if(UTIL::RegexUtils::doesPatternMatch("^poo.*",   trait_data.expr))
  {
    ParentOfOriginCalculator::createPooStatusTrait(mp, trait_data);
  }
  
  // I guess it's a custom trait...
  else
  {
    Function::create_trait(errors, mp, trait_data);
  }
}
              
} // End namespace FUNC
} // End namespace SAGE
