#include "freq/Parser.h"

namespace SAGE {
namespace FREQ {

//========================================
//
//  parseFreqBlock(...)
//
//========================================
Configuration 
Parser::parseFreqBlock(size_t analysis_num, const LSFBase * param_block, const FPED::FilteredMultipedigree & mp, cerrorstream & errors)
{
  // Create configuration:
  Configuration config;
  
  // Check the 'out' attribute for the freq block:

  AttrVal a = attr_value(param_block, "out");

  if(a.has_value() && (a.String() != ""))
  {
    config.ofilename = a.String();
  }
  else
  {
    std::ostringstream s;

    s << "freq";

    if(analysis_num)
      s << analysis_num;

    config.ofilename = s.str();
  }

  // Check all the freq block parameters:

  for(LSFList::const_iterator param = param_block->List()->begin(); param != param_block->List()->end(); ++param)
  {
    if(!*param)
      continue;

    if(!(*param)->name().size()) 
      continue;

    string name = toUpper((*param)->name());

    if(name == "MARKER")
    {
      AttrVal a = attr_value(*param, 0);

      if(a.has_value())
      {
        std::string marker_name = strip_ws(a.String());
        size_t      marker_id   = mp.info().marker_find(marker_name);

        if(marker_id < mp.info().marker_count())
        {
          config.addMarker(marker_id);
        }
        else
        {
          errors << priority(error) << "Marker '" << marker_name << "' not found." << endl;
        }
      }
      else
      {
        errors << priority(error) << "Marker parameter is missing name.  Skipping..." << endl;
      }
    }
    else if(name == "FOUNDER_WEIGHT" || name == "WEIGHT")
    {
      AttrVal a = attr_value(*param, 0);

      if(a.has_value())
      {
        if(toUpper(a.String()) == "NONE")
        {
          config.use_founder_weight = false;
        }
        else if(finite(a.Real()) && (a.Real() >= 0) && (a.Real() <= 1))
        {
          config.use_founder_weight = true;
          config.founder_weight     = a.Real();
        }
        else
        {
          errors << priority(error) << "Invalid value '" << a.Real() << "' for parameter 'founder_weight'" << endl;
        }
      }
      else
      {
        errors << priority(warning) << "Parameter 'founder_weight' requires a value. Parameter will be ignored." << endl;
      }
    }
    else if(name == "SKIP_MLE")
    {
      APP::ParsingFunctions::parse_boolean(*param, config.skip_mle);
    }
    else if(name == "INBREEDING")
    {
      APP::ParsingFunctions::parse_boolean(*param, config.estimate_inbreeding);
    }
    else if(name == "MAXFUN_DEBUG")
    {
      config.maxfun_debug = true;
    }
  }
  
  // If no markers were specified, add ALL of them to the list!
  if(config.getMarkers().size() == 0)
  {
    for(size_t i = 0; i < mp.info().marker_count(); ++i)
      config.addMarker(i);
  }
  
  return config;
}

//========================================
//
//  createDefaultConfig(...)
//
//========================================
Configuration 
Parser::createDefaultConfig(const FPED::FilteredMultipedigree & mp)
{
  Configuration config;

  for(size_t i = 0; i < mp.info().marker_count(); ++i)
    config.addMarker(i);

  config.skip_mle           = false;
  config.use_founder_weight = false;
  config.ofilename          = "freq";
  
  return config;
}
    

} // End namespace FREQ
} // End namespace SAGE

