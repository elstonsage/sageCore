#include "mfsubmodels/Transformation.h"

namespace SAGE        {
namespace MFSUBMODELS {

/// Parses an LSFBase object and returns the Configuration object
///
/// \param param The LSFBase object to be parsed
/// \param info  Where to send information messages
/// \param error Where to send error messages
Transformation::Configuration 
Transformation::Parser::parse_block(const LSFBase * param, std::ostream & info, cerrorstream & errors)
{
  Configuration config;
        
  if(!param->List())
    return config;
          
  for(LSFList::const_iterator iter = param->List()->begin(); iter != param->List()->end(); ++iter)
  {
    std::string param_name = toUpper((*iter)->name());

         if(param_name == "OPTION")  parse_option  (config, *iter, info, errors);
    else if(param_name == "LAMBDA1") parse_lambda1 (config, *iter, info, errors);
    else if(param_name == "LAMBDA2") parse_lambda2 (config, *iter, info, errors);
  } 
    
  return config;
}

void Transformation::Parser::parse_option
    (Configuration& config,
     const LSFBase* param,
     std::ostream & info,
     cerrorstream & errors)
{
  if(!param->attrs())
  {
    errors << priority(warning) << "No option specified for option parameter of transformation sub-block." << std::endl;
    return;
  }

  AttrList::const_iterator attr_iter = param->attrs()->find("VALUE");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    std::string option_string = toUpper(attr_iter->second.String());

         if(option_string == "NONE")          config.set_type(NONE);
    else if(option_string == "BOX_COX")       config.set_type(BOX_COX);
    else if(option_string == "GEORGE_ELSTON") config.set_type(GEORGE_ELSTON);
    else
    {
      errors << priority(error) << "Value '" << option_string << "' for option parameter of transformation sub-block not recognized.  Ignoring ... " << std::endl;
    }
  }
}

void Transformation::Parser::parse_lambda1
    (Configuration& config,
     const LSFBase* param,
     std::ostream & info,
     cerrorstream & errors)
{
  if(!param->attrs())
  {
    errors << priority(warning) << "No attributes specified for lambda1 parameter of transformation sub-block." << std::endl;
    return;
  }

  AttrList::const_iterator attr_iter;
    
  // Val
    
  attr_iter = param->attrs()->find("VAL");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.set_lambda1_init_est(attr_iter->second.Real());
    }
    else
    {   
      errors << priority(error) << "Value for 'val' attribute of lambda1 not understood.  Ignoring ..." << endl;
    }
  }  
    
  // Fixed

  attr_iter = param->attrs()->find("FIXED");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(toUpper(attr_iter->second.String()) == "TRUE")
    {
      config.set_lambda1_fixed(true);
    }
    else if(toUpper(attr_iter->second.String()) == "FALSE")
    {
      config.set_lambda1_fixed(false);
    }
    else
    {   
      errors << priority(error) << "Value for 'fixed' attribute of lambda1 not understood.  Ignoring ..." << endl;
    }
  }  

  // Lower bound

  attr_iter = param->attrs()->find("LOWER_BOUND");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.my_lambda1.lower_bound = attr_iter->second.Real();
    }
    else
    {   
      errors << priority(error) << "Value for 'lower_bound' attribute of lambda1 not understood.  Ignoring ..." << endl;
    }
  }  

  // Upper bound

  attr_iter = param->attrs()->find("UPPER_BOUND");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.my_lambda1.upper_bound = attr_iter->second.Real();
    }
    else
    {   
      errors << priority(error) << "Value for 'upper_bound' attribute of lambda1 not understood.  Ignoring ..." << endl;
    }
  }
}
  
void Transformation::Parser::parse_lambda2
    (Configuration& config,
     const LSFBase* param,
     std::ostream & info,
     cerrorstream & errors)
{
  if(!param->attrs())
  {
    errors << priority(warning) << "No attributes specified for lambda2 parameter of transformation sub-block." << std::endl;
    return;
  }

  AttrList::const_iterator attr_iter;
    
  // Val
    
  attr_iter = param->attrs()->find("VAL");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.set_lambda2_init_est(attr_iter->second.Real());
    }
    else
    {   
      errors << priority(error) << "Value for 'val' attribute of lambda2 not understood.  Ignoring ..." << endl;
    }
  }  
    
  // Fixed

  attr_iter = param->attrs()->find("FIXED");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(toUpper(attr_iter->second.String()) == "TRUE")
    {
      config.set_lambda2_fixed(true);
    }
    else if(toUpper(attr_iter->second.String()) == "FALSE")
    {
      config.set_lambda2_fixed(false);
    }
    else
    {   
      errors << priority(error) << "Value for 'fixed' attribute of lambda2 not understood.  Ignoring ..." << endl;
    }
  }  

  // Lower bound

  attr_iter = param->attrs()->find("LOWER_BOUND");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.my_lambda2.lower_bound = attr_iter->second.Real();
    }
    else
    {   
      errors << priority(error) << "Value for 'lower_bound' attribute of lambda2 not understood.  Ignoring ..." << endl;
    }
  }  

  // Upper bound

  attr_iter = param->attrs()->find("UPPER_BOUND");

  if(attr_iter != param->attrs()->end() && attr_iter->second.has_value())
  {
    if(finite(attr_iter->second.Real()))
    {
      config.my_lambda2.upper_bound = attr_iter->second.Real();
    }
    else
    {   
      errors << priority(error) << "Value for 'upper_bound' attribute of lambda2 not understood.  Ignoring ..." << endl;
    }
  }
}

OUTPUT::Table 
Transformation::Configuration::summarize_as_table() const
{
  OUTPUT::Table t("Transformation configuration");

  if(get_type() == NONE)
  {
    t << OUTPUT::NamedString("Note", "No transformation applied.");
  }
  else
  {
    t << (OUTPUT::TableRow() << "Method" << (get_type() == BOX_COX ? "Box-Cox" : "George-Elston"));
        
    if(get_lambda1().initial_type == MAXFUN::Parameter::FIXED)
    {
      t << (OUTPUT::TableRow() << "Lambda1 is fixed at " << get_lambda1().initial_estimate);
    }
    else
    {
      t << (OUTPUT::TableRow() << "Lambda1 is estimated");
    }
      
    if(get_lambda2().initial_type == MAXFUN::Parameter::FIXED)
    {
      t << (OUTPUT::TableRow() << "Lambda2 is fixed at " << get_lambda2().initial_estimate);
    }
    else
    {
      t << (OUTPUT::TableRow() << "Lambda2 is estimated");
    }
  }
  
  return t;
}

OUTPUT::Section 
Transformation::Configuration::dump() const
{
  return (OUTPUT::Section("Transformation config") 
    << OUTPUT::NamedString("Mode", my_type == BOX_COX ? "Box-Cox" : (my_type == NONE ? "None" : "George-Elston")) 
    << get_lambda1().dump() 
    << get_lambda2().dump());
}

/// Calculates a geometric mean, if required by the model
///
/// \param traits The vector of values on the basis of which the geometric mean will be calculated.
/// \retval true Mean was successfully calculated or mean not required (none model)
/// \retval false Mean was \b not successfully calculated.
bool Transformation::Facade::calculate_geometric_mean(const std::vector<double> & traits, bool calculate)
{
  if(calculate)
  {
    return (this->*my_calc_fcn)(traits);
  }
  else
  {
    my_geom_mean = 1.0;
    return  true;
  }
}

double
Transformation::Facade::geometric_mean() const
{
  return  my_geom_mean;
}


///
/// Transforms a vector of values.
/// NOTE: This function only works if the transformation type is NOT 'NONE'.
///
/// \param traits The vector of values which will be transformed.
/// \retval true Transformation was successful.
/// \retval false Transformation was \b not successful.
bool Transformation::Facade::transform(std::vector<double>& traits) const
{
  for(std::vector<double>::iterator i = traits.begin(); i != traits.end(); ++i)
    if(!transform(*i))
      return false;
      
  return true;
}

///
/// Transforms a single value. 
/// NOTE: This function only works if the transformation type is NOT 'NONE'.
///
/// \param trait The value which will be transformed.
/// \retval true Transformation was successful.
/// \retval false Transformation was \b not successful.
bool Transformation::Facade::transform(double & trait) const
{
  return (this->*my_transf_fcn)(trait);
}

bool Transformation::Facade::calculate_geom_mean_none
      (const std::vector<double>& trait)
{
  return true;
}
bool Transformation::Facade::calculate_geom_mean_box_cox
      (const std::vector<double>& traits)
{
  double lambda2 = get_lambda2();

  // Reset the geometric mean to QNAN:
  my_geom_mean = QNAN;
  
  // Calculate the number of non-NaN attributes and the product of their
  // shifted values
  size_t     adj_trait_count =  0;
  log_double G                 (1.0);
  
  for(std::vector<double>::const_iterator iter = traits.begin(); iter != traits.end(); ++iter)
  {
    if(!SAGE::isnan(*iter))
    { 
      ++adj_trait_count;
  
      double shifted_t = *iter + lambda2;
  
      if(shifted_t <= 0)
      {
        // If any value is negative after shifting, the geometric mean can't 
        // be calculated
        return false;
      }
      else
      {
        G *= shifted_t;
      }
    }
  }
  
  // If there were elements, calculate the geometric mean
  if(adj_trait_count) 
  {
    G.pow(1.0 / adj_trait_count);
    my_geom_mean = G.get_double();
  }
  
  // Return whether or not the geometric mean is QNAN:
  return !SAGE::isnan(my_geom_mean);        
}

bool Transformation::Facade::calculate_geom_mean_george_elston
      (const std::vector<double>& traits)
{
  double lambda2 = get_lambda2();

  // Reset the geometric mean to QNAN:
  my_geom_mean = QNAN;
  
  // Calculate the number of non-NaN attributes and the product of their
  // shifted absolute values
  size_t     adj_trait_count = 0;
  log_double G(1.0);

  for(std::vector<double>::const_iterator iter = traits.begin(); iter != traits.end(); ++iter)
  {
    if(!SAGE::isnan(*iter))
    {
      ++adj_trait_count;

      G *= fabs(*iter + lambda2) + 1;
    }
  }

  if(adj_trait_count)
  {
    G.pow(1.0 / (double)adj_trait_count);
    my_geom_mean = G.get_double();
  }
  
  // Return whether or not the geometric mean is QNAN:
  return !SAGE::isnan(my_geom_mean);        
}

bool Transformation::Facade::transform_none          (double& trait) const
{
  return true;
}

bool Transformation::Facade::transform_box_cox       (double& trait) const
{
  if(SAGE::isnan(trait))
    return false;

  double shifted_trait = trait + get_lambda2();
    
  if(shifted_trait < 0.0) 
  {
    return false;
  }
  
  double lambda1 = get_lambda1();
  
  if(lambda1 < 1e-10)  // Changed from lambda1 == 0 by djb per RCE 6-15-10
  {
    trait = my_geom_mean * log(shifted_trait);
  }
  else // lambda1 != 0
  {
    trait = (pow(shifted_trait, lambda1) - 1.0) / (lambda1 * pow(my_geom_mean, lambda1 - 1.0));
  }

  return !SAGE::isnan(trait);
}

bool Transformation::Facade::transform_george_elston (double& trait) const
{
  if(SAGE::isnan(trait))
    return false;

  double shifted_trait = trait + get_lambda2();
    
  double lambda1 = get_lambda1();
  
#define JACOBIAN 
#if defined(JACOBIAN)
  if(abs(lambda1) < 1e-10)     // Changed from lambda1 == 0 by djb per RCE 6-15-10
  {
    trait = (shifted_trait < 0.0 ? -1.0 : (shifted_trait == 0 ? 0.0 : 1.0))  * my_geom_mean * log(abs(shifted_trait) + 1);
  }
  else                         // lambda1 != 0
  {
    trait = (shifted_trait < 0.0 ? -1.0 : (shifted_trait == 0 ? 0.0 : 1.0)) * (pow(abs(shifted_trait) + 1, lambda1) - 1) /
            (lambda1 * pow(my_geom_mean, lambda1 - 1.0));
  }
#else
  if(abs(lambda1) < 1e-10)     // Changed from lambda1 == 0 by djb per RCE 6-15-10
  {
    trait = (shifted_trait < 0.0 ? -1.0 : (shifted_trait == 0 ? 0.0 : 1.0))  * log(abs(shifted_trait) + 1);
  }
  else                         // lambda1 != 0
  {
    trait = (shifted_trait < 0.0 ? -1.0 : (shifted_trait == 0 ? 0.0 : 1.0)) * (pow(abs(shifted_trait) + 1, lambda1) - 1) /
            lambda1;
  }
#endif
  
  return !SAGE::isnan(trait);
}  

} // End namespace MFSUBMODELS
} // End namespace SAGE

