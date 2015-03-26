#include "mfsubmodels/TransformationParser.h"

using namespace std;

namespace SAGE   {
namespace MAXFUN {

//========================================================================
//
//  parse(...)
//
//=======================================================================
bool
TransformationParser::parse(
     TransformationSubmodel&           tsm,
     const LSFBase*                      param,
     TransformationSubmodel::sm_option option,
     cerrorstream&                       errors)
{
  static char exp_list[][13] = 
    { "VAL", "FIXED", "VALUE", "INTERACTION", "LOWER_BOUND", "UPPER_BOUND" };

  //bool  option_specified = false;                             
  
  model_input  lambda_one(QNAN, LAMBDA_ONE_DEFAULT_FIXED);
  model_input  lambda_two(QNAN, LAMBDA_TWO_DEFAULT_FIXED);
  double       lambda_one_lb = QNAN;
  double       lambda_one_ub = QNAN;

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
          if(option_string == "NONE")
          {
            option = TransformationSubmodel::no_trans;
            //option_specified = true;
          }
          else if(option_string == "BOX_COX")
          {
            option = TransformationSubmodel::box_cox;
            //option_specified = true;
          }
          else if(option_string == "GEORGE_ELSTON")
          {
            option = TransformationSubmodel::george_elston;
            //option_specified = true;
          }
          else
          {
            errors << priority(error) << "Value '" << option_string << "' for "
                   << "option parameter of transformation sub-block not recognized.  "
                   << "Ignoring ... " << endl;
          }
        }
      }
      else if(param_name == "LAMBDA1")
      {
        APP::check_for_incorrect_attributes(*iter, exp_list, 6, errors);
        
        if(! get_lambda_one(lambda_one, lambda_one_lb, lambda_one_ub, *iter, errors))
        {
          return false;
        }
        
      }
      else if(param_name == "LAMBDA2")
      {
        // Note: Checks only the first two attributes, because the others
        //       don't belong for lambda_two
        APP::check_for_incorrect_attributes(*iter, exp_list, 2, errors);

        if(! get_lambda_two(lambda_two, *iter, errors))
        {
          return false;
        }
      }
      else
      {
        errors << priority(error) << "Parameter '" << param_name << "' in transformation sub-block "
               << "not recognized.  Ignoring ..." << endl;
      }
    }
  }
  else
  {
    errors << priority(warning) << "transformation sub-block contains no parameters."
           << endl;
  }
  

  if(! tsm.set(option, lambda_one, lambda_two, lambda_one_lb, lambda_one_ub))
  {
    return false;
  }

  return true;

}

//======================================================================
//
//  get_value(...)
//
//======================================================================
bool 
TransformationParser::get_value(
     model_input&   mi,
     double         def,
     const LSFBase* param,
     const string&  name_phrase,
     cerrorstream&  errors)
{
  double value = QNAN;
  double val   = QNAN;

  APP::LSFConvert::error_t err1 = APP::LSFConvert::to_real(param, value);
  APP::LSFConvert::error_t err2 = APP::LSFConvert::to_real(param, "val", val);

  if(err1 == APP::LSFConvert::GOOD && err2 == APP::LSFConvert::GOOD)
  {
    errors << priority(error) << "Both value and val attribute set for "
           << name_phrase << ".  Will ignore value." << endl;
    err1 = APP::LSFConvert::INVALID;
  }

  switch(err1)
  {
    case APP::LSFConvert::GOOD          :
      mi.value = value;
      break;

    case APP::LSFConvert::ATTR_INVALID  :
      errors << priority(error) << "Value of " << name_phrase
             << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }

  switch(err2)
  {
    case APP::LSFConvert::GOOD          :
      mi.value = val;
      break;

    case APP::LSFConvert::ATTR_INVALID  :
      errors << priority(error) << "'val' of " << name_phrase
             << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }

  // If we haven't got a value by now, we can use the default.
  // If the default is nan, there's no real change.   
 
  if(SAGE::isnan(mi.value))
    mi.value = def;

  // Parse fixed attribute

  bool fixed;

  APP::LSFConvert::error_t fixed_err = APP::LSFConvert::to_boolean(param, "FIXED", fixed);

  switch(fixed_err)
  {
    case APP::LSFConvert::GOOD :
      mi.fixed = fixed;

      // - User specified value for fixed, but not for val.  Use default
      //   for val.  If left as QNAN, sub-model will ignore entire model_input.
      if(fixed && SAGE::isnan(mi.value))
      {
        errors << priority(critical) << "No value given for " << name_phrase
               << " with attribute, fixed, equal to 'true'.  Skipping analysis ..." << endl;
        return false;
      }
      break;

    case APP::LSFConvert::ATTR_INVALID :
      errors << priority(error) << "Value of 'fixed' attribute of " 
             << name_phrase << " not understood.  Ignoring ..." << endl;

    default :
      break;
  }
  
  return true;
}

//====================================================================================
//
//  get_lambda_one(...)
//
//====================================================================================
bool  
TransformationParser::get_lambda_one
  (model_input&   lambda_one,
   double&        lambda_one_lb, 
   double&        lambda_one_ub,
   const LSFBase* param,
   cerrorstream&  errors)
{
  string  name_phrase = "lambda_one in transformation sub-block";
  
  if(! get_value(lambda_one, QNAN, param, name_phrase, errors))
  {
    return false;
  }

  get_lambda_one_bounds(lambda_one_lb, lambda_one_ub, param, name_phrase, 
                        lambda_one.fixed, errors);
                        
  return true;
}


//==================================================================================
//
//  get_lambda_two(...)
//
//==================================================================================
bool  
TransformationParser::get_lambda_two
  (model_input&   lambda_two,
   const LSFBase* param,
   cerrorstream&  errors)
{
  string  name_phrase = "lambda_two in transformation sub-block";
  
  // - fixed default available for this parameter.
  //
  if(! get_value(lambda_two, LAMBDA_TWO_DEFAULT_VALUE, param, name_phrase, errors))
  {
    return false;
  }
  
  return true;
}

//==================================================================================
//
//  get_lambda_one_bounds(...)
//
//==================================================================================
void  
TransformationParser::get_lambda_one_bounds(
      double&        lambda_one_lb,
      double&        lambda_one_ub,
      const LSFBase* param,
      const string&  name_phrase,
      bool           fixed,
      cerrorstream&  errors)
{
  AttrList* a_list = param->attrs();
  AttrList::const_iterator iter;
  if(a_list)
  {
    // - Lower bound.
    //
    iter = a_list->find("LOWER_BOUND");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          lambda_one_lb = a.Real();
        }
        else
        {   
          errors << priority(error) << "Value for 'lower_bound' attribute of " << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }  
    }    
         
    // - Upper bound.
    //
    iter = a_list->find("UPPER_BOUND");
    if(iter != a_list->end())
    {
      AttrVal a  = iter->second;
      if(a.has_value())
      {
        if(finite(a.Real()))
        {
          lambda_one_ub = a.Real();
        }
        else
        {   
          errors << priority(error) << "Value for 'upper_bound' attribute of" << name_phrase
                 << " not understood.  Ignoring ..." << endl;
        }
      }  
    }    
  }      
}   

} // End namespace MFSUBMODELS
} // End namespace SAGE

