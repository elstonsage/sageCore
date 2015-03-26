#include "func/FunctionParser.h"

namespace SAGE {
namespace FUNC {

//===================================================
//  create_trait_data()
//===================================================
FunctionParser::TraitData 
FunctionParser::create_trait_data(const LSFBase* function_block, cerrorstream errors)
{
  FunctionParser  f;

  f.parse(function_block, errors);

  return f.get_data();
}

//===================================================
//  parse()
//===================================================
void
FunctionParser::parse(const LSFBase* function, cerrorstream errors)
{
  // Make sure that the LSFBase object is usable:
  assert(function != NULL);
  assert(function->List() != NULL);
  
  for(LSFList::const_iterator iter = function->List()->begin(); iter != function->List()->end(); ++iter)
  {
    if(!*iter)
      continue;
    
    string sub_name = toUpper((*iter)->name());
    
    if(sub_name == "CONSTANT")
    {
      parse_constant(iter, errors);
    }
    else if(sub_name == "TIME_LIMIT")
    {
      parse_time_limit(iter, errors);
    }    
    else if(sub_name == "TRAIT" || sub_name == "PHENOTYPE" || sub_name == "COVARIATE") 
    {
      // Has a valid trait, covariate, or phenotype
      // parameter already been found?
      if(my_data.trait_name.size() || my_data.expr.size())
      {
        write_error_msg(duplicate_trait, errors);
      }
      else 
      {
        parse_trait(iter, errors);
      }
    }
    else if(sub_name == "VERBOSE")
    {
      my_data.verbose = true;
    }
    else // Invalid parameter!
    {
      write_error_msg(invalid_param, errors, (*iter)->name());
    }
  }
}

//====================================
//  parse_constant()
//====================================
void
FunctionParser::parse_constant(const LSFList::const_iterator & iter, cerrorstream errors)
{
  // Make sure that CONSTANT has an assigned value and that the value is non-empty:
  if(attr_value(*iter, "CONSTANT", 0).has_value() && attr_value(*iter, "CONSTANT", 0).String() != "")
  {
    // Store the CONSTANT's value as 'name':
    string name = attr_value(*iter, "CONSTANT", 0).String();

    // Make sure that EXPRESSION attribute exists and has a non-empty value:
    if(attr_value(*iter, "CONSTANT", "EXPRESSION").has_value() && attr_value(*iter, "CONSTANT", "EXPRESSION").String() != "")
    {
      // Store the EXPRESSION's value as 'expr':
      string expr = attr_value(*iter, "CONSTANT", "EXPRESSION").String();
      
      // Make sure that this constant hasn't been added yet:
      bool name_present = false;

      for(TraitData::constant_vector::const_iterator pos = my_data.constants.begin(); pos != my_data.constants.end(); ++pos)
        name_present |= (pos->first == name);
        
      if(name_present) // Oops! name already added once...
      {
        write_error_msg(constant_redef, errors, name);
      }
      else // Let's add the name & expression:
      {
        my_data.constants.push_back(make_pair(name, expr));
      }
    }
    else // Oops! EXPRESSION is missing or has an empty value:
    {
      write_error_msg(no_constant_expr, errors);
    }
  }
  else // Oops! CONSTANT is missing or has an empty value:
  {
    write_error_msg(no_constant_name, errors);
  }
}

//============================================
//  parse_trait()
//============================================
void
FunctionParser::parse_trait(const LSFList::const_iterator& iter, cerrorstream errors)
{
  // The name of the trait type (TRAIT, PHENOTYPE, COVARIATE)
  string trait_type_name = toUpper((*iter)->name());
  
  if(trait_type_name == "TRAIT")
  {
    my_data.usage = RPED::RefTraitInfo::trait_variate;
  }
  else if(trait_type_name == "PHENOTYPE")
  {
    my_data.usage = RPED::RefTraitInfo::unknown_use;
  }
  else if(trait_type_name == "COVARIATE")
  {
    my_data.usage = RPED::RefTraitInfo::trait_covariate;
  }

  // First, make sure that there's a value associated with the trait_type_name and it's non-empty:
  if(attr_value(*iter, trait_type_name, 0).has_value() && attr_value(*iter, trait_type_name, 0).String() != "")
  {
    // Set the name that we've found:    
    my_data.trait_name = attr_value(*iter, trait_type_name, 0).String();

    // Now, make sure there's an expression and it's non-empty:
    if(attr_value(*iter, trait_type_name, "EXPRESSION").has_value() && attr_value(*iter, trait_type_name, "EXPRESSION").String() != "")
    {
      my_data.expr = attr_value(*iter, trait_type_name, "EXPRESSION").String();
    }
    else    // No expression.
    {
      write_error_msg(no_trait_expr, errors, trait_type_name);
      return;
    }
  }
  else     // No trait name.
  {
    write_error_msg(no_trait_name, errors, trait_type_name);
    return;
  }

  // Get the missing code:
  my_data.missing = (*iter)->attrs()->find("missing") == (*iter)->attrs()->end() ? "" : (*iter)->attrs()->find("missing")->second.String();

  // Check for binary trait specification:
  if((*iter)->attrs()->has_attr("binary"))
  {
    // Get affected & unaffected codes:
    string affected   = (*iter)->attrs()->find("affected")   == (*iter)->attrs()->end() ? "1" : (*iter)->attrs()->find("affected")   ->second.String(),
                unaffected = (*iter)->attrs()->find("unaffected") == (*iter)->attrs()->end() ? "0" : (*iter)->attrs()->find("unaffected") ->second.String();

    // Oops! Unaffected = affected!
    if(affected == unaffected)
    {
      errors << priority(error) << "Affected and unaffected codes in function block may not be the same.  Function skipped." << endl;

      my_data.skip = true;
      return;
    }
    
    // Oops! Affected = missing or Unaffected = missing:
    else if(affected == my_data.missing || unaffected == my_data.missing)
    {
      errors << priority(error) << "Affected and unaffected codes "
             << "in function block must not be the same as the missing value code.  " 
             << "Function skipped." << endl;

      my_data.skip = true;
      return;
    }

    // Everything's ok:
    else
    {
      my_data.binary     = true;
      my_data.affected   = affected;
      my_data.unaffected = unaffected;
    }
  }
} 

//==================================
//  parse_time_limit()
//==================================
void
FunctionParser::parse_time_limit(const LSFList::const_iterator& iter,
                                  cerrorstream errors                 )
{
  int  time_limit;
  AttrVal  value = attr_value(*iter, "TIME_LIMIT", 0);
  if(value.has_value())
  {
    time_limit = value.Int();
    if(time_limit < 1)
    {
      write_error_msg(bad_time_limit, errors);
      my_data.time_limit = DEFAULT_TIME_LIMIT;
    }
    else
    {
      my_data.time_limit = time_limit;
    }
  }
  else
  {
    my_data.time_limit = DEFAULT_TIME_LIMIT;
  }
}

void  
FunctionParser::write_error_msg(error_msg msg, cerrorstream errors, 
                                 string string_var ) const
{
  switch(msg)
  {
    case duplicate_trait:
      errors << priority(error) << "More than one function block parameter from the group: "
             << "trait, covariate, or phenotype was found."
             << "\nDuplicate parameter ignored." << endl;
      break;    
    case invalid_param:
      errors << priority(error) << "Invalid parameter in function block, '"
             << string_var << "'.\nParameter ignored." << endl;
      break;
    case no_constant_name:
      errors << priority(error) << "No name given for function block parameter, 'constant'."
             << "\nParameter ignored." << endl;
      break;
    case no_constant_expr:
      errors << priority(error) << "No expression given for function block parameter, 'constant'."
             << "\nParameter ignored." << endl;
      break;
    case no_trait_name:
      errors << priority(error) << "No name given for function block parameter, '" << string_var 
             << "'."
             << "\nParameter ignored." << endl;
      break;
    case no_trait_expr:
      errors << priority(error) << "No expression given for function block parameter, '" << string_var 
             << "'."
             << "\nParameter ignored." << endl;
      break;
    case bad_time_limit:
      errors << priority(warning) << "Invalid value for function block parameter, 'time limit'. "
             << "\nDefault value of " << DEFAULT_TIME_LIMIT << " used." << endl;
      break;
    case constant_redef:
      errors << priority(error) << "Function block constant, '" << string_var << "', already "
             << "defined."
             << "\nSecond definition ignored." << endl;
    default:
      ;
  }
}

} // End namespace FUNC
} // End namespace SAGE
