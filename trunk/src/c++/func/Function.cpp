#include "func/Function.h"
#include <sstream>

namespace SAGE {
namespace FUNC {

//============================================================================
//
//  create_trait(...) #1
//
//============================================================================
void 
Function::create_trait(
  cerrorstream errors,
  RPED::MultiPedigree & mp,
  const std::string & trait_name,
  const std::string & function_code)
{
  std::istringstream code_str("trait=" + trait_name + ",expression=\"" + function_code + "\"");
  
  LSFBase * new_params = new LSFBase ("function");
  
  LSF_input(code_str).input_to(new_params, false);
    
  create_trait(errors, mp, FunctionParser::create_trait_data(new_params, errors));
}

//============================================================================
//
//  create_trait(...) #2
//
//============================================================================
void 
Function::create_trait(
  cerrorstream errors,  
  RPED::MultiPedigree & mp,
  const FunctionParser::TraitData & par_data)
{
  Function(errors).create(mp, par_data);
}

//============================================================================
//
// Constructor
//
//============================================================================
Function::Function(cerrorstream errors) : my_errors(errors) {}

//============================================================================
//
//  create(...)
//
//============================================================================
void
Function::create(RPED::RefMultiPedigree & mp, const FunctionParser::TraitData & par_data)
{
  // Check for empty name or expression.

  if(!par_data.trait_name.size() || !par_data.expr.size())
  {
    write_error_msg(missing_info_msg, "", "");
    return;
  }

  // Check that user supplied name for variable to be created is not that
  // of an existing variable.

  if(mp.info().trait_exists(par_data.trait_name))
  {
    write_error_msg(existing_name_msg, par_data.trait_name, "");
    return;
  }
  
  // Create python interface and do a few more tests:
  PythonInterface calculator(my_errors);

  if(!set_constants(calculator, par_data.constants, par_data.time_limit))
  {
    return;
  }
  
  if(!find_variables_used(mp, par_data.expr, calculator))
  {
    write_error_msg(bad_syntax_msg, "", par_data.expr);
    return;
  }
  
  if(!calculator.compile(par_data.expr))
  {
    write_error_msg(bad_syntax_msg, "", par_data.expr);
    return;
  }
  
  // Create the trait in the multipedigree:

  size_t trait_num = 0;

  if(par_data.binary)
  {
    trait_num = mp.info().add_binary_trait(par_data.trait_name, par_data.usage);

    mp.info().trait_info(trait_num).set_string_affected_code             (par_data.affected);        
    mp.info().trait_info(trait_num).set_numeric_affected_code   (str2doub(par_data.affected));
    mp.info().trait_info(trait_num).set_string_unaffected_code           (par_data.unaffected);        
    mp.info().trait_info(trait_num).set_numeric_unaffected_code (str2doub(par_data.unaffected));
  }
  else // Continuous trait
  {
    trait_num = mp.info().add_continuous_trait(par_data.trait_name, par_data.usage);
  }

  mp.info().trait_info(trait_num).set_string_missing_code           (par_data.missing);        
  mp.info().trait_info(trait_num).set_numeric_missing_code (str2doub(par_data.missing));

  // Resize the vectors in each RefPedInfo object to accomodate the new trait:
  for(RPED::PedigreeIterator p_iter = mp.pedigree_begin(); p_iter != mp.pedigree_end(); ++p_iter)
    p_iter->info().resize_traits(p_iter->info().trait_count() + 1);
    
  // Calculate and set values for the new trait:

  bool runtime_error_written = false;

  for(size_t i = 0; i < mp.member_count(); ++i)
  {
    RPED::Member & member      = mp.member_index(i);
    double         trait_value = numeric_limits<double>::quiet_NaN();
      
    if(setMemberSpecificVariables(calculator, member))
    {
      trait_value = calculator.run(par_data.expr, par_data.time_limit);

      if(SAGE::isnan(trait_value) && (!runtime_error_written))
      {
        my_errors << priority(error) 
                  << "Could not evaluate function block expression '" << par_data.expr 
                  << "' for one or more individuals.  Missing value(s) generated." << endl;

        runtime_error_written = true;
      }
    }
    
    // Set the trait value:
    member.pedigree()->info().set_trait(
      member.index(), 
      trait_num, 
      SAGE::isnan(trait_value) ? mp.info().trait_info(trait_num).string_missing_code() : doub2str(trait_value), 
      mp.info().trait_info(trait_num));
  }
}

//============================================================================
//
//  find_variables_used(...)
//
//============================================================================
bool                
Function::find_variables_used(
    const RPED::RefMultiPedigree & mp, 
    const std::string            & expr, 
          PythonInterface        & calculator)
{
  std::vector<std::string> expr_names;

  try
  {
    expr_names = calculator.getNameList(expr);
  }
  catch(const PythonInterface::PythonException & e)
  {
    return false;
  }
  
  my_traits_used  . clear();
  my_strings_used . clear();
  my_markers_used . clear();

  for(std::vector<std::string>::iterator name_itr = expr_names.begin(); name_itr != expr_names.end(); ++name_itr)
  {
    if(name_itr->size() == 0)
      continue;

    if((*name_itr)[0] == '_') // Names begining w. an underscore are reserved for internal use.
      return false;

    if(valid_trait(mp, *name_itr))
      my_traits_used.insert(make_pair(*name_itr, mp.info().trait_find(*name_itr)));

    // If the string is in the multipedigree, add it:
    if(mp.info().string_find(*name_itr) != (size_t)-1)
      my_strings_used.insert(make_pair(*name_itr, mp.info().string_find(*name_itr)));
          
    // If the marker is in the multipedigree and is codominant, add it:
    size_t marker_id = mp.info().marker_find(*name_itr);
    
    if(marker_id != (size_t)-1 && mp.info().marker_info(marker_id).codominant())
      my_markers_used.insert(make_pair(*name_itr, marker_id));
  }

  return no_duplicates();
}

//============================================================================
//
//  valid_trait
//
//============================================================================
bool
Function::valid_trait(const RPED::RefMultiPedigree& mp, const string& name) const
{
  size_t trait_index = mp.info().trait_find(name);
  
  if(trait_index != (size_t)(-1))
  {
    RPED::RefTraitInfo::trait_t type = mp.info().trait_info(trait_index).type();
    
    return (type == RPED::RefTraitInfo::continuous_trait ||
            type == RPED::RefTraitInfo::binary_trait     ||
            type == RPED::RefTraitInfo::discrete_trait   ||
            type == RPED::RefTraitInfo::categorical_trait);
  }
  else
  {
    return false;
  }
}

// - Are any names common to my_traits_used and my_markers_used?
//
bool
Function::no_duplicates() const
{
  typedef std::list<const input_variable_map*> variable_map_list;
  variable_map_list maps;

  maps.push_back(&my_traits_used);
  maps.push_back(&my_strings_used);
  maps.push_back(&my_markers_used);

  variable_map_list::const_iterator   i,j;
  input_variable_map::const_iterator  i_iter, j_iter;

  for(i = maps.begin(); i != maps.end(); ++i)
    for(j = i,++j; j != maps.end(); ++j)
      for(i_iter = (*i)->begin(); i_iter != (*i)->end(); ++i_iter)
        for(j_iter = (*j)->begin(); j_iter != (*j)->end(); ++j_iter)
          if(toUpper(i_iter->first) == toUpper(j_iter->first))
            return false;

  return true;
}

bool
Function::set_constants(PythonInterface& calculator, 
  const FunctionParser::TraitData::constant_vector& constants,
  unsigned int time_limit) 
{
  for(FunctionParser::TraitData::constant_vector::const_iterator iter = constants.begin(); iter != constants.end(); ++iter)
  {
    if(!calculator.calculateValueInEnvironment(iter->first, iter->second, time_limit))
    {
      write_error_msg(bad_const_expr_msg, "", iter->second);
      return false;
    }
  }
  
  return true;
}

//=======================================
//
//  setMemberSpecificVariables(...)
//
//=======================================
bool 
Function::setMemberSpecificVariables(const PythonInterface & interface, const RPED::Member & member) const
{
  return setMemberSpecificTraits  (interface, member) && 
         setMemberSpecificStrings (interface, member) && 
         setMemberSpecificMarkers (interface, member);
}

//=======================================
//
//  setMemberSpecificTraits(...)
//
//=======================================
bool                
Function::setMemberSpecificTraits(const PythonInterface& calculator, const RPED::Member & member) const
{
  for(input_variable_map::const_iterator iter = my_traits_used.begin(); iter != my_traits_used.end(); ++iter)
  {
    double value = member.pedigree()->info().trait(member.index(), iter->second);

    if(SAGE::isnan(value))
      return false;
    else 
      calculator.addDoubleToEnvironment(iter->first, value);
  }

  return true;
}

//=========================================
//
//  setMemberSpecificStrings(...)
//
//=========================================
bool                
Function::setMemberSpecificStrings(const PythonInterface& calculator, const RPED::Member & member) const
{
  for(input_variable_map::const_iterator iter = my_strings_used.begin(); iter != my_strings_used.end(); ++iter)
    if(!calculator.addStringToEnvironment(iter->first, member.pedigree()->info().get_string(member.index(), iter->second)))
      return false;

  return true;
}
                       
//===========================================
//
//  setMemberSpecificMarkers(...)
//
//===========================================
bool                
Function::setMemberSpecificMarkers(const PythonInterface & calculator, const RPED::Member & member) const
{
  for(input_variable_map::const_iterator m_iter = my_markers_used.begin(); m_iter != my_markers_used.end(); ++m_iter)
  {
    size_t marker_id   = m_iter->second;
    string marker_name = m_iter->first;

    if(marker_id == (size_t)(-1) || marker_id >= member.multipedigree()->info().marker_count())
      return false;
    
    RPED::RefMarkerInfo mkri = member.multipedigree()->info().marker_info(marker_id);

    if(member.pedigree()->info().phenotype_missing(member.index(), marker_id, mkri))
      return false;
    
    size_t pheno_id = member.pedigree()->info().phenotype(member.index(), marker_id);

    if(pheno_id == MLOCUS::NPOS)
      return false;
    
    // Find the alleles corresponding to the phenotype.  Only codominant markers
    // were put into my_markers_used list.

    MLOCUS::penetrance_model::unphased_penetrance_iterator up_iter = mkri.unphased_penetrance_begin(pheno_id);

    std::string allele_one = up_iter.unphased_geno().allele1().name(),
                allele_two = up_iter.unphased_geno().allele2().name(); 
    
    if(!calculator.addAlleleListToEnvironment(marker_name, allele_one, allele_two))
      return false;
  }
  
  return true;
}

// - Undo a partial trait addition.
//
void
Function::rollback(RPED::RefMultiPedigree& mp, RPED::RefMultiPedigree::pedigree_iterator& l_iter)
{
  RPED::RefMultiPedigree::pedigree_iterator stop = l_iter;
  ++stop;
  
  for(RPED::PedigreeIterator iter = mp.pedigree_begin(); iter != stop; ++iter)
    iter->info().resize_traits(iter->info().trait_count() - 1);
  
  mp.info().remove_last_trait();
}

void
Function::write_error_msg(Function::error_msg msg, const string& name, const string& expr)
{
  switch(msg)
  {
    case existing_name_msg: 

      my_errors << priority(error) 
                << "Name for new variable '" 
                << name
                << "', in function block is that of an existing variable, "
                << "Function parameter ignored." 
                << endl;
      break;
      
    case bad_syntax_msg: 

      my_errors << priority(error) 
                << "Function block expression, '" 
                << expr 
                << "', contains a syntax error, an ambiguous name "
                << "(used for both a trait and a marker), or a name "
                << "beginning with an underscore.  Function parameter ignored."
                << endl;
      break;  
      
    case eval_failure_msg: 

      my_errors << priority(error) 
                << "Could not evaluate function block "
                << "expression '" 
                << expr 
                << "' for " + name + "."
                << "  No value generated."
                << endl;
      break;  
      
    case internal_error_msg:      

      my_errors << priority(error)
                << "Internal error occured during "
                << "processing of function block expression, '" 
                << expr 
                << "'.  Variable creation aborted." 
                << endl;
      break;
      
    case bad_const_expr_msg:

      my_errors << priority(error) 
                << "Function block constant expression, '" 
                << expr
                << "', not understood.  "
                << "Function parameter ignored."
                << endl;
                
      break;
      
    case missing_info_msg:

      my_errors << priority(error)
                << "Function block parameter is missing a "
                << "variable name and/or expression. "
                << "Function parameter ignored."
                << endl;
      break;
      
    default:
      ;
  }
}

void
Function::write_error_msg(const string& message)
{
  my_errors << message << endl;
}

} // End namespace FUNC
} // End namespace SAGE
