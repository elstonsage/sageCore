#include <string>
#include "LSF/LSF.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "rped/rped.h"
#include "tdtex/Parser.h"

namespace SAGE  {
namespace TDTEX {

//================================================================
//
//  parse_parameters(...)
//
//================================================================
Configuration 
Parser::parse_parameters(const RPED::RefMPedInfo & mped_info, const LSFBase * params, cerrorstream & err)
{
  Parser p(mped_info, err);
  
  p.parse_test_parameter_section(params);
  
  return p.get_configuration();
}
  
//======================================================
//
//  CONSTRUCTOR
//
//======================================================
Parser::Parser(const RPED::RefMPedInfo& mped_info, cerrorstream& err) : 
  APP::BasicParser  (err), 
  my_mped_info      (mped_info)
{ }
            

//==================================================================
//
//  parse_test_parameter_section(...)
//
//==================================================================
void
Parser::parse_test_parameter_section(const LSFBase* params)
{
  // Make sure the params are ok:
  if(!params)
    throw std::exception();

  if(!params->List())
    throw std::exception();

  // Grab the 'out' attribute (if it was set):
  AttrVal a = attr_value(params, "out");

  // Set the output filename:
  my_config.set_ofilename(a.has_value() && a.String().size() ? a.String() : "tdtex.out");

  // Process each parameter:
  for(LSFList::const_iterator i = params->List()->begin(); i != params->List()->end(); ++i)
    parse_test_parameter(*i);
}

//==================================================================
//
//  parse_test_parameter(...)
//
//==================================================================
void Parser::parse_test_parameter(const LSFBase* param)
{
  if(!param || !param->name().size()) return;

  AttrVal     a;
  std::string name = toUpper(param->name());

  if(name == "MARKER")
  {
    parse_marker(param);
  }
  else if(name == "TRAIT")
  {
    parse_trait(param);
  }
  else if(name == "PARENT_TRAIT" || name == "PARENTAL_TRAIT")
  {
    parse_parent_trait(param);
  }
  else if(name == "MAX_CHILDREN" || name == "MAX_CHILD")
  {
    size_t children;

    parse_limit(param, children, "max_children");

    my_config.set_max_children(children);
  }
  else if(name == "MAX_PAIRS" || name == "MAX_SIB_PAIRS" || name == "MAX_SIBS")
  {
    size_t pairs;
    parse_limit(param, pairs, "max_sib_pairs");
    my_config.set_max_sib_pairs(pairs);
  }
  else if(name == "SAMPLE")
  {
    parse_sample(param);
  }
  else if(name == "SKIP_EXACT" || name == "SKIP_EXACT_TEST" || name == "SKIP_EXACT_TESTS")
  {
    bool skip = my_config.get_skip_mc_test() && my_config.get_skip_mcmh_test() && my_config.get_skip_permutation_test();

    parse_boolean(param, skip);

    my_config.set_skip_mc_test          (skip);
    my_config.set_skip_mcmh_test        (skip);
    my_config.set_skip_permutation_test (skip);

  }
  else if(name == "SEX_DIFF" || name == "SEX_DIFFERENTIAL")
  {
    bool diff = my_config.get_sex_differential();

    parse_boolean(param, diff);

    my_config.set_sex_differential(diff);
  }
  else if(name == "SKIP_PERM" || name == "SKIP_PERMUTATION_TEST")
  {
    bool skip = my_config.get_skip_permutation_test();

    parse_boolean(param, skip);

    my_config.set_skip_permutation_test(skip);
  }
  else if(name == "SKIP_MC" || name == "SKIP_MC_TEST")
  {
    bool skip = my_config.get_skip_mc_test();

    parse_boolean(param, skip);

    my_config.set_skip_mc_test(skip);
  }
  else if(name == "SKIP_MCMH" || name == "SKIP_MCMH_TEST")
  {
    bool skip = my_config.get_skip_mcmh_test();

    parse_boolean(param, skip);

    my_config.set_skip_mcmh_test(skip);
  }
}

//================================================================
//
//  parse_limit(...)
//
//================================================================
void
Parser::parse_limit(
  const LSFBase * param, 
        attr_id   id,
        size_t  & value, 
  const string  & name)
{
  if(!param)
    return;

  AttrVal a = attr_value(param, id);

  if(a.has_value())
  {
    std::string n = toUpper(a.String());

    if(finite(a.Real()))
    {
      value = (size_t)a.Real();
    }
    else if(n == "UNLIMITED" || n == "ALL" || n == "TRUE")
    {
      value = (size_t)-1;
    }
    else if(n == "NONE" || n == "NO" || n == "FALSE")
    {
      value = 0;
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter '" << name << "'" << endl;
    }
  }
  else
  {
    errors << priority(warning) << "Parameter '" << name << "' requires a value.  Parameter will be ignored." << endl;
  }
}

inline void
Parser::parse_limit(const LSFBase* param, size_t& value, const std::string& name)
{
  parse_limit(param, 0, value, name);
}
   
inline void
Parser::parse_limit(const LSFBase* param, const std::string& attr, size_t& value, const std::string& name)
{
  parse_limit(param, AttrNameMgr.query(attr), value, name);
}
     

//==================================================
//
//  parse_marker(...)
//
//==================================================
void
Parser::parse_marker(const LSFBase* param)
{
  AttrVal a = attr_value(param, 0);

  if(a.has_value())
  {
    std::string mname = strip_ws(a.String());
    size_t      m     = my_mped_info.marker_find(mname);

    if(m == (size_t)-1)
    {
      errors << priority(error) << "Marker '" << mname << "' not found." << endl;
    }
    else
    {
      my_config.set_marker(m);
    }
  }
  else
  {
    errors << priority(error) << "Marker parameter is missing name.  Skipping..." << endl;
  }
}

//=====================================================
//
//  parse_trait(...)
//
//=====================================================
void
Parser::parse_trait(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if(a.has_value())
  {
    std::string tname = strip_ws(a.String());
    size_t      t     = my_mped_info.trait_find(tname);

    if(t == (size_t)-1)
    {
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
    }
    else
    {
      my_config.set_trait(t);
    }
  }
  else
  {
    errors << priority(error) << "Trait parameter is missing name.  Skipping..."  << endl;
  }
}

//==================================================
//
// parse_parent_trait(...)
//
//==================================================
void
Parser::parse_parent_trait(const LSFBase* param)
{
  AttrVal a = attr_value(param,0);

  if(a.has_value())
  {
    std::string tname = strip_ws(a.String());
    size_t      t     = my_mped_info.trait_find(tname);

    if(t == (size_t)-1)
    {
      errors << priority(error) << "Trait '" << tname << "' not found." << endl;
    }
    else
    {
      my_config.set_parent_trait(t);
    }
  }
  else
  {
    errors << priority(error) << "Parent trait parameter is missing name.  Skipping..." << endl;
  }
}


//===============================================
//
//  parse_sample(...)
//
//===============================================
void
Parser::parse_sample(const LSFBase* param)
{
  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    std::string s = strip_ws(toUpper(a.String()));

    if( s == "ALLELE" || s == "ALLELES" )
    {
      my_config.set_method(Configuration::ALLELES);
    }
    else if ( s == "GENO" || s == "GENOTYPES" )
    {
      my_config.set_method(Configuration::GENOTYPES);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'sample'." << endl;
    }
  }
  else
  {
    errors << priority(error) << "No value given for parameter 'sample'." << std::endl;
  }
}

} // End namespace TDTEX
} // End namespace SAGE
