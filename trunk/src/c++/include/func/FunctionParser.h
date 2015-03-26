#ifndef FUNCPARSER_H
#define FUNCPARSER_H

#include <iostream>
#include <string>
#include <algorithm>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "LSF/LSF.h"
#include "rped/rped.h"

#define DEFAULT_TIME_LIMIT 30

namespace SAGE {
namespace FUNC {

class FunctionParser
{
  public:

    struct TraitData
    {
      // A vector of pairs, where the first string is the name of the constant, and the second string
      // is the python expression that will calculate it.
      typedef vector<pair<string, string> >  constant_vector;

      TraitData();
      TraitData(const TraitData &);
      TraitData& operator= (const TraitData &);
      
      constant_vector  constants;
      string           trait_name;
      RPED::RefTraitInfo::trait_use  usage;  
      string           expr;
      unsigned int     time_limit;
      string           missing;
      bool             binary;
      string           affected;
      string           unaffected;
      bool             skip;
      bool             verbose;
    };
  
    // Reads in a function block, parses it, and returns a TraitData instance describing
    // that block.
    static TraitData create_trait_data(const LSFBase* function_block, cerrorstream errors);
    
  private:
    enum error_msg { duplicate_trait, invalid_param, no_constant_name, 
                     no_constant_expr, no_trait_name, no_trait_expr, 
                     bad_time_limit, constant_redef };  

    const TraitData&  get_data()  const { return my_data; }
    
    void  parse(const LSFBase* function, cerrorstream errors);
    void  parse_constant(const LSFList::const_iterator& iter, cerrorstream errors);
    void  parse_trait(const LSFList::const_iterator& iter, cerrorstream errors);
    void  parse_time_limit(const LSFList::const_iterator& iter, cerrorstream errors);
    void  write_error_msg(error_msg msg, cerrorstream errors, string string_var = "" ) const;
  
    // Data members.
    TraitData  my_data;
};


inline
FunctionParser::TraitData::TraitData() 
    : trait_name(""), usage(RPED::RefTraitInfo::unknown_use), expr(""), 
      time_limit(DEFAULT_TIME_LIMIT), missing(""), binary(false), 
      affected(""), unaffected(""), skip(false), verbose(false)
{ 
  constants.clear();
}

inline 
FunctionParser::TraitData::TraitData(const TraitData& other) 
    : constants(other.constants), trait_name(other.trait_name), 
      usage(other.usage), expr(other.expr), time_limit(other.time_limit), 
      missing(other.missing), binary(other.binary), affected(other.affected), 
      unaffected(other.unaffected), skip(other.skip), verbose(other.verbose)
{}


inline
FunctionParser::TraitData&
FunctionParser::TraitData::operator=(const TraitData& other)
{
  if(this != &other)
  {
    constants  = other.constants;
    trait_name = other.trait_name;
    usage      = other.usage;
    expr       = other.expr;
    time_limit = other.time_limit;
    missing    = other.missing;
    binary     = other.binary;
    affected   = other.affected;
    unaffected = other.unaffected;
    skip       = other.skip;
    verbose    = other.verbose;
  }
  
  return *this;
}

} // End namespace FUNC
} // End namespace SAGE

#endif
