//==========================================================================
//  File:       test_mcmc_parser.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              May. 04
//
//  Notes:      This file implements a parser for test_mcmc analysis.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/test_mcmc_parser.h"

namespace SAGE
{

namespace MCMC
{

test_mcmc_parser::test_mcmc_parser(cerrorstream&  err)
                : errors(err)
{}

test_mcmc_parser::test_mcmc_parser(const test_mcmc_parser& p)
                : errors(p.errors)
{
  my_parameters = p.my_parameters;
  my_regions = p.my_regions;
}

test_mcmc_parser::~test_mcmc_parser()
{}

void
test_mcmc_parser::parse_test_parameter_section(const LSFBase* params)
{
  if( !params || !params->List() )
  {
    return;
  }

  LSFList::const_iterator i = params->List()->begin();
  for( ; i != params->List()->end(); ++i )
  {
    const LSFBase* param = *i;

    if( !param || !param->name().size() ) continue;

    parse_test_parameter(param);
  }
}


void
test_mcmc_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size() )
    return;

  string n = toUpper( param->name() );

  if( n == "MODE" || n == "IBD_MODE" || n == "ANALYSIS_MODE" )
  {
    parse_mode(param);
  }
  else if( n == "REGION" )
  {
    parse_region(param);
  }

  return;
}

void
test_mcmc_parser::parse_region(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  string s;

  parse_string(param, s, errors);

  if( !s.size() ) return;

  AttrVal v = attr_value( param, "OUTPUT" );

  my_regions.push_back(test_mcmc_region_type(s, v.String()));
}        

void
test_mcmc_parser::parse_mode(const LSFBase* param)
{
  AttrVal v = attr_value(param, 0);

  if( v.has_value() )
  {
    string s = toUpper(v.String());

    if( s == "MULTI" || s == "MULTIPOINT" || s == "MULTI-POINT" || s == "MULTI_POINT" )
    {
      my_parameters.set_multipoint(true);
    }
    else if(s == "SINGLE" || s == "SINGLEPOINT" || s != "SINGLE-POINT" || s == "SINGLE_POINT" )
    {
      my_parameters.set_multipoint(false);
    }
    else
    {
      my_parameters.set_multipoint(true);

      errors << priority(warning) << "Mode " << v.String() 
             << " is unknown.  Mode will be ignored."
             << endl;
    }
  }
}

} // end of namespace MCMC

} // end of namespace SAGE

