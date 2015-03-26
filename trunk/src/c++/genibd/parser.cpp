//==========================================================================
//  File:       parser.cpp
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              Nov. 03
//
//  Notes:      This file implements a parser for genibd analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/parser.h"

namespace SAGE
{

namespace GENIBD
{

genibd_parser::genibd_parser(cerrorstream&  err)
             : my_mcmc_parser(my_parameters.get_sim_parameters()), errors(err)
{}

genibd_parser::~genibd_parser()
{}

void
genibd_parser::parse_test_parameter_section(const LSFBase* params)
{
  if( !params || !params->List() )
  {
    return;
  }

  // Get the output attribute

  if( params->attrs() )
  {
    AttrVal v = attr_value(params,"OUT");

    if( !v.has_value() || v.String().empty() )
    {
      v = attr_value(params,"OUTPUT");
    }
     
    if( v.has_value() && !v.String().empty() )
    {
      string s = v.String();

      my_parameters.set_output(s);
    }
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
genibd_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size() )
    return;

  string n = toUpper( param->name() );

  if( n == "TITLE" )
  {
    string s = my_parameters.title();

    parse_string(param, s);

    my_parameters.set_title(s);
  }
  else if( n == "MAX_PEDIGREE" || n == "MAX_CUTOFF" || n == "MAX_SIZE"
                               || n == "MAX_BITS" )
  {
    int i = my_parameters.max_exact_size();

    parse_integer(param, i);

    my_parameters.set_max_exact_size(i);
  }
  else if( n == "MODE" || n == "IBD_MODE" || n == "ANALYSIS_MODE" )
  {
    parse_mode(param);
  }
  else if( n == "SCAN_TYPE" || n == "SCAN" )
  {
    parse_scan_type(param);
  }    
  else if( n == "ALLOW_LOOPS" )
  {
    parse_loops(param);
  }    
  else if( n == "OUTPUT_IBD_STATE" || n == "IBD_STATE" || n == "IBD_STATE_OUT"  )
  {
    parse_output_ibd_state(param);
  }    
  else if( n == "USE_SIMULATION" || n == "SIMULATION" )
  {
    parse_simulation(param);
  }    
  else if( n == "SPLIT_PEDIGREES" )
  {
    parse_family(param);
  }    
  else if( n == "OUTPUT_PAIR_TYPES" )
  {
    parse_pair_types(param);
  }
  else if( n == "REGION" )
  {
    parse_region(param);
  }
  else
  {
    my_mcmc_parser.parse_test_parameter(param);
//    my_parameters.get_sim_parameters() = my_mcmc_parser.parameters();
  }

  return;
}

void
genibd_parser::parse_region(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  string s;

  parse_string(param, s, errors);

  if( !s.size() ) return;

  AttrVal v = attr_value( param, "OUTPUT" );

  my_regions.push_back(genibd_region_type(s, v.String()));
}        

void
genibd_parser::parse_mode(const LSFBase* param)
{
  AttrVal v = attr_value(param, 0);

  if( v.has_value() )
  {
    string s = toUpper(v.String());

    if( s == "MULTI" || s == "MULTIPOINT" || s == "MULTI-POINT" || s == "MULTI_POINT" )
    {
      my_parameters.set_multipoint(true);
      my_parameters.get_sim_parameters().set_multipoint(true);
    }
    else if(s == "SINGLE" || s == "SINGLEPOINT" || s != "SINGLE-POINT" || s == "SINGLE_POINT" )
    {
      my_parameters.set_multipoint(false);
      my_parameters.get_sim_parameters().set_multipoint(false);
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

void
genibd_parser::parse_scan_type(const LSFBase* param)
{
  AttrVal v = attr_value(param, 0);

  if( v.has_value() )
  {
    string s = toUpper(v.String());

    bool intervals = s == "I" || s == "INTERVAL" || s == "INTERVALS";

    my_parameters.set_scan_interval(intervals);

    if( param->attrs() )
    {
      AttrVal v = attr_value(param,"DISTANCE"); // "SCAN_DISTANCE", "INTERVAL_DISTANCE"

      if( !v.has_value() || v.String().empty() )
      {
        v = attr_value(param,"INTERVAL_DISTANCE");
      }
       
      if( v.has_value() && !v.String().empty() )
      {
        double distance = v.Real();

        my_parameters.set_interval_distance(distance);
      }
    }
  }
}

void genibd_parser::parse_loops(const LSFBase* param)
{
  bool b = my_parameters.allow_loops();

  parse_boolean(param, b);

  my_parameters.set_allow_loops(b);
}

void genibd_parser::parse_output_ibd_state(const LSFBase* param)
{
  bool b = my_parameters.output_ibd_state();

  parse_boolean(param, b);

  my_parameters.set_output_ibd_state(b);
}

void
genibd_parser::parse_simulation(const LSFBase* param)
{
  choice_type c = my_parameters.allow_simulation();

  parse_choice(param, c, "USE_SIMULATION");

  my_parameters.set_allow_simulation(c);
}

void genibd_parser::parse_family(const LSFBase* param)
{
  choice_type c = my_parameters.allow_family_splitting();

  parse_choice(param, c, "SPLIT_PEDIGREES");

  my_parameters.set_allow_family_splitting(c);
}

void genibd_parser::parse_pair_types(const LSFBase* param)
{
  AttrVal v = attr_value(param, 0);

  if( v.has_value() )
  {
    string s = toUpper(v.String());
    
    if( s == "SIBLINGS" || s == "SIBS" || s == "SIB" || s == "S" )
    {
      my_parameters.set_pair_category(SIB);
    }
    else if( s == "ALL SIBS" || s == "ALL_SIBS" || s == "H" || s == "A" )
    {
      my_parameters.set_pair_category(ALL_SIB);
    }
    else if( s == "RELATIVE" || s == "REL" || s == "R" || s == "RELATIVES" )
    {
      my_parameters.set_pair_category(RELATIVE);
    }
    else if( s == "ALL PAIRS" || s == "ALL_PAIRS" || s == "ALLPAIRS")
    {
      my_parameters.set_pair_category(ALL);
    }
    else
    {
      errors << priority(warning) << "Unable to evaluate pair types \"" << v.String()
             << "\".  Type will be ignored." << endl;
    }
  }
}

void
genibd_parser::parse_choice(const LSFBase* param,
                            choice_type&   value,
                            const string&  name)
{
  AttrVal a=attr_value(param, 0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());
    if( n == "TRUE" || n == "YES" )
      value = YES;
    else if( n == "FALSE" || n == "NO" )
      value = NO;
    else if( n == "ALWAYS" )
      value = ALWAYS;
    else
      errors << priority(error) << "Unknown value for parameter '"
                                << name << "'" << endl;
  }
}

} // end of namespace GENIBD

} // end of namespace SAGE
