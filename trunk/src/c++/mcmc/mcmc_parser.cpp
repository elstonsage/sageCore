//==========================================================================
//  File:       mcmc_parser.cpp
//
//  Author:     Qing Sun
//              Geoff Wedig
//
//  History:    Initial implementation.                        Aug. 28, 1998
//              1.0 Completely rewritten so makes sense        Feb   4, 2000
//              2.0 Updated to new libraries.                  yjs May. 2004
//
//  Notes:      This file implements a parser for mcmc analysis.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/mcmc_parser.h"

namespace SAGE
{

namespace MCMC
{

mcmc_parser::mcmc_parser(mcmc_parameters& p, SAGE::cerrorstream&  err)
           : SAGE::APP::BasicParser(err), my_parameters(p)
{
  clear();
}

mcmc_parser::mcmc_parser(mcmc_parameters& p, const SymbolTable* syms, SAGE::cerrorstream& err)
           : SAGE::APP::BasicParser(err), my_parameters(p)
{
  clear();
  parse_symbols(syms);
}

mcmc_parser::mcmc_parser(mcmc_parameters& p, const LSFBase* params, SAGE::cerrorstream& err)
           : SAGE::APP::BasicParser(err), my_parameters(p)
{
  clear();
  parse_test_parameter_section(params);
}

mcmc_parser::mcmc_parser(mcmc_parameters& p, const SymbolTable* syms,
                         const LSFBase* params, SAGE::cerrorstream& err)
           : SAGE::APP::BasicParser(err), my_parameters(p)
{
  clear();
  parse_symbols(syms);
  parse_test_parameter_section(params);
}

void
mcmc_parser::parse_symbols(const SymbolTable* syms)
{ 
  // This code needs adding eventually
}

void
mcmc_parser::parse_parameter(const LSFBase* param)
{ }

void
mcmc_parser::parse_test_parameter_section(const LSFBase* params)
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

  my_parameters.normalize_weights();
}

void mcmc_parser::parse_test_parameter(const LSFBase* param)
{
  if( !param || !param->name().size() )
    return;

  AttrVal a;

  string n = toUpper( param->name() );

  if(n == "MODE" || n == "IBD_MODE")
  {
    parse_mode(param);
  }
  else if(n == "DEMEMORIZATION_STEPS" || n == "DEMEM_STEPS")
  {
    parse_dememorization(param);
  }
  else if(n == "SIMULATION_STEPS" || n == "SIM_STEPS")
  {
    parse_sim_steps(param);
  }
  else if(n == "SIMULATION_BATCHES" || n == "SIM_BATCHES" || n == "BATCH_COUNT")
  {
    parse_batches(param);
  }
  else if(n == "USE_FACTOR" || n == "USE_FACTORING" || n == "USE_SCALING")
  {
    parse_factor(param);
  }
  else if(n == "BASE_FACTOR" || n == "BASE_SCALE" || n == "BASE_SCALING_FACTOR")
  {
    parse_base_factor(param);
  }
  else if(n == "DEMEMORIZATION_FACTOR" || n == "DEMEM_FACTOR" ||
          n == "DEMEMORIZATION_SCALING_FACTOR" || n == "DEMEM_SCALE")
  {
    parse_demem_factor(param);
  }
  else if(n == "SIMULATION_FACTOR" || n == "SIM_FACTOR" ||
          n == "SIMULATION_SCALING_FACTOR" || n == "SIM_SCALE")
  {
    parse_sim_factor(param);
  }
  else if(n == "SIMULATION_BATCH_FACTOR" || n == "SIM_BATCH_FACTOR" ||
          n == "BATCH_SCALING_FACTOR" || n == "BATCH_SCALE")
  {
    parse_batch_factor(param);
  }
  else if(n == "SIMULATION_LOCAL_MARKER" || n == "SIM_LOCAL_MARKER" ||
          n == "LOCAL_MARKER" )
  {
    parse_local_marker(param);
  }
  else if(n == "SIMULATION_LOCAL_INDIVIDUAL" || n == "SIM_LOCAL_IND" ||
          n == "LOCAL_IND" || n == "LOCAL_IND")
  {
    parse_local_ind(param);
  }
  else if(  n == "T0_WEIGHT" || n == "T0"
         || n == "T1_WEIGHT" || n == "T1"
         || n == "T2_WEIGHT" || n == "T2")
  {
    char c = n[1];

    int i = c - '0';

    parse_weights(i, param, n);
  }
  else if(n == "MAX_TUNNEL")
  {
    parse_tunnel(param);
  }
  else if(n == "SINGLE_MARKER_TUNNEL_WEIGHT" || n == "SINGLE_MARKER_WEIGHT")
  {
    double d;

    d = my_parameters.get_single_marker_tunneling();

    parse_real(param, d);

    my_parameters.set_single_marker_tunneling(d);
  }
  else if(n == "RANDOM_SEED")
  {
    parse_random_seed(param);
  }
}

void
mcmc_parser::parse_mode(const LSFBase* param)
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

void mcmc_parser::parse_dememorization(const LSFBase* param)
{
  int i = my_parameters.get_dememorization_step();

  parse_integer(param, i);

  my_parameters.set_dememorization_step(i);
}

void mcmc_parser::parse_sim_steps(const LSFBase* param)
{
  int i = my_parameters.get_simulation_step();

  parse_integer(param, i);

  my_parameters.set_simulation_step(i);
}

void mcmc_parser::parse_batches(const LSFBase* param)
{
  int i = my_parameters.get_batch_count();

  parse_integer(param, i);

  my_parameters.set_batch_count(i);
}

void mcmc_parser::parse_factor(const LSFBase* param)
{
  bool b = my_parameters.get_use_factor();

  parse_boolean(param, b);

  my_parameters.set_use_factor(b);
}

void mcmc_parser::parse_base_factor(const LSFBase* param)
{  
  double i = 0;

  parse_real(param, i);

  if(i > 0)
  {
    my_parameters.set_dememorization_factor(i * 0.5);
    my_parameters.set_simulation_factor(i * 10);
    my_parameters.set_batch_factor(i);
  }
}

void mcmc_parser::parse_demem_factor(const LSFBase* param)
{
  double i = my_parameters.get_dememorization_factor();

  parse_real(param, i);

  my_parameters.set_dememorization_factor(i);
}

void mcmc_parser::parse_sim_factor(const LSFBase* param)
{
  double i = my_parameters.get_simulation_factor();

  parse_real(param, i);

  my_parameters.set_simulation_factor(i);
}

void mcmc_parser::parse_batch_factor(const LSFBase* param)
{
  double i = my_parameters.get_batch_factor();

  parse_real(param, i);

  my_parameters.set_batch_factor(i);
}

void mcmc_parser::parse_local_marker(const LSFBase* param)
{
  double i = my_parameters.get_local_marker();

  parse_real(param, i);

  my_parameters.set_local_marker(i);
}

void mcmc_parser::parse_local_ind(const LSFBase* param)
{
  double i = my_parameters.get_local_individual();

  parse_real(param, i);

  my_parameters.set_local_individual(i);
}

void mcmc_parser::parse_tunnel(const LSFBase* param)
{
  int i = my_parameters.get_max_tunnel();

  parse_integer(param, i);

  my_parameters.set_max_tunnel(i);
}

void mcmc_parser::parse_weights(int i, const LSFBase* param, const string& name)
{
  double d;

  d = my_parameters.get_transition_weight(i);

  parse_real(param, d);

  my_parameters.set_transition_weight(i, d);
}
void mcmc_parser::parse_random_seed(const LSFBase* param)
{
  int i = my_parameters.get_random_seed();

  parse_integer(param, i);

  my_parameters.set_random_seed(i);
}

} // end of namespace MCMC

} // end of namespace SAGE
