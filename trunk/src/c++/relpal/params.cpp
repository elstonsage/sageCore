//=============================================================================
// File:    params.cpp
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                  yjs Mar. 07
//
// Notes:   This file contains implementation for following data structures.
//            struct dependent_variable
//            struct covariate_type
//            struct marker_type
//            struct independent_variable
//
// Copyright (c) 2007 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/params.h"

using namespace std;

namespace SAGE   {
namespace RELPAL {

dependent_variable::dependent_variable(size_t t, bool bin)
{
  trait_index = t;
  binary      = bin;

  valid = false;
}

string
dependent_variable::name(const relative_pairs& pairs) const
{
  string t = pairs.trait_name(trait_index);

  if( !t.size() )
    t = "Trait #" + long2str(trait_index);

  return t;
}

//
//------------------------------------------------------------------------
//

covariate_type::covariate_type(size_t c, size_t t, bool test)
              : covariate_index(c), adj_trait_index(t),
                test_variable(test), valid(false)
{}

string
covariate_type::name(const relative_pairs& pairs) const
{
  string n = pairs.trait_name(covariate_index);

  if( !n.size() )
    n = "Covariate #" + long2str(covariate_index);

  return n;
}

string
covariate_type::effect_name() const
{
  string n = "Cov";

  return n;
}

string
covariate_type::short_effect_name() const
{
  string n = "C";

  return n;
}

//
//------------------------------------------------------------------------
//

marker_type::marker_type(size_t m, bool test)
           : marker_index(m), effect(TOTAL), test_variable(test), x_linked(false), valid(false)
{}

string
marker_type::name(const relative_pairs& pairs) const
{
  string n = pairs.marker_name(marker_index);

  if( !n.size() )
    n = "Marker #" + long2str(marker_index);

  return n;
}

string
marker_type::effect_name() const
{
  string n = "(A+D)GenVar";

  return n;
}


string
marker_type::short_effect_name() const
{
  string n = "A+D";

  return n;
}

//
//------------------------------------------------------------------------
//

string covariate_name(const relative_pairs& pairs, const vector<covariate_type>& covariates)
{
  string n;

  if( covariates.size() )
    n = covariates[0].name(pairs);

  for(size_t i = 1; i < covariates.size(); ++i)
  {
    n += " x " + covariates[i].name(pairs);
  }
  return n;
}

string marker_name(const relative_pairs& pairs, const vector<marker_type>& markers)
{
  if( !markers.size() )
    return "";

  string n = markers[0].name(pairs);
  for( size_t i = 1; i < markers.size(); ++i )
  {
    n += " x " + markers[i].name(pairs);
  }

  return n;
}

string covariate_effect_name(bool shortname, const vector<covariate_type>& covariates)
{
  if( !covariates.size() )
    return "";

  if( !shortname && covariates.size() == 1)
    return covariates[0].effect_name();

  string n = covariates[0].short_effect_name();
  for( size_t i = 1; i < covariates.size(); ++i )
    n += " x " + covariates[i].short_effect_name();

  return n;
}

string marker_effect_name(bool shortname, const vector<marker_type>& markers)
{
  if( !markers.size() )
    return "";

  if( !shortname && markers.size() == 1)
    return markers[0].effect_name();

  string n = markers[0].short_effect_name();
  for( size_t i = 1; i < markers.size(); ++i )
    n += " x " + markers[i].short_effect_name();

  return n;
}

//
//------------------------------------------------------------------------
//

interaction_type::interaction_type()
                : type(CC), batch(false), test_variable(false), valid(false)
{ 
  markers.resize(0);
  covariates.resize(0);
}

string
interaction_type::name(const relative_pairs& pairs) const
{
  string n = "";

  switch( type )
  {
    case MM    : n = marker_name(pairs, markers); break;
    case CC    : n = covariate_name(pairs, covariates); break;
    case MC    : n = marker_name(pairs, markers) + " x " + covariate_name(pairs, covariates); break;
    default    : n = ""; break;    
  }

  return n;
}

string
interaction_type::effect_name() const
{
  string n = "";

  switch( type )
  {
    case MM    : n = marker_effect_name(true, markers); break;
    case CC    : n = covariate_effect_name(true, covariates); break;
    case MC    : n = marker_effect_name(true, markers) + " x " + covariate_effect_name(true, covariates); break;
    default    : n = ""; break;    
  }

  return n;
}

string
interaction_type::short_effect_name() const
{
  string n = "";

  switch( type )
  {
    case MM    : n = marker_effect_name(true, markers); break;
    case CC    : n = covariate_effect_name(true, covariates); break;
    case MC    : n = marker_effect_name(true, markers) + " x " + covariate_effect_name(true, covariates); break;
    default    : n = ""; break;    
  }

  return n;
}

//
//------------------------------------------------------------------------
//

independent_variable::independent_variable()
{ 
  clear();
}

void
independent_variable::clear()
{
  markers.resize(0);
  covariates.resize(0);

  type          = COVARIATE;
  test_variable = false;
  valid         = false;
}

string
independent_variable::dump(const relative_pairs& pairs, bool with_t) const
{
  return name(pairs) + " (" + effect_name(with_t) + ")";
}

string
independent_variable::name(const relative_pairs& pairs) const
{
  string n = "";
  switch( type )
  {
    case INTERCEPT   : n = "INTERCEPT"; break;
    case RANDOM_ERR  : n = "RANDOM_EFF"; break;
    case COMMON_ENV  : n = "COMMON_ENV_EFF"; break;
    case POLYGENIC   : n = "POLYGENIC_EFF"; break;
    case MARKER      : n = marker_name(pairs, markers); break;
    case COVARIATE   : n = covariate_name(pairs, covariates); break;
    case MM_INTER    : n = marker_name(pairs, markers); break;
    case CC_INTER    : n = covariate_name(pairs, covariates); break;
    case MC_INTER    : n = marker_name(pairs, markers) + " x " + covariate_name(pairs, covariates); break;
    default          : n = ""; break;    
  }

  return n;
}

string
independent_variable::effect_name(bool with_t) const
{
  string n = "";

  switch( type )
  {
    case INTERCEPT   : n = "INTERCEPT"; break;
    case RANDOM_ERR  : n = "RANDOM_EFF"; break;
    case COMMON_ENV  : n = "COMMON_ENV_EFF"; break;
    case POLYGENIC   : n = "POLYGENIC_EFF"; break;
    case MARKER      : n = marker_effect_name(false, markers); break;
    case COVARIATE   : n = covariate_effect_name(false, covariates); break;
    case MM_INTER    : n = marker_effect_name(true, markers); break;
    case CC_INTER    : n = covariate_effect_name(true, covariates); break;
    case MC_INTER    : n = marker_effect_name(true, markers) + " x " + covariate_effect_name(true, covariates); break;
    default          : n = ""; break;    
  }

  if( with_t )
    n = n + "(t" + long2str(t1) + ",t" + long2str(t2) + ")";

  return n;
}

string
independent_variable::short_effect_name(bool with_t) const
{
  string n = "";

  switch( type )
  {
    case INTERCEPT   : n = "INTERCEPT"; break;
    case RANDOM_ERR  : n = "RANDOM_EFF"; break;
    case COMMON_ENV  : n = "COMMON_ENV_EFF"; break;
    case POLYGENIC   : n = "POLYGENIC_EFF"; break;
    case MARKER      : n = marker_effect_name(true, markers); break;
    case COVARIATE   : n = covariate_effect_name(true, covariates); break;
    case MM_INTER    : n = marker_effect_name(true, markers); break;
    case CC_INTER    : n = covariate_effect_name(true, covariates); break;
    case MC_INTER    : n = marker_effect_name(true, markers) + " x " + covariate_effect_name(true, covariates); break;
    default          : n = ""; break;    
  }

  if( with_t )
    n = n + "(t" + long2str(t1) + ",t" + long2str(t2) + ")";

  return n;
}

//
// ----------------------------------------------------------------------
//

filtering_type::filtering_type(size_t s, string n)
              : subset_index(s), subset_name(n)
{}

string
filtering_type::name() const
{
  return subset_name;
}

//
// ----------------------------------------------------------------------
//

data_options::data_options()
            : use_members(DEFAULT), use_pairs(ALL)
{
  subsets.resize(0);
}

string
data_options::dump(string pre_space) const
{
  string n = "";
  if( subsets.size() )
  {
    n += pre_space + "subset:";
    for( size_t i = 0; i < subsets.size(); ++i )
      n += " " + subsets[i].name();
    n += "\n";
  }

  n += pre_space + "use_members = ";
  if( use_members == INFORMATIVE_LOCAL )
    n += "locally informative\n";
  else if( use_members == INFORMATIVE_REGION )
    n += "regionally informative\n";
  else
    n += "every\n";

  n += pre_space + "use_pairs = ";
  if( use_pairs == FSIB )
    n += "full sib\n";
  else if( use_pairs == HSIB )
    n += "half sib\n";
  else if( use_pairs == SIB )
    n += "both sib\n";
  else
    n += "all\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

output_options::output_options()
              : detailed_out(false), export_out(false), data_out(false),
                debug_out(false), residual_out(false)
{
  param_name_max     = DEFAULT_PARAM_NAME_MAX;
  param_eff_name_max = DEFAULT_PARAM_NAME_MAX;
}

string
output_options::dump(string pre_space) const
{
  string n = pre_space + "detailed_out = ";

  if( detailed_out )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "export_out = ";
  if( export_out )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "data_out = ";
  if( data_out )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "debug_out = ";
  if( debug_out )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "residual_out = ";
  if( residual_out )
    n += "yes\n";
  else
    n += "no\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

analysis_options::analysis_options()
                : transform_residuals(true), naive_variance(true), sandwich_variance(true),
                  alternative_variance(true), IBD_variance(true), state_file_name("")
{}

string
analysis_options::dump(string pre_space) const
{
  string n = pre_space + "transform_residuals = ";

  if( transform_residuals )
    n += "yes\n";
  else
    n += "no\n";

  n += dump_var(pre_space);

  return n;
}

string
analysis_options::dump_var(string pre_space) const
{
  string n = pre_space + "naive_variance = ";

  if( naive_variance )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "sandwich_variance = ";
  if( sandwich_variance )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "alternative_variance = ";
  if( alternative_variance )
    n += "yes\n";
  else
    n += "no\n";

  n += pre_space + "IBD_variance = ";
  if( IBD_variance )
  {
    n += "yes";
    if( state_file_name.size() )
      n += (", file name = " + state_file_name);
    n += "\n";
  }
  else
    n += "no\n";

  return n;
}

//
// ----------------------------------------------------------------------
//

asymptotic_pvalue_options::asymptotic_pvalue_options()
                         : valid(false), seed(0)
{
  replicates     = 0;
  min_replicates = 20;
  max_replicates = 10000;

  threshold  = 0.05;
  width      = 0.2;
  confidence = 0.95;
}

string
asymptotic_pvalue_options::dump(string pre_space) const
{
  string n = pre_space + "seed = " + long2str(seed) + "\n";

  n += pre_space + "replicates = " + long2str(replicates) + "\n";
  n += pre_space + "min_replicates = " + long2str(min_replicates) + "\n";
  n += pre_space + "max_replicates = " + long2str(max_replicates) + "\n";

  n += pre_space + "threshold = " + doub2str(threshold) + "\n";
  n += pre_space + "width = " + doub2str(width) + "\n";
  n += pre_space + "confidence = " + doub2str(confidence) + "\n";

  return n;
}

} // end of namespace RELPAL
} // end of namespace SAGE
