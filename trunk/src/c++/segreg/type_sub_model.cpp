//============================================================================
// File:      type_sub_model.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   djb created.                                         6/13/01
//            djb renamed, factored out genotype_sub_model  
//                and added genotype_specific_variance_sub_model.  7-30-01               
//                                                          
// Notes:     implementation of genotype_mean_specific_sub_model 
//            and genotype_specific_variance_sub_model classes.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "segreg/type_sub_model.h"

using namespace std;
namespace SAGE
{
namespace SEGREG
{
const std::string  TYPE_MEAN_NAME        = "Type means";
const std::string  TYPE_SUSCEPT_NAME     = "Type susceptibilities";
const std::string  TYPE_VAR_NAME         = "Type variances";

const double  MEAN_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN();  // 0;
const bool    MEAN_DEFAULT_FIXED = false;

const double  VAR_DEFAULT_VALUE = numeric_limits<double>::quiet_NaN();  // 1;
const bool    VAR_DEFAULT_FIXED = false;
const double  VAR_EPSILON = .00001;
const double  VAR_LB = 0;

//============================================================================
// IMPLEMENTATION:  genotype_specific_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
// - All three types have same value.
//
bool
genotype_specific_sub_model::set_one
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB)
{
  // Initially set parameters to defaults
  
  set_defaults();

  // Determine how much info we have.
  
  genotype_info  info = total_info(type_AA, type_AB, type_BB);
  
  bool return_value = false;
  switch(info)
  {
    case no_geno:
      return_value = true;
      break;
      
    case AA:
      return_value = set_one_one_value(type_AA);
      break;
      
    case AB:
      return_value = set_one_one_value(type_AB);
      break;
      
    case BB:
      return_value = set_one_one_value(type_BB);
      break;
      
    case AA_AB:
      return_value = set_one_two_values(type_AA, type_AB);
      break;
      
    case AA_BB:
      return_value = set_one_two_values(type_AA, type_BB);
      break;
      
    case AB_BB:
      return_value = set_one_two_values(type_AB, type_BB);
      break;
      
    case all:
      if(/*lint --e(777) */ type_AA.value == type_AB.value && type_AA.value == type_BB.value)
      {
        if(type_AA.fixed == type_AB.fixed && type_AA.fixed == type_BB.fixed)
        {
          my_default = false;

          my_types       [index_AA] = type_AA.value;
          my_types_fixed [index_AA] = type_AA.fixed;

          return_value = true;
        }
        else
        {
          my_errors << priority(critical) << "Conflicting data for genotype specific "
                    << long_plural_type() << ".  Skipping analysis ... " << endl;
        }
      }
      else
      {
        if(! (type_AA.fixed || type_AB.fixed || type_BB.fixed))
        {
          my_default = false;

          my_types       [index_AA] = type_AA.value;
          my_types_fixed [index_AA] = false;

          my_errors << priority(error) << "Conflicting data for genotype specific "
                    << long_plural_type() << ".  Making arbitrary choice (AA) ..." << endl;

          return_value = true;
        }
        else
        {
          my_errors << priority(critical) << "Conflicting data for genotype specific "
                    << long_plural_type() << ".  Skipping analysis ... " << endl;
        }
      }
      break;
      
    default:
      SAGE_internal_error();
  }
  
  // - Sync. the types with the estimate

  my_types[index_AB] =
  my_types[index_BB] = my_types[index_AA];

  my_types_fixed[index_AB] =
  my_types_fixed[index_BB] = my_types_fixed[index_AA];
  
  return return_value;
}

// - Types AA and AB are equal.
//
bool
genotype_specific_sub_model::set_two_dom
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB)
{
  // Assume we're good, unless told otherwise.

  bool return_value = true;

  // Initially set parameters to defaults
  
  set_defaults();

  // Determine how much info we have
  
  genotype_info  info = total_info(type_AA, type_AB, type_BB);
  
  // First, lets deal with BB, the easy case
  
  if(info & BB)
  {
    my_default = false;
    my_types      [2] = type_BB.value;
    my_types_fixed[2] = type_BB.fixed;
  }

  // Deal with AA and AB cases.  If just one set, use it, if both set, resolve
  // the conflict.
  
  switch(info & AA_AB)
  {
    case no_geno:
      break;
      
    case AA:
      my_default = false;

      my_types[0]       = type_AA.value;
      my_types_fixed[0] = type_AA.fixed;
      break;
      
    case AB:
      my_default = false;

      my_types[0]       = type_AB.value;
      my_types_fixed[0] = type_AB.fixed;
      break;
      
    case AA_AB:
      if(/*lint --e(777) */ type_AA.value == type_AB.value)
      {
        if(type_AA.fixed == type_AB.fixed)
        {
          my_default = false;

          my_types       [0] = type_AA.value;
          my_types_fixed [0] = type_AA.fixed;
        }
        else
        {
          my_errors << priority(critical) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AA and AB.  Skipping analysis ... " << endl;
          return_value = false;
        }
      }
      else
      {
        if(! (type_AA.fixed || type_AB.fixed))
        {
          my_default = false;

          my_types[0]       = type_AA.value;
          my_types_fixed[0] = type_AA.fixed;
          my_errors << priority(error) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AA and AB.  Making arbitrary choice ..." << endl;
        }
        else
        {
          my_errors << priority(critical) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AA and AB.  Skipping analysis ... " << endl;
          return_value = false;
        }
      }
      break;
    default:
      SAGE_internal_error();
  }
  
  // Copy the AA status into AB

  my_types      [index_AB] = my_types      [index_AA];
  my_types_fixed[index_AB] = my_types_fixed[index_AA];
  
  return return_value;
}

// - May be interpreted as either A dominant or A recessive depending
//   on the value of the flag, my_two_is_dom.
//
bool
genotype_specific_sub_model::set_two
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB)
{
  model_input tab = type_AB;

  if(! SAGE::isnan(type_AB.value))
  {
    double alternate_value = my_two_is_dom ? type_AA.value : type_BB.value;

    if(SAGE::isnan(alternate_value) || type_AB.value != alternate_value)
    {
      my_errors << priority(warning) << "Genotype specific " << long_sing_type()
                << " AB specified in " << type_param() 
                << " sub-block with option two.  Ignoring ..." << endl;
      tab.value = QNAN;
    }
  }

  return  my_two_is_dom ? set_two_dom(type_AA, tab, type_BB) :
                          set_two_rec(type_AA, tab, type_BB)  ;
}


bool
genotype_specific_sub_model::set_three(const model_input& type_AA, const model_input& type_AB,
                                            const model_input& type_BB)
{
  bool return_value =  set_three(type_AA, type_AB, type_BB, MAXFUN::Parameter::INDEPENDENT);
  
  return return_value;
}

// - Types BB and AB are equal.
//
bool
genotype_specific_sub_model::set_two_rec
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB)
{
  // Assume we're good unless told otherwise.
  
  bool return_value = true;

  // Initially set parameters to defaults
  
  set_defaults();

  // Determine how much info we have
  
  genotype_info  info = total_info(type_AA, type_AB, type_BB);
  
  // First, lets deal with AA, the easy case
  
  if(info & AA)
  {
    my_default = false;
    my_types      [index_AA] = type_AA.value;
    my_types_fixed[index_AA] = type_AA.fixed;
  }

  switch(info & AB_BB)
  {
    case no_geno:
      break;
      
    case AB:
      my_default = false;

      my_types      [index_BB] = type_AB.value;
      my_types_fixed[index_BB] = type_AB.fixed;
      break;
      
    case BB:
      my_default = false;

      my_types      [index_BB] = type_BB.value;
      my_types_fixed[index_BB] = type_BB.fixed;
      break;
      
    case AB_BB:
      if(/*lint --e(777) */ type_AB.value == type_BB.value)
      {
        if(type_AB.fixed == type_BB.fixed)
        {
          my_default = false;

          my_types      [index_BB] = type_AB.value;
          my_types_fixed[index_BB] = type_AB.fixed;
        }
        else
        {
          my_errors << priority(error) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AB and BB.  Skipping analysis ... " << endl;
          return_value = false;
        }
      }
      else
      {
        if(! (type_AB.fixed || type_BB.fixed))
        {
          my_default = false;

          my_types[index_BB] = type_AB.value;
          my_types_fixed[index_BB] = type_AB.fixed ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT;
          my_errors << priority(error) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AB and BB.  Making arbitrary choice ..." << endl;
        }
        else
        {
          my_errors << priority(critical) << "Conflicting data for genotype specific "
                    << long_plural_type() << " AB and BB.  Skipping analysis ... " << endl;
          return_value = false;
        }
      }
      break;
      
    default:
      SAGE_internal_error();
  }
  
  // Copy the BB status into AB

  my_types      [index_AB] = my_types      [index_BB];
  my_types_fixed[index_AB] = my_types_fixed[index_BB];
  
  return return_value;
}

// - 'AB' parameter is the average of the other two parameters.
//
bool
genotype_specific_sub_model::set_three_add
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB)
{
  bool return_value = true;

  set_defaults();

  // Determine our info
  
  genotype_info  info = total_info(type_AA, type_AB, type_BB);

  if(info & AA)
  {
    my_default = false;
    my_types[index_AA] = type_AA.value;
    my_types_fixed[index_AA] = type_AA.fixed;
  }

  if(info & BB)
  {
    my_default = false;
    my_types[index_BB] = type_BB.value;
    my_types_fixed[index_BB] = type_BB.fixed;
  }

  if(info & AB)
  {
    my_errors << priority(error) <<  "Genotype specific " << long_sing_type()
              << " AB specified in " << type_param() << " sub-block with option three_add.  Ignoring ..."
              << endl;
  }

  if(info & AA_BB)
  {
    set_three_add_set_AB();
  }
  
  return return_value;
}

//
// ===============  Ancillary Functions
//
// - Return a shorter version of type word.
//
string
genotype_specific_sub_model::short_sing_type() const
{
  switch(my_derived_type)
  {
    case  mean     : return "mean";
    case  variance : return "variance";
    case  susc     : return "susc";
    case  invalid  :
    default        : return "";
  }
}

// - Return plural form of type word.
//
string
genotype_specific_sub_model::long_plural_type() const
{
  switch(my_derived_type)
  {
    case  mean     : return "means";
    case  variance : return "variances";
    case  susc     : return "susceptibilities";
    case  invalid  :
    default        : return "";
  }
}

// - Return type_parameter corresponding to my_derived_type.
//
string 
genotype_specific_sub_model::type_param() const
{
  return type_param(my_derived_type); 
}

// - Return type_parameter corresponding to derived_type.
//
string 
genotype_specific_sub_model::type_param(derived_type t)
{
  string  param("type_");
  switch(t)
  {
    case  mean     : return param += "mean";
    case  variance : return param += "var";
    case  susc     : return param += "suscept";
    case  invalid  :
    default        : SAGE_internal_error();
  }

  // This can never happen
  return param; 
}

// - If genotype variance sub-model, set default value and lower bound.
//   VAR_EPSILON needed because maxfun treats bounds as inclusive and 
//   here we want the bound to be treated as exclusive.
// - If genotype mean sub-model, set default value.
//
void
genotype_specific_sub_model::set_defaults()
{
  for(int i = 0; i < NUM_OF_TYPES; ++i)
  {
    my_types[i] = QNAN;
    my_types_fixed[i] = false;
  }
}

// - Set parameters when option is one and one model input has a 
//   value.
//
bool
genotype_specific_sub_model::set_one_one_value(const model_input& mi)
{
  my_default = false;

  my_types       [index_AA] = mi.value;
  my_types_fixed [index_AA] = mi.fixed;
  return true;
}

// - Set parameters when option is one and two model inputs have
//   values.
//
bool
genotype_specific_sub_model::set_one_two_values(const model_input& mi_one,
                                                     const model_input& mi_two )
{
  bool  return_value = false;
  
  if(/*lint --e(777) */ mi_one.value == mi_two.value)
  {
    if(mi_one.fixed == mi_two.fixed)
    {
      my_default = false;

      my_types[0] = mi_one.value;
      my_types_fixed[0] = mi_one.fixed;

      return_value = true;
    }
    else
    {
      my_errors << priority(critical) << "Conflicting data for genotype specific "
                << long_plural_type() << ".  Skipping analysis ... " << endl;
    }
  }
  else
  {
    if(! (mi_one.fixed || mi_two.fixed))
    {
      my_default = false;

      my_types[0] = mi_one.value;
      my_types_fixed[0] = mi_one.fixed;

      return_value = true;

      my_errors << priority(error) << "Conflicting data for genotype specific "
                << long_plural_type() << ".  Making arbitrary choice ... " << endl;
    }
    else
    {
      my_errors << priority(critical) << "Conflicting data for genotype specific "
                << long_plural_type() << ".  Skipping analysis ... " << endl;
    }
  }
  
  return return_value;
}

// - Used by set three add to set value of 'AB' parameter as 
//   a function of the other two parameters.
//
void
genotype_specific_sub_model::set_three_add_set_AB()
{
  my_types[index_AB] = (my_types[index_AA] + my_types[index_BB]) / 2.0;
  my_types_fixed[index_AB] = my_types_fixed[index_AA] && my_types_fixed[index_BB];
}


// - Handles options three, three_dec, and three_inc which differ
//   only in their constraints.
//
bool
genotype_specific_sub_model::set_three
      (const model_input& type_AA, const model_input& type_AB, const model_input& type_BB,
       bool non_fixed_status)
{
  set_defaults();

  genotype_info  info = total_info(type_AA, type_AB, type_BB);
  
  if(info != no_geno)
    my_default = false;

  if(! SAGE::isnan(type_AA.value))
  {
    my_types[0] =  type_AA.value;
    my_types_fixed[0] = type_AA.fixed;
  }

  if(! SAGE::isnan(type_AB.value))
  {
    my_types[1]  = type_AB.value;
    my_types_fixed[1] = type_AB.fixed;
  }

  if(! SAGE::isnan(type_BB.value))
  {
    my_types[2] =  type_BB.value;
    my_types_fixed[2] = type_BB.fixed;
  }
  
  return true;
}

int
genotype_specific_sub_model::finalizeConfiguration()
{
  switch(my_option)
  {
    case one:
      initialize_one();
      break;
    case two:
      if(is_two_dom())
        initialize_two_dom();
      else
        initialize_two_rec();
      break;
    case three:
      initialize_three();
      break;
    case two_dom:
      initialize_two_dom();
      break;
    case two_rec:
      initialize_two_rec();
      break;
    case three_add:
      initialize_three_add();
      break;
    case three_dec:
    case three_inc:
      // For three_dec and three_inc, they have the same values as the basic
      // three model, but have a functional relationship to one another.
      initialize_three();
      
      // Redefine any non-fixed as independent functional.
      
      my_parameters[0].initial_type = my_types_fixed[0] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[1].initial_type = my_types_fixed[1] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;
      my_parameters[2].initial_type = my_types_fixed[2] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL;

      break;
  }

  return 0;
}

void
genotype_specific_sub_model::initialize_one()
{
  // - Set default values for option one.
  //
  my_parameters.resize(1);
  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type(),
      my_types_fixed[index_AA] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_AA],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
}

void
genotype_specific_sub_model::initialize_two_dom()
{
  // - Set default values for option two.
  //
  my_parameters.resize(2);
  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AA_AB",
      my_types_fixed[index_AA] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_AA],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);

  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_BB",
      my_types_fixed[index_BB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_BB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
}

void
genotype_specific_sub_model::initialize_two_rec()
{ 
  // - Set default values for option two_rec.
  //
  my_parameters.resize(2);
  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AA",
      my_types_fixed[index_AA] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_AA],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AB_BB",
      my_types_fixed[index_BB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_BB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
}

// - 'AB' parameter is the average of the other two parameters.
//
void
genotype_specific_sub_model::initialize_three_add()
{
  // - Set default values for option three_add.
  //
  my_parameters.resize(3);
  my_parameters[index_AA] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AA",
      my_types_fixed[index_AA] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[index_AA],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
  my_parameters[index_AB] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AB",
      my_types_fixed[index_AB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::DEPENDENT,
      my_types[index_AB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
  my_parameters[index_BB] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_BB",
      my_types_fixed[index_BB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL,
      my_types[index_BB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
}

void
genotype_specific_sub_model::initialize_three()
{
  // - Set default values for option three.
  //
  my_parameters.resize(3);
  my_parameters[0] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AA",
      my_types_fixed[index_AA] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_AA],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
  my_parameters[1] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_AB",
      my_types_fixed[index_AB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_AB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
  my_parameters[2] = MAXFUN::ParameterInput
    ( toUpper(long_plural_type()),
      short_sing_type() + "_BB",
      my_types_fixed[index_BB] ? MAXFUN::Parameter::FIXED : MAXFUN::Parameter::INDEPENDENT,
      my_types[index_BB],
      (my_derived_type == variance)? VAR_LB + VAR_EPSILON : NEGATIVE_INF,
      POSITIVE_INF);
}

//
// ==============  Synchronization w. Maxfun
//
int
genotype_specific_sub_model::synchronize_one()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(0);
  my_types[index_BB] = getParam(0);

  return 0;
}

int  
genotype_specific_sub_model::synchronize_two()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(my_two_is_dom ? 0 : 1);
  my_types[index_BB] = getParam(1);

  return 0;
}
 
int  
genotype_specific_sub_model::synchronize_three()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(1);
  my_types[index_BB] = getParam(2);

  return 0;
}

int  
genotype_specific_sub_model::synchronize_two_dom()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(0);
  my_types[index_BB] = getParam(1);

  return 0;
}
 
int  
genotype_specific_sub_model::synchronize_two_rec()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(1);
  my_types[index_BB] = getParam(1);

  return 0;
}

int  
genotype_specific_sub_model::synchronize_three_add()
{
  my_types[index_AA] = getParam(0);
  my_types[index_BB] = getParam(2);

  // - For this option one parameter is a function of the other two.

  my_types[index_AB] = (my_types[index_AA] + my_types[index_BB]) / 2.0;
  
  // - Update maxfun w. newly calculated value of 'AB' parameter.
  //
  getParam(1) = my_types[index_AB];
  
  return 0;
}
 
//============================================================================
// IMPLEMENTATION:  genotype_specific_mean_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
bool
genotype_specific_mean_susc_sub_model::set
      (sm_option opt, const model_input& mean_AA, const model_input& mean_AB, 
       const model_input& mean_BB, primary_type pt)
{
  
  // - Parser revised 6-27-2 so that it should catch this situation before
  //   execution gets here.
  //
  if(! primary_trait_type_ok(pt))
  {
    my_errors << priority(critical) << "Primary trait type inconsistent "
              << "with specification of " << type_param()
              << " sub-block.  Skipping analysis ..." << endl;
    return false;
  }
  
  return set(opt, mean_AA, mean_AB, mean_BB);
}

// - Call the appropriate sub-model setting function for the 
//   specified sub-model option.  A return value of false means
//   that sub-model set failed.
//
bool
genotype_specific_mean_susc_sub_model::set(sm_option opt, const model_input& mean_AA,
                                      const model_input& mean_AB, const model_input& mean_BB)
{
  if(isLinked())
  {
    return false;
  }
  
  my_option = opt;

  bool return_value = false;
  switch(opt)
  {
    case one:
      my_type_count = 1;
      return_value = set_one(mean_AA, mean_AB, mean_BB);
      break;
    case two:
      my_type_count = 2;
      return_value = set_two(mean_AA, mean_AB, mean_BB);
      break;
    case three:
      my_type_count = 3;
      return_value = set_three(mean_AA, mean_AB, mean_BB);
      break;
    case two_dom:
      my_type_count = 2;
      return_value = set_two_dom(mean_AA, mean_AB, mean_BB);
      break;
    case two_rec:
      my_type_count = 2;
      return_value = set_two_rec(mean_AA, mean_AB, mean_BB);
      break;
    case three_add:
      my_type_count = 3;
      return_value = set_three_add(mean_AA, mean_AB, mean_BB);
      break;
    case three_dec:
      my_type_count = 3;
      return_value = set_three_dec(mean_AA, mean_AB, mean_BB);
      break;
    case three_inc:
      my_type_count = 3;
      return_value = set_three_inc(mean_AA, mean_AB, mean_BB);
      break;

    default:
      SAGE_internal_error();  
  }
  
  return return_value;
}

bool
genotype_specific_sub_model::set_three_dec
      (const model_input& mean_AA, const model_input& mean_AB, const model_input& mean_BB)
{
  bool  return_value = set_three(mean_AA, mean_AB, mean_BB, MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL);
     
  return_value = return_value &&
              constraint_met(&genotype_specific_sub_model::three_dec_constraint_met);
                           
  return return_value;
}

bool
genotype_specific_sub_model::set_three_inc
      (const model_input& mean_AA, const model_input& mean_AB, const model_input& mean_BB)
{
  bool  return_value =  set_three(mean_AA, mean_AB, mean_BB, MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL);
                           
  return_value = return_value && 
                 constraint_met(&genotype_specific_sub_model::three_inc_constraint_met);
  
  return return_value;
}      

//
// ===============  Ancillary Functions
//
bool
genotype_specific_sub_model::constraint_met(constraint_met_ptr constraint)
{
  if(! (this->*constraint)())
  {
    my_errors << priority(critical) << "Genotype specific " << long_plural_type()
              <<" violate constraints.  Skipping analysis ..." << endl;
    return false;
  }
  else
  {
    return true;
  }
}

// - Values of genotype specific means must be in decreasing order. 
//   Comparisons involving a QNAN always return true.
//
bool
genotype_specific_sub_model::three_dec_constraint_met() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_gte_AB = AA_or_AB_nan || my_types[0] >= my_types[1];
  bool  AB_gte_BB = AB_or_BB_nan || my_types[1] >= my_types[2];
  bool  AA_gte_BB = AA_or_BB_nan || my_types[0] >= my_types[2];

  if(AA_gte_AB && AB_gte_BB && AA_gte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

// - Values of genotype specific means must be in increasing order.
//   Comparisons involving a QNAN always return true.
//
bool
genotype_specific_sub_model::three_inc_constraint_met() const
{
  bool  AA_or_AB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[1]);
  bool  AB_or_BB_nan = SAGE::isnan(my_types[1]) || SAGE::isnan(my_types[2]);
  bool  AA_or_BB_nan = SAGE::isnan(my_types[0]) || SAGE::isnan(my_types[2]);
  
  bool  AA_lte_AB = AA_or_AB_nan || my_types[0] <= my_types[1];
  bool  AB_lte_BB = AB_or_BB_nan || my_types[1] <= my_types[2];
  bool  AA_lte_BB = AA_or_BB_nan || my_types[0] <= my_types[2];

  if(AA_lte_AB && AB_lte_BB && AA_lte_BB)
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool
genotype_specific_mean_susc_sub_model::primary_trait_type_ok(primary_type pt) const
{
  return (my_derived_type == mean && pt == pt_CONTINUOUS) ||
         (my_derived_type == susc && pt == pt_BINARY)     ||
          pt == pt_ONSET;
}


//
// ==============  Synchronization w. Maxfun
//
// - Call the appropriate sub-model synchronizing function for the 
//   current sub-model option. 
//
int
genotype_specific_mean_susc_sub_model::update() 
{
  int  return_value = 1;

  switch(my_option)
  {
    case one:
      return_value = synchronize_one();
      break;
    case two:
      return_value = synchronize_two();
      break;
    case three:
      return_value = synchronize_three();
      break;
    case two_dom:
      return_value = synchronize_two_dom();
      break;
    case two_rec:
      return_value = synchronize_two_rec();
      break;
    case three_add:
      return_value = synchronize_three_add();
      break;
    case three_dec:
      return_value = synchronize_three_dec();
      break;
    case three_inc:
      return_value = synchronize_three_inc();
      break;
    default:
      SAGE_internal_error();
  }
  
  return return_value;
}


// - If constraints are met, maxfun requires a return of 0.  If not, maxfunc
//   requires a return of 1. 
//
int  genotype_specific_sub_model::synchronize_three_dec()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(1);
  my_types[index_BB] = getParam(2);

  if(three_dec_constraint_met())
  {
    return 0;
  }
  else
  {
    // - sub-model will now be in an inconsistent state, but this will be 
    //   corrected on the next successful synchronization.  In the interim
    //   maxfun will not call evaluate() so we do not have a problem 
    //   (per discussion w. gcw).
    //
    return 1;    
  }
}
 
int  genotype_specific_sub_model::synchronize_three_inc()
{
  my_types[index_AA] = getParam(0);
  my_types[index_AB] = getParam(1);
  my_types[index_BB] = getParam(2);

  if(three_inc_constraint_met())
  {
    return 0;
  }
  else
  {
    return 1;   
  }
}

//============================================================================
// IMPLEMENTATION:  genotype_specific_susceptibility_sub_model
//============================================================================
//
// - All member functions that are distinct from mean sub_model are inline.


//============================================================================
// IMPLEMENTATION:  genotype_specific_variance_sub_model
//============================================================================
//
// ===============  Setting the sub-model
//
bool
genotype_specific_variance_sub_model::set
      (sm_option opt, const model_input& var_AA, const model_input& var_AB, 
       const model_input& var_BB, bool type_missing, bool trans_missing, primary_type pt)
{
  // - Parser revised 6-27-2 so that it should catch this situation before
  //   execution gets here.
  // 
  if(! primary_trait_type_ok(pt))
  {
    my_errors << priority(warning) << "Specification of type_var sub-block "
              << "inconsistent with primary trait type.  Ignoring ..." << endl;
    return true;
  }
  
  bool set_successful = set(opt, var_AA, var_AB, var_BB, type_missing, trans_missing);

  return set_successful;
}

// - Call the appropriate sub-model setting function for the 
//   specified sub-model option.  A return value of false means
//   that sub-model set failed.
//
bool
genotype_specific_variance_sub_model::set
      (sm_option opt, const model_input& vAA,const model_input& vAB, 
       const model_input& vBB, bool type_missing, bool trans_missing)
{
  if(isLinked())
  {
    return false;
  }

  set_defaults();
  
  // - Not specifying mean sub-model is a way of specifying commingling
  //   analysis, but certain constraints apply.  See also set functions
  //   in transmission and genotype frequency sub-models.
  //
  if(type_missing && opt != one)
  {
    my_errors << priority(error) << "type_var sub-block specified without "
              << "specifying a mean sub-block.  Ignoring type_var sub-block ..." << endl;
    return true;
  }
  
  // - Specifying mean sub-model and not specifying transmission model 
  //   is a way of specifying transmission analysis, but certain constraints 
  //   apply.  See also set function in the genotype frequency sub-model.
  //
  assert(my_m_ptr != 0);
  if(! type_missing && trans_missing && opt != one && my_m_ptr->option() != one)
  {
    my_errors << priority(error) << "type_var sub-block specified with "
              << "mean sub-block specified and transmission sub-block not "
              << "specified.  Ignoring type_var sub-block ..." << endl;
    return true;
  }
  
  if(mv_types(opt))
  {
    mv_types_message(opt, my_m_ptr->option());
    return false;
  }

  model_input var_AA = vAA;
  model_input var_AB = vAB;
  model_input var_BB = vBB;
  
  if(! (var_meets_constraints(var_AA) &&
        var_meets_constraints(var_AB) &&
        var_meets_constraints(var_BB)    ))
  {
    return false;
  }
  
  my_option = opt;

  bool  return_value = false;
  switch(opt)
  {  
    case one:
      my_type_count = 1;
      return_value = set_one(var_AA, var_AB, var_BB);
      break;
    case two:
      my_type_count = 2;
      return_value = set_two(var_AA, var_AB, var_BB);
      break;
    case three:
      my_type_count = 3;
      return_value = set_three(var_AA, var_AB, var_BB);
      break;
    case two_dom:
      my_type_count = 2;
      return_value = set_two_dom(var_AA, var_AB, var_BB);
      break;
    case two_rec:
      my_type_count = 2;
      return_value = set_two_rec(var_AA, var_AB, var_BB);
      break;
    case three_add:
      my_type_count = 3;
      return_value = set_three_add(var_AA, var_AB, var_BB);
      break;

    case three_dec:  // Should never happen!
    case three_inc: 
    default:
      SAGE_internal_error();
  }

  return return_value;
}

//
// ===============  Ancillary Functions
//

bool
genotype_specific_variance_sub_model::primary_trait_type_ok(primary_type pt) const
{
  return pt == pt_CONTINUOUS || pt == pt_ONSET;
}

// - Check to see that the type_var option is consistent w. the type_mean
//   option.
//
void
genotype_specific_variance_sub_model::mv_types_message
      (sm_option opt, sm_option mean_opt)
{
  //lint --e{1705}

  genotype_specific_sub_model* m = my_m_ptr;

  my_errors << priority(critical) << "Option, " << option_2_description(opt)
            << ", specified in type_var sub-block with option, "
            << m->static_option_description(mean_opt, m->my_derived_type) << ", specified "
            << "in " << type_param(m->my_derived_type) << " sub-block.  Skipping analysis ..." << endl;
}

// - Is there an inconsistency between mean and variance sub-models in terms of number
//   of types?
//
bool  
genotype_specific_variance_sub_model::mv_types(sm_option var_option) const
{
  bool types_inconsistent = true;
  
  assert(my_m_ptr != 0);
  genotype_specific_mean_sub_model::sm_option  mean_option = my_m_ptr->option();
  
  switch(mean_option)
  {
    case one:
        types_inconsistent = var_option != genotype_specific_variance_sub_model::one;
      break;
      
    case two:
        types_inconsistent = (var_option != genotype_specific_variance_sub_model::one  &&
                              var_option != genotype_specific_variance_sub_model::two);
      break;
      
    case two_dom:
        types_inconsistent = (var_option != genotype_specific_variance_sub_model::one      &&
                              var_option != genotype_specific_variance_sub_model::two_dom);
      break;
      
    case two_rec:
        types_inconsistent = (var_option != genotype_specific_variance_sub_model::one     &&
                              var_option != genotype_specific_variance_sub_model::two_rec);
      break;
      
    case three:
    case three_add:
    case three_dec:
    case three_inc:
      types_inconsistent = false;
      break;
      
    default:
      SAGE_internal_error();
  }
  
  return types_inconsistent; 
}

// - If a fixed input value is not greater than lb (or QNAN), return false.
//   If a non-fixed input value is not greater than lb (or QNAN), set value to
//   QNAN.  In either case, write an appropriate message.
//
bool
genotype_specific_variance_sub_model::var_meets_constraints(model_input& var)
{
  bool  return_value = false;
  if(SAGE::isnan(var.value))
  {
    return_value = true;
  }
  else
  {
    if(VAR_LB + VAR_EPSILON <= var.value)
    {
      return_value = true;
    }
    else
    {
      if(var.fixed == true)
      {
        my_errors << priority(critical) << "Non-positive value specified for parameter, var, "
                  << "in type_var sub-block.  Skipping analysis ..." << endl;
      }
      else
      {
        var.value = QNAN;
        my_errors << priority(error) << "Non-positive value specified for parameter, var, "
                  << "in type_var sub-block.  Ignoring ..." << endl;
        return_value = true;
      }
    }
  }
  
  return return_value;  
}

//
// ==============  Synchronization w. Maxfun
//
// - Call the appropriate sub-model synchronizing function for the 
//   current sub-model option. 
//
int
genotype_specific_variance_sub_model::update() 
{
  int  return_value = 1;

  switch(my_option)
  {
    case one:
      return_value = synchronize_one();
      break;
    case two:
      return_value = synchronize_two();
      break;
    case three:
      return_value = synchronize_three();
      break;
    case two_dom:
      return_value = synchronize_two_dom();
      break;
    case two_rec:
      return_value = synchronize_two_rec();
      break;
    case three_add:
      return_value = synchronize_three_add();
      break;

    case three_inc:  // Should never happen
    case three_dec:

    default:
      SAGE_internal_error();
  }
  
  return return_value;
}

}
}

