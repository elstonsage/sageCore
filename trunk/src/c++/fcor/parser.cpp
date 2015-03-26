//****************************************************************************
//* File:      parser.cpp                                                    *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                    1.0  yjs  Added parsing for var-cov option  Jan 23 01 *
//*                                                                          *
//* Notes:     This source file defines fuctions to parse fcor parameters.   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/parser.h"

using namespace std;

namespace SAGE {
namespace FCOR {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     FcorParser                                                  ~
// ~                                                                        ~
// ~ Purpose:   Parse fcor parameters.                                      ~
// ~                                                                        ~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                                              
//---------------------------------------------------------------------------
// Out-of-Line Implementation of FcorParser                                  
//---------------------------------------------------------------------------

FcorParser::FcorParser(const RPED::RefMultiPedigree& mp, cerrorstream& err)
{
  errors = err;

  my_multipedigree = &mp;
  my_trait_list.resize(0);
}

FcorParser::FcorParser(const FcorParser& f)
{
  my_multipedigree         = f.my_multipedigree;
  my_trait_list            = f.my_trait_list;
  my_analysis_options      = f.my_analysis_options;
  my_relationship_name_map = f.my_relationship_name_map;
  errors                   = f.errors;
}

void
FcorParser::build_default_fcor_analysis()
{
  build_trait_list();
}

void
FcorParser::parse_fcor_analysis(const LSFBase* fcor_params)
{
  if( !fcor_params || !fcor_params->List() )
  {
    build_trait_list();
    return;
  }

  my_trait_list.resize(0);

  LSFList::const_iterator i;
  AttrVal a;

  for( i = fcor_params->List()->begin(); i != fcor_params->List()->end(); ++i)
  {
    if( !*i ) continue;

    a = attr_value(*i, "TRAIT", 0);
    if( a.has_value() )
    {
      size_t trait_index = my_multipedigree->info().trait_find(a.String());

      if( trait_index == (size_t) - 1 )
        errors << priority(information)
               << "Invalid trait name for fcor_analysis : " << a.String()
               << "\n             Skipped... " << endl;
      else
        my_trait_list.push_back(make_pair(a.String(), trait_index));
    }
  }
  
  if( !my_trait_list.size() )
    build_trait_list();

  for( i = fcor_params->List()->begin(); i != fcor_params->List()->end(); ++i )
  {
    const LSFBase* param = *i;

    if( !param || !(param->name().size()) ) continue;

    string name = toUpper( param->name() );

    if( name == "CLASS_WEIGHT" || name == "WEIGHT" )
      parse_class_weight(attr_value(param, 0), my_analysis_options.class_weight);
    else if( name == "SEX_NAME" || name == "GENDER_NAME" || name == "VERBOSE_NAME" )
      parse_boolean_value(attr_value(param, 0), my_analysis_options.gender_name);
    else if( name == "STANDARD_ERROR" )
      parse_boolean_value(attr_value(param, 0), my_analysis_options.standard_error);
    else if( name == "CONSERVATIVE" )
      parse_boolean_value(attr_value(param, 0), my_analysis_options.conservative);
    else if( name == "GENERATION_LIMIT" )
      parse_generation_limit(attr_value(param, "GENERATION_LIMIT", 0));
    else if( name == "CORRELATIONS" || name == "TYPE" )
      parse_pairset(param);
    else if( name == "HOMOGENEITY_TEST" )
      parse_homogeneity_test(param);
    else if( name == "VAR_COV" )
      parse_var_cov(param);
    else if( name == "OUTPUT_OPTIONS" || name == "OUTPUT_OPTION" || name == "OUTPUT" )
      parse_output_options(param);
    else if( name == "KE" )
      parse_boolean_value(attr_value(param, 0), my_analysis_options.KE_with_optimal);
    else if( name == "TRAIT" )
    {}
    else
      errors << priority(information)
             << "Invalid parameter for fcor_analysis : " << param->name()
             << "\n             Skipped..." << endl;
  }

  if(     my_analysis_options.class_weight == WEIGHT_COUNT
      && !my_analysis_options.standard_error )
  {
    //errors << priority(information)
    //       << "Standard errors need to be calculated to find optimal weights. "
    //       << " Setting it to true..." << endl;
    my_analysis_options.standard_error = true;
  }

  if(   !(my_analysis_options.class_weight == WEIGHT_COUNT || get_trait_count() > 1)
      &&  my_analysis_options.detailed_output )
  {
    errors << priority(information)
           << "No detailed results available. "
           << " Setting it to false..." << endl;
    my_analysis_options.detailed_output = false;
  }
}

void
FcorParser::build_trait_list()
{
  if( !my_multipedigree )
    return;

  size_t trait_count = my_multipedigree->info().trait_count();

  my_trait_list.resize(0);
  for( size_t t = 0; t < trait_count; ++t )
  {
    if( my_multipedigree->info().trait_info(t).type() == RPED::RefTraitInfo::binary_trait )
      continue;

    name_index_type a_trait;
    a_trait.first  = toUpper(my_multipedigree->info().trait_info(t).name());
    a_trait.second = t;

    my_trait_list.push_back(a_trait);
  }
}

void
FcorParser::parse_class_weight(const AttrVal& a, weight_type& class_weight)
{
  if( a.has_value() )
  {
    if(         toUpper(a.String()) == "PAIR_WISE"
             || toUpper(a.String()) == "PAIR"      )
      class_weight = PAIR_WISE;

    else if(    toUpper(a.String()) == "UNIFORM" )
      class_weight = UNIFORM;

    else
    {
      class_weight = WEIGHT_COUNT;
      errors << priority(information)
             << "Invalid value for weight : " << a.String()  
             << ".  The default is used." << endl;
    }
  }
  else
  {
    class_weight = WEIGHT_COUNT;
    errors << priority(information)
           << "No value for weight : " << a.String()  
           << ".  The default is used." << endl;
  }
}

void
FcorParser::parse_generation_limit(const AttrVal& a)
{
  if( a.has_value() )
  {
    if( !finite(a.Real()) || a.Int() < 0 )
      errors << priority(information)
             << "Invalid value for generation_limit : " << a.String()
             << "\n            Skipped..." << endl;
    else
      my_analysis_options.generation_limit = a.Int();
  }
}

void
FcorParser::parse_pairset(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    if(         toUpper(a.String()) == "SUBTYPES"
             || toUpper(a.String()) == "SUB"      )
      my_analysis_options.pairset = SUBTYPES;

    else if(    toUpper(a.String()) == "MAINTYPES"
             || toUpper(a.String()) == "MAIN"       )
      my_analysis_options.pairset = MAINTYPES;

    else if(    toUpper(a.String()) == "BOTH" )
      my_analysis_options.pairset = BOTH;

    else
    {
      my_analysis_options.pairset = SUBTYPES;
      errors << priority(information)
             << "Invalid value for correlations : " << a.String()  
             << "\n             The default subtypes is used." << endl;
    }
  }
  else
  {
    my_analysis_options.pairset = SUBTYPES;
    errors << priority(information)
           << "No value for correlations is specified."
           << "\n             The default subtypes is used." << endl;
  }
}

void
FcorParser::parse_homogeneity_test(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  if( param->attrs()->has_attr("valid_only") )
    my_analysis_options.valid_only_homog = true;

  if(    param->attrs()->has_attr("ind")
      || param->attrs()->has_attr("individual")
      || param->attrs()->has_attr("individual_trait")
      || param->attrs()->has_attr("single")
      || param->attrs()->has_attr("single_trait")
      || param->attrs()->has_attr("each")
      || param->attrs()->has_attr("each_trait") )
    my_analysis_options.individual_homog = true;

  parse_boolean_value(attr_value(param, 0), my_analysis_options.homogeneity_test);
}

void
FcorParser::parse_var_cov(const LSFBase* vc_params)
{
  if( !vc_params || !vc_params->List() )
    return;

  var_cov_param::matrix_type m;
  if( vc_params->attrs() && vc_params->attrs()->has_attr("joint") )
    m = var_cov_param::JOINT;
  else
    m = var_cov_param::SINGLE;

  vector<name_index_type> vc_traits;
  vc_traits.resize(0);

  vector<name_index_type> sorted_vc_traits;
  sorted_vc_traits.resize(0);

  LSFList::const_iterator i;
  AttrVal a;

  for( i = vc_params->List()->begin(); i != vc_params->List()->end(); ++i)
  {
    if( !*i ) continue;

    a = attr_value(*i, "TRAIT", 0);
    if( a.has_value() )
    {
      size_t trait_index = my_multipedigree->info().trait_find(a.String());

      bool trait_exist = false;
      for( size_t t = 0; t < get_trait_list().size(); ++t )
        if( trait_index == get_trait_list()[t].second )
          trait_exist = true;

      if( trait_exist )
        vc_traits.push_back(make_pair(a.String(), trait_index));
      else
        errors << priority(information) << "Invalid trait name for var_cov block : "
               << a.String() << "\n             Skipped..." << endl;
    }
  }

  if( !vc_traits.size() )
    sorted_vc_traits = my_trait_list;
  else
  {
    for( size_t anal_t = 0; anal_t < get_trait_list().size(); ++anal_t )
      for( size_t t = 0; t < vc_traits.size(); ++t )
        if( get_trait_list()[anal_t].second == vc_traits[t].second )      
        {
          sorted_vc_traits.push_back(get_trait_list()[anal_t]);
          break;
        }
  }

  vector<string> co_name;
  for( i = vc_params->List()->begin(); i != vc_params->List()->end(); ++i )
  {
    if( !*i ) continue;

    if( toUpper((*i)->name() ) == "CORRELATION" )
    {
      a = attr_value(*i, "CORRELATION", 0);
      if( a.has_value() )
        co_name.push_back(toUpper(a.String()));
    }
    else if( toUpper((*i)->name() ) == "TRAIT" )
    {}
    else
      errors << priority(information) << "Invalid parameter for var_cov block : "
             << (*i)->name() << "\n            Skipped..." << endl;
  }
  if( !co_name.size() )
  {
    errors << priority(information) << "No correlation type(s) in var_cov block. "
           << "\n             Skipped..." << endl;
    return;
  }
  else if( co_name.size() == 1 )
  {
    name_name_type cn = make_pair(co_name[0], co_name[0]);
    my_analysis_options.var_covs.push_back(var_cov_param(m, cn, sorted_vc_traits));
  }
  else
  {
    name_name_type cn = make_pair(co_name[0], co_name[1]);
    my_analysis_options.var_covs.push_back(var_cov_param(m, cn, sorted_vc_traits)); 
  }
}

void
FcorParser::parse_output_options(const LSFBase* params)
{
  if( !params || !params->List() )
    return;

  LSFList::const_iterator i;
  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    if( !*i  || !(*i)->name().size() )
      continue;

    string name = toUpper( (*i)->name() );

    if( name == "DETAILED_OUT" || name == "DETAILED" )
      parse_boolean_value(attr_value(*i, 0), my_analysis_options.detailed_output);
    else if( name == "TABULAR_OUT" || name == "TABULAR" )
      parse_boolean_value(attr_value(*i, 0), my_analysis_options.xls_output);
    else if( name == "PAIRS_OUT" || name == "PAIR_OUT" || name == "PAIRS" || name == "PAIR" )
      parse_boolean_value(attr_value(*i, 0), my_analysis_options.pair_output);
    else if( name == "SEX_NAME" || name == "GENDER_NAME" || name == "VERBOSE_NAME" )
      parse_boolean_value(attr_value(*i, 0), my_analysis_options.gender_name);
    else
      errors << priority(information) << "Invalid parameter for output_options block : "
             << (*i)->name() << "\n            Skipped..." << endl;
  }
}

void
FcorParser::parse_boolean_value(const AttrVal& a, bool& value)
{
  if( a.has_value() )
  {
    if( toUpper(a.String()) == "YES" || toUpper(a.String()) == "TRUE" )
      value = true;
    else
      value = false;
  }
  else
    value = false;

  return;
}

void
FcorParser::view_parameter() const
{
  cout << endl;
  cout << "===================================================" << endl;

  if( my_analysis_options.pairset == MAINTYPES )
    cout << "# Main correlation within generation ";
  else if( my_analysis_options.pairset == BOTH )
    cout << "# Sub & Main correlation within generation ";
  else
    cout << "# Sub correlation within generation ";
  cout << my_analysis_options.generation_limit << endl;

  cout << "# Trait to be used for analysis :" << endl;
  for( size_t t = 0; t < my_trait_list.size(); ++t )
    cout << "     Trait " << t << " : " << my_trait_list[t].first << endl;

  cout << "# Class weight : " << get_weight(false) << endl;

  if(    my_analysis_options.standard_error
      && my_analysis_options.conservative )
    cout << "# Conservative calculation of standard errors" << endl;
  else if(     my_analysis_options.standard_error
           && !my_analysis_options.conservative )
    cout << "# Robust calculation of standard errors" << endl;

  if( my_analysis_options.gender_name )
    cout << "# Long name" << endl;
  else
    cout << "# Short name" << endl;

  if( my_analysis_options.homogeneity_test )
    cout << "# Homogeneity test within generation " << my_analysis_options.generation_limit << endl;

  for( size_t v = 0; v < my_analysis_options.var_covs.size(); ++v )
    cout << "# " << my_analysis_options.var_covs[v].name();

  cout << "===================================================" << endl;
}

void
FcorParser::build_relationship_name_map()
{
  my_relationship_name_map["0"]     = "SELF";
  my_relationship_name_map["1"]     = "MOTHER:FATHER";
  my_relationship_name_map["10"]    = "PARENT:OFFSPRING";
  my_relationship_name_map["11"]    = "SIBLING";
  my_relationship_name_map["11h"]   = "HALF-SIBLING";
  my_relationship_name_map["20"]    = "GRANDPARENTAL";
  my_relationship_name_map["21"]    = "AVUNCULAR";
  my_relationship_name_map["22"]    = "COUSIN";
  my_relationship_name_map["M"]     = toUpper("male-self");
  my_relationship_name_map["F"]     = toUpper("female-self");
  my_relationship_name_map["M,F"]   = toUpper("mother:father");
  my_relationship_name_map["MM"]    = toUpper("father:son");
  my_relationship_name_map["FM"]    = toUpper("mother:son");
  my_relationship_name_map["MF"]    = toUpper("father:daughter");
  my_relationship_name_map["FF"]    = toUpper("mother:daughter");
  my_relationship_name_map["M,M"]   = toUpper("brother:brother");
  my_relationship_name_map["F,M"]   = toUpper("sister:brother");
  my_relationship_name_map["F,F"]   = toUpper("sister:sister");
  my_relationship_name_map["M,M,M"] = toUpper("paternal-half-brother:half-brother");
  my_relationship_name_map["F,M,M"] = toUpper("paternal-half-sister:half-brother");
  my_relationship_name_map["F,M,F"] = toUpper("paternal-half-sister:half-sister");
  my_relationship_name_map["M,F,M"] = toUpper("maternal-half-brother:half-brother");
  my_relationship_name_map["F,F,M"] = toUpper("maternal-half-sister:half-brother");
  my_relationship_name_map["F,F,F"] = toUpper("maternal-half-sister:half-sister");
  my_relationship_name_map["MMM"]   = toUpper("grandfather-through-father:grandson");
  my_relationship_name_map["FMM"]   = toUpper("grandmother-through-father:grandson");
  my_relationship_name_map["MFM"]   = toUpper("grandfather-through-mother:grandson");
  my_relationship_name_map["FFM"]   = toUpper("grandmother-through-mother:grandson");
  my_relationship_name_map["MMF"]   = toUpper("grandfather-through-father:granddaughter");
  my_relationship_name_map["FMF"]   = toUpper("grandmother-through-father:granddaughter");
  my_relationship_name_map["MFF"]   = toUpper("grandfather-through-mother:granddaughter");
  my_relationship_name_map["FFF"]   = toUpper("grandmother-through-mother:granddaughter");
  my_relationship_name_map["M,MM"]  = toUpper("uncle-through-father:nephew");
  my_relationship_name_map["F,MM"]  = toUpper("aunt-through-father:nephew");
  my_relationship_name_map["M,FM"]  = toUpper("uncle-through-mother:nephew");
  my_relationship_name_map["F,FM"]  = toUpper("aunt-through-mother:nephew");
  my_relationship_name_map["M,MF"]  = toUpper("uncle-through-father:niece");
  my_relationship_name_map["F,MF"]  = toUpper("aunt-through-father:niece");
  my_relationship_name_map["M,FF"]  = toUpper("uncle-through-mother:niece");
  my_relationship_name_map["F,FF"]  = toUpper("aunt-through-mother:niece");
  my_relationship_name_map["MM,MM"] = toUpper("male-cousin-through-father:male-cousin-through-father");
  my_relationship_name_map["MF,MM"] = toUpper("male-cousin-through-mother:male-cousin-through-father");
  my_relationship_name_map["MF,FM"] = toUpper("male-cousin-through-mother:male-cousin-through-mother");
  my_relationship_name_map["FM,MM"] = toUpper("female-cousin-through-father:male-cousin-through-father");
  my_relationship_name_map["FM,FM"] = toUpper("female-cousin-through-father:male-cousin-through-mother");
  my_relationship_name_map["FF,MM"] = toUpper("female-cousin-through-mother:male-cousin-through-father");
  my_relationship_name_map["FF,FM"] = toUpper("female-cousin-through-mother:male-cousin-through-mother");
  my_relationship_name_map["FM,MF"] = toUpper("female-cousin-through-father:female-cousin-through-father");
  my_relationship_name_map["FF,MF"] = toUpper("female-cousin-through-mother:female-cousin-through-father");
  my_relationship_name_map["FF,FF"] = toUpper("female-cousin-through-mother:female-cousin-through-mother");
}

// end of FcorParser Implementation

} // end of namespace FCOR
} // end of namespace SAGE
