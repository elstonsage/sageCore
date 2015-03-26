//==========================================================================
//  File:       parser.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial implementation.                              Jul. 03
//
//  Notes:      This file implements a parser for reltest analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "reltest/parser.h"

namespace SAGE
{

namespace RELTEST
{

reltest_parser::reltest_parser()
{
  my_genome = NULL;
  my_analysis_name = "default_analysis";

  my_nucfam_out          = false;
  my_detailed_out        = false;
  my_calculate_cutpoints = true;

  my_analysis_pairtypes.resize(0);
  my_analysis_regions.resize(0);
  my_preset_cutpoints.resize(0);
}

reltest_parser::reltest_parser(const reltest_parser& rp)
{
  my_genome              = rp.my_genome;
  my_analysis_name       = rp.my_analysis_name;

  my_nucfam_out          = rp.my_nucfam_out;
  my_detailed_out        = rp.my_detailed_out;
  my_calculate_cutpoints = rp.my_calculate_cutpoints;

  my_analysis_pairtypes  = rp.my_analysis_pairtypes;
  my_analysis_regions    = rp.my_analysis_regions;
  my_preset_cutpoints    = rp.my_preset_cutpoints;
}

reltest_parser::~reltest_parser()
{}

void
reltest_parser::build_region_list()
{
  if( !my_genome )
    return;

  size_t region_count = my_genome->region_count();

  my_analysis_regions.resize(0);
  for( size_t r = 0; r < region_count; ++r )
  {
    if( my_genome->region(r).valid() && my_genome->region(r).locus_count() > 0 )
      my_analysis_regions.push_back(my_genome->region(r).index());
  }
}

void
reltest_parser::build_default_reltest_analysis(reltest_data& rd)
{
  my_genome = rd.genome();
  build_region_list();

  my_analysis_pairtypes.resize(5, true);
  my_preset_cutpoints.resize(5, std::numeric_limits<double>::quiet_NaN());
}

void
reltest_parser::parse_reltest_analysis(reltest_data&       rd,
                                       const LSFBase*      params,
                                       cerrorstream&       err)
{
  my_genome = rd.genome();
  errors = err;

  //handling analysis options
  if( !params || !params->List() )
  {
    build_region_list();
    return;
  }

  my_analysis_pairtypes.resize(0);
  my_analysis_regions.resize(0);
  my_preset_cutpoints.resize(0);

  my_analysis_pairtypes.resize(5, false);
  my_preset_cutpoints.resize(5, std::numeric_limits<double>::quiet_NaN());

  LSFList::const_iterator i;
  AttrVal a;

  for( i = params->List()->begin(); i != params->List()->end(); ++i )
  {
    const LSFBase* param = *i;

    if( !param || !(param->name().size()) ) continue;

    if( param->attrs() == NULL || param->attrs()->empty() ) continue;

    string name = toUpper( param->name() );

    //check if there is any prefered region.
    if( name == "REGION" )
    {
      parse_region(param);
    }
    else if( name == "PAIR_TYPE" )
    {
      parse_pairtype(param);
    }
    else if( name == "CUT_POINTS" || name == "CUTPOINTS" )
    {
      parse_cutpoints(param);
    }
    else if( name == "NUCFAM" || name == "NUCFAM_FILE" )
    {
      parse_nucfam(param);
    }
    else if( name == "DETAILED" || name == "DETAILED_FILE" )
    {
      parse_detailed(param);
    }
  }

  if( !my_analysis_regions.size() )
    build_region_list();

  bool pairtype_specified = false;

  my_calculate_cutpoints = false;
  for( size_t i = 0; i < my_preset_cutpoints.size(); ++i )
  {
    if( SAGE::isnan(my_preset_cutpoints[i]) )
      my_calculate_cutpoints = true;

    if( my_analysis_pairtypes[i] )
      pairtype_specified = true;
  }

  if( !pairtype_specified )
  {
    my_analysis_pairtypes.resize(0);
    my_analysis_pairtypes.resize(5, true);
  }
}

void
reltest_parser::parse_region(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    string rname = a.String();

    if( !my_genome->region(rname).valid() || my_genome->region(rname).locus_count() == 0 )
      errors << priority(warning) << "Region name '" << rname
                << "' is invalid.  It will be ignored."  << endl;
    else
      my_analysis_regions.push_back(my_genome->region(rname).index());
  }
  else
  {
    errors << priority(information)
           << "No value for 'region' is specified.  Skipping..."
           << endl;
  }
}

void
reltest_parser::parse_pairtype(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());

    if( n == "SIB" )
      my_analysis_pairtypes[SIB] = true;
    else if( n == "HSIB" )
      my_analysis_pairtypes[HSIB] = true;
    else if( n == "MZTWIN" )
      my_analysis_pairtypes[MZTWIN] = false; //if supported, turn it to 'true'
    else if( n == "PARENT_OFFSPRING" || n=="PO" )
      my_analysis_pairtypes[PARENT_OFFSPRING] = true;
    else if( n == "MARITAL" )
      my_analysis_pairtypes[MARITAL] = true;
    else
      errors << priority(warning) << "Pair type '" << n
             << "' unrecognized.  It will be ignored." << endl;
  }
  else
  {
    errors << priority(information)
           << "No value for 'pair_type' is specified.  Skipping..."
           << endl;
  }
}

void
reltest_parser::parse_cutpoints(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrList::const_iterator a;

  if( (a = param->attrs()->find("UNRELATED")) != param->attrs()->end() )
  { 
    if( finite(a->second.Real()) )
    {
      my_preset_cutpoints[Cu] = a->second.Real();
    }
    else
      errors << priority(error)
             << "Unrelated cutpoint is invalid. "
             << "Cutpoint will be ignored." << endl;
  }

  if( (a = param->attrs()->find("HSIB")) != param->attrs()->end() )
  { 
    if( finite(a->second.Real()) )
    {
      my_preset_cutpoints[Ch] = a->second.Real();
    }
    else
      errors << priority(error)
             << "Half-sibling cutpoint is invalid. "
             << "Cutpoint will be ignored." << endl;
  }

  if( (a = param->attrs()->find("MZTWIN")) != param->attrs()->end() )
  { 
    if( finite(a->second.Real()) )
    {
      my_preset_cutpoints[Cm] = a->second.Real();
    }
    else
      errors << priority(error)
             << "Monozygotic twin cutpoint is invalid. "
             << "Cutpoint will be ignored." << endl;
  }

  if( (a = param->attrs()->find("PARENT_OFFSPRING")) != param->attrs()->end() )
  { 
    if( finite(a->second.Real()) )
    {
      my_preset_cutpoints[Cp] = a->second.Real();
    }
    else
      errors << priority(error)
             << "Parent/offspring cutpoint is invalid. "
             << "Cutpoint will be ignored." << endl;
  }
}

void
reltest_parser::parse_nucfam(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());

    if( n=="YES" || n == "TRUE" )
      my_nucfam_out = true;
    else if( n=="NO" || n == "FALSE" )
      my_nucfam_out = false;
    else
      errors << priority(warning)
             << "Invalid value for 'nucfam'. "
             << "Nuclear family file won't be generated." << endl;
  }
  else
  {
    errors << priority(information)
           << "No value for 'nucfam' is specified.  Skipping..."
           << endl;
  }
}

void
reltest_parser::parse_detailed(const LSFBase* param)
{
  if( !param || !param->attrs() )
    return;

  AttrVal a = attr_value(param, 0);

  if( a.has_value() )
  {
    string n = toUpper(a.String());

    if( n=="YES" || n == "TRUE" )
      my_detailed_out = true;
    else if( n=="NO" || n == "FALSE" )
      my_detailed_out = false;
    else
      errors << priority(warning)
             << "Invalid value for 'detailed'. "
             << "Detailed output file won't be generated." << endl;
  }
  else
  {
    errors << priority(information)
           << "No value for 'detailed' is specified.  Skipping..."
           << endl;
  }
}

void
reltest_parser::view_parameter() const
{
  cout << endl;
  cout << "===================================================" << endl;

  cout << "  Analysis name : " << get_analysis_name() << endl;
  cout << "  Pair Types    :";

  if( get_analysis_pairtypes()[SIB] )
    cout << " SIB";
  if( get_analysis_pairtypes()[HSIB] )
    cout << " HSIB";
  if( get_analysis_pairtypes()[MZTWIN] )
    cout << " MZTWIN";
  if( get_analysis_pairtypes()[PARENT_OFFSPRING] )
    cout << " PARENT_OFFSPRING";
  if( get_analysis_pairtypes()[MARITAL] )
    cout << " MARITAL";
  cout << endl << endl;

  cout << "  Preset Cutpoints :" << endl;
  
  cout << "            UNRELATED : ";
  if( SAGE::isnan(get_preset_cutpoints()[Cu]) )
    cout << "******";
  else
    cout << fp(get_preset_cutpoints()[Cu],5,2);
  cout << endl;

  cout << "            HSIB      : ";
  if( SAGE::isnan(get_preset_cutpoints()[Ch]) )
    cout << "******";
  else
    cout << fp(get_preset_cutpoints()[Ch],5,2);
  cout << endl;

  cout << "            MZTWIN    : ";
  if( SAGE::isnan(get_preset_cutpoints()[Cm]) )
    cout << "******";
  else
    cout << fp(get_preset_cutpoints()[Cm],5,2);
  cout << endl;

  cout << "            PARENT_OFFSPRING : ";
  if( SAGE::isnan(get_preset_cutpoints()[Cp]) )
    cout << "******";
  else
    cout << fp(get_preset_cutpoints()[Cp],5,2);
  cout << endl;

  cout << endl;

  if( calculate_cutpoints() )
    cout << "  Calculate cutpoints!" << endl;

  if( generate_nucfam_output() )
    cout << "  Generate nucfam output!" << endl;

  if( generate_detailed_output() )
    cout << "  Generate detailed output!" << endl;

  cout << endl;
  cout << "  Region Count : " << get_analysis_regions().size();

  for( size_t r = 0; r < get_analysis_regions().size(); ++r )
  {
    if( !(r % 5) )
      cout << endl << "         ";

    cout << "  " << my_genome->region_name(get_analysis_regions()[r]);
  }
  cout << endl;

  cout << endl;
  cout << "===================================================" << endl;
  cout << endl;
}

} // end of namespace RELTEST

} // end of namespace SAGE
