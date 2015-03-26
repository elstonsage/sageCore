//============================================================================
// File:      parser.cpp                      
//                                                                          
// Author:    Dan Baechle, Geoff Wedig                                    
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "boost/bind.hpp"
#include "mlod/data.h"
#include "mlod/parser.h"

using namespace std;

namespace SAGE
{

namespace MLOD
{

//============================================================================
// IMPLEMENTATION:  parser
//============================================================================
//
/// Sole purpose of this function is to provide a reasonable alternative to
/// the function name, 'parse_test_parameter_section'.
///
void
Parser::parse(const LSFBase* params)
{
  parse_test_parameter_section(params);
}

///  Iterate through parameters of mlod block.  Client code must insure that 
///  params->name() is 'mlod'or 'mlod_analysis'.
///
void
Parser::parse_test_parameter_section(const LSFBase* params)
{
  //lint --e{613}  Params will never be NULL during analysis thanks to the
  //               assertion
  
  assert(params != 0);
  string  params_name = toUpper(params->name());
  assert(params_name == "MLOD" || params_name == "MLOD_ANALYSIS");
  
  // Setup for parsing an alaysis.
  
  reset_parameters();
  print_header();
  
  // Make sure there are options specified.  If not, then it can't be a valid
  // analysis.
  //lint --e{666}  Side effects of List() are not an issue
  if(params->List() == 0)
  {
      errors << priority(error) << "No parameters detected for MLOD analysis. "
             << "Starting with version 5.1, MLOD analyses "
             << "require parameters to be specified in a block like other S.A.G.E. "
             << "programs.  Please see MLOD documentation for details. "
             << "Skipping analysis ..." << endl;
  }
  else
  {
    // Parse the ,out attribute if it exists
    parse_out(params, params_name);
    
    // Iterate through the block and parse each parameter
    
    LSFList::const_iterator  iter;
    for(iter = params->List()->begin(); iter != params->List()->end(); ++iter)
    {
      string  param_name = toUpper((*iter)->name());
      if(param_name == "TITLE")
      {
        string s;
        
        //lint -e{534} Ignoring return type.
        parse_string(*iter, s);
        
        my_parameters.set_title(s);
      }
      else if(param_name == "TRAIT_MARKER")
      {
        parse_trait_marker(*iter);
      }
      else if(param_name == "REGION")
      {
        parse_region(*iter);
      }
      else if(param_name == "MAX_SIZE")
      {
        parse_max_ped_size(*iter);
      }
      else if(param_name == "SCAN_TYPE")
      {
        parse_scan_type(*iter);
      }
      else if(param_name == "DISTANCE")
      {
        parse_distance(*iter);    
      }
      else if(param_name == "OUTPUT_PEDIGREES")
      {
        parse_output_pedigrees_option(*iter);
      }
      else if(param_name == "SAMPLE_DETAIL")
      {
        parse_ind_sample_table_option(*iter);
      }
      else
      {
        produce_bad_parameter_warning(param_name);    
      }
    }
    
    validate();
  }

  print_footer();
}

/// Parses the ',out' attribute of the mlod analysis block
///
void
Parser::parse_out(const LSFBase* params, const string& )
{
  //lint --e{613}  Params will never be NULL
  if(params->attrs())
  {
      string s;
      
      //lint -e{534} Ignoring return type.
      parse_string(params, "OUT", s);

      my_parameters.set_base_output_filename(s);
  }
}

/// Parses the 'TRAIT_MARKER' parameter in the mlod analysis block
///
void  
Parser::parse_trait_marker(const LSFBase* param)
{
  string  trait_marker_name;
  //lint -e{534} Ignoring return type.
  parse_string(param, trait_marker_name);
  if(trait_marker_name.size())
  {
    const MLOCUS::inheritance_model_map& trt_imm = my_data.get_trait_loci();
    const MLOCUS::inheritance_model_map& all_imm = my_data.markers();
    
    size_t  trait_marker_id = trt_imm.index(trait_marker_name);
    if(trait_marker_id != (size_t)(-1))
    {
      // Get the trait marker from the global set.  We use the global index to look
      // up our traits in general.
      size_t global_id = all_imm.index(trait_marker_name);
      
      if(is_trait_marker_new(global_id))
      {
        my_parameters.get_trait_list().push_back
            (std::make_pair(global_id, all_imm[global_id]));
      }
      else
      {
        errors << priority(warning) << "Repetition of trait marker " << trait_marker_name
               << " ignored." << endl;
      }
    }
    else
    {
      errors << priority(error) << "'" << trait_marker_name << "' not a valid trait marker.  "
             << "Ignoring ..." << endl;
    }
  }
}

/// Parses the 'REGION parameter in the mlod analysis block
///
void  
Parser::parse_region(const LSFBase* param)
{
  string  region_name;
  //lint -e{534} Ignoring return type.
  parse_string(param, region_name);

  if(region_name.size())
  {
    // - genome_description interface lacks some needed const accessors.  Quick
    //   attempt to add them produced a 'const cascade.'
    //
    RPED::genome_description* gd = const_cast<RPED::genome_description*>(my_data.genome());
    assert(gd != 0);
    //lint --e{613}  --e{666} bd will never be NULL; built doesn't have side effects to worry about
    assert(gd->built());
    
    if(gd->region(region_name).index() != (size_t)(-1)) 
    {
      my_parameters.set_region(gd->region(region_name));
    }
    else
    {
      errors << priority(error) << "'" << region_name << "' not a valid region.  "
             << "Ignoring ..." << endl;
    }
  }
}
  
/// Parses the 'MAX_SIZE' parameter in the mlod analysis block
///
void  
Parser::parse_max_ped_size(const LSFBase* param)
{
  int  value(-1); 
  
  size_t current_value = my_parameters.get_max_ped_size();
  
  //lint -e{534} Ignoring return type.
  parse_integer(param, value);
  
  if(value > 0)
  {
    if(2 <= value && value <= MLOD_MAX_MAX_PED_SIZE)
    {
      my_parameters.set_max_ped_size((size_t) value);
    }
    else if(value <= 2)
    {
      errors << priority(error) << "Maximum 2n-f pedigree size must be "
             << "at least 2.  Value of " << value << " will be ignored.  Current value of "
             << current_value << " will be used instead." << endl;
    }
    else
    {
      errors << priority(error) << "Maximum 2n-f pedigree size must be "
             << " less than or equal to " << MLOD_MAX_MAX_PED_SIZE << ".  Value of " 
             << value << " will be ignored.  Current value of "
             << current_value << " will be used instead." << endl;
    }
  }
}

/// Parses the 'SCAN_TYPE' parameter in the mlod analysis block
///
void  
Parser::parse_scan_type(const LSFBase* param)
{
  string  value;
  //lint -e{534} Ignoring return type.
  parse_string(param, value);

  if(value.size())
  {
    string cap_value = toUpper(value);

    if(cap_value == "MARKER" || cap_value == "MARKERS")
    {
      my_parameters.set_scan_type(AnalysisParameters::ST_MARKER);
    }
    else if(cap_value == "INTERVAL" || cap_value == "INTERVALS")
    {
      my_parameters.set_scan_type(AnalysisParameters::ST_INTERVAL);
    }
    else if(cap_value == "ALL" || cap_value == "BOTH")
    {
      my_parameters.set_scan_type(AnalysisParameters::ST_BOTH);
    }
    else
    {
      errors << priority(error) << "Value for parameter, scan_type, "
             << "not recognized.  Ignoring ..." << endl;
    }

    if( my_parameters.get_scan_type() != AnalysisParameters::ST_MARKER && param->attrs() )
    {
      AttrVal v = attr_value(param,"DISTANCE");

      if( !v.has_value() || v.String().empty() )
      {
        v = attr_value(param,"INTERVAL_DISTANCE");
      }

      if( !v.has_value() || v.String().empty() )
      {
        v = attr_value(param,"SCAN_DISTANCE");
      }

      if( v.has_value() && !v.String().empty() )
      {
        check_distance_value(v.Real());
      }
    }
  }
}
 
/// Parses the 'DISTANCE' parameter in the mlod analysis block
///
void  
Parser::parse_distance(const LSFBase* param)
{
  double  value = QNAN; 

  //lint -e{534} Ignoring return type.
  parse_real(param, value);

  check_distance_value(value);
}

void  
Parser::check_distance_value(double value)
{
  if(! SAGE::isnan(value))
  {
    if(DISTANCE_WARNING_LIMIT < value && value <= DISTANCE_UPPER_LIMIT)
    {
      my_parameters.set_scan_distance(value);
    }
    else if(DISTANCE_LOWER_LIMIT < value && value <= DISTANCE_WARNING_LIMIT)
    {
      errors << priority(warning) << "Specified distance '" << value 
             << "' is less than or equal to " << DISTANCE_WARNING_LIMIT 
             << " centimorgans.  Program run time may be unacceptably long." << endl;

      my_parameters.set_scan_distance(value);
    }
    else if(value <= DISTANCE_LOWER_LIMIT)
    {
      errors << priority(error) << "Specified distance '" << value 
             << "' is less than or equal to " << DISTANCE_LOWER_LIMIT 
             << " centimorgans.  Ignoring ..." << endl;
    }
    else
    {
      errors << priority(error) << "Specified distance '" << value 
             << "' is greater than  " << DISTANCE_UPPER_LIMIT 
             << " centimorgans.  Ignoring ..." << endl;
    }
  }
}

/// Parses the 'OUTPUT_PEDIGREES' parameter in the mlod analysis block
///
void
Parser::parse_output_pedigrees_option(const LSFBase* param)
{
  string  value;
  //lint -e{534} Ignoring return type.
  parse_string(param, value);

  if(value.size())
  {
    string cap_value = toUpper(value);
    if(cap_value == "NONE")
    {
      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_NONE);
    }
    else if(cap_value == "MARKERS" || cap_value == "MARKER")
    {
      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_MARKERS);
    }
    else if(cap_value == "INTERVALS" || cap_value == "INTERVAL")
    {
      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_INTERVALS);
    }
    else if(cap_value == "ALL" || cap_value == "BOTH")
    {
      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_ALL);
    }
    else
    {
      errors << priority(error) << "Value for parameter, output_pedigrees, "
             << "not recognized.  Ignoring ..." << endl;
    }
  }
}

/// Parses the 'SAMPLE_DETAIL' parameter in the mlod analysis block
///
void
Parser::parse_ind_sample_table_option(const LSFBase* param)
{
  string  value;
  //lint -e{534} Ignoring return type.
  parse_string(param, value);

  if(value.size())
  {
    string cap_value = toUpper(value);
    if(cap_value == "NONE")
    {
      my_parameters.set_ind_sample_table_option(AnalysisParameters::IS_NONE);
    }
    else if(cap_value == "REMOVED")
    {
      my_parameters.set_ind_sample_table_option(AnalysisParameters::IS_REMOVED);
    }
    else if(cap_value == "ALL")
    {
      my_parameters.set_ind_sample_table_option(AnalysisParameters::IS_ALL);
    }
    else
    {
      errors << priority(error) << "Value for parameter, output_pedigrees, "
             << "not recognized.  Ignoring ..." << endl;
    }
  }
}

/// Prints the header of a new analysis block
///
void
Parser::print_header() const
{
  my_messages << "Beginning new analysis block ...\n" << endl;
}

/// Prints the final status of the parse of the analysis block and
/// a final footer.
void
Parser::print_footer() const
{
  my_messages << "\nAnalysis parsing complete.  ";

  if(my_parameters.is_valid())
    my_messages << "Analysis valid." << endl;
  else
    my_messages << "Analysis is invalid.  Please check your file." << endl;

  my_messages << "------------------------------------------------------------------------------\n"
           << endl;
}

void
Parser::produce_bad_parameter_warning(const string& param_name)
{
  errors << priority(error) << "Parameter '" << param_name
         << "' not recognized.  Ignoring ..." << endl;
}

bool
Parser::is_trait_marker_new(size_t trait_marker_id) const
{
  const AnalysisParameters::TraitModelList& tr = my_parameters.get_trait_list();

  AnalysisParameters::TraitModelList::const_iterator  iter;

  iter = find_if(tr.begin(), tr.end(), boost::bind(is_same_index, _1, trait_marker_id));
  
  return iter == tr.end();
}

// - Check constraints between parameters.
//
void
Parser::validate()
{
  validate_ped_detail();
  validate_trait_count();
  validate_region();
}

// - scan_type must be 'interval' if output_pedigrees is 'all'.
//
void
Parser::validate_ped_detail()
{
  if( my_parameters.get_ped_output_detail_option() == AnalysisParameters::PD_NONE )
    return;

  if( my_parameters.get_scan_type() == AnalysisParameters::ST_MARKER )
  {
    if( my_parameters.get_ped_output_detail_option() == AnalysisParameters::PD_ALL )
    {
      errors << priority(error) << "output_pedigrees may not be 'both' when "
             << "scan_type is 'marker'.  Setting output_pedigrees to 'marker' ..." 
             << endl; 

      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_MARKERS);
    }
    else if( my_parameters.get_ped_output_detail_option() == AnalysisParameters::PD_INTERVALS )
    {
      errors << priority(error) << "output_pedigrees may not be 'interval' when "
             << "scan_type is 'marker'.  Setting output_pedigrees to 'marker' ..." 
             << endl; 

      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_MARKERS);
    }
  }
  else if( my_parameters.get_scan_type() == AnalysisParameters::ST_INTERVAL )
  {
    if( my_parameters.get_ped_output_detail_option() == AnalysisParameters::PD_ALL )
    {
      errors << priority(error) << "output_pedigrees may not be 'both' when "
             << "scan_type is 'interval'.  Setting output_pedigrees to 'interval' ..." 
             << endl; 

      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_INTERVALS);
    }
    else if( my_parameters.get_ped_output_detail_option() == AnalysisParameters::PD_MARKERS )
    {
      errors << priority(error) << "output_pedigrees may not be 'marker' when "
             << "scan_type is 'interval'.  Setting output_pedigrees to 'interval' ..." 
             << endl; 

      my_parameters.set_ped_output_detail_option(AnalysisParameters::PD_INTERVALS);
    }
  }
}

// - There MUST be at least one trait_marker specified.
//
void
Parser::validate_trait_count()
{
  if(my_parameters.get_trait_list().empty())
  {
    errors << priority(critical) << "No valid trait markers specified.  "
           << "Skipping analysis ..." << endl;
  }
}

// - There MUST be at least one region specified.
//
void
Parser::validate_region()
{
  if(!my_parameters.get_region().valid())
  {
    errors << priority(critical) << "No valid region specified.  "
           << "Skipping analysis ..." << endl;
  }
}

}
}
