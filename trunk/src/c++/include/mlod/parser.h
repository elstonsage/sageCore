#ifndef MLOD_NEW_PARSER_H
#define MLOD_NEW_PARSER_H

//============================================================================
// File:      parser.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/28/02 created                     - djb
//                                                                          
// Notes:    Defines class, parser, for parsing mlod analysis block of
//            a SAGE parameter file.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <list>
#include <algorithm>
#include <string>
#include <iostream>
#include <cassert>
#include "rped/genome_description.h"
#include "rped/rped.h"
#include "app/aparser.h"
#include "LSF/LSF.h"
#include "LSF/Attr.h"
#include "mlod/definitions.h"
#include "mlod/analysis_parameters.h"

namespace SAGE
{
namespace MLOD
{

class Data;
  
///  Extracts information from mlod analysis block of SAGE parameter
///  file and creates the AnalysisParameters object
class Parser : public APP::BasicParser
{
  public:
  
  /// \name Object Management
  ///
  //@{
    /// Constructor.
    ///
    /// This is based upon the constructor for the APP::BasicParser
    Parser(const Data& data, ostream& messages = std::cout,
           cerrorstream& errors = sage_cerr);
  //@}
           
  /// \name Standard Interface
  //@{
    /// This function is provided to simplify use, and is not in the 
    /// standard interface.
    void  parse(const LSFBase* params);
  
    /// Declared in SAGE::BasicParser as virtual base, so must be overridden here.
    ///
    virtual void  parse_symbols(const SymbolTable* syms);
    /// Declared in SAGE::BasicParser as virtual base, so must be overridden here.
    ///
    virtual void  parse_parameter(const LSFBase* param);
    /// Declared in SAGE::BasicParser as virtual base, so must be overridden here.
    ///
    virtual void  parse_test_parameter_section(const LSFBase* params);
    /// Declared in SAGE::BasicParser as virtual base, so must be overridden here.
    ///
    virtual void  parse_test_parameter(const LSFBase* param);
  //@}
  
  /// \name Access to parsed data
  //@{
    /// Returns the most recently parsed analysis
    ///
    const AnalysisParameters&  parameters() const;  
    
    /// Returns the current analysis id.  This is a counter which starts at
    /// 1 and counts for each parsed analysis.
    size_t  analysis_id() const;
  //@}
    
  private:
  
    /// \internal
    /// Default Contructor (disabled) 
    ///
    //lint -e{1704} Default construction not allowed.
    Parser();

    /// Resets to a default analysis state.
    ///
    void  reset_parameters();

    void  parse_out(const LSFBase* params ,const string& params_name);
    void  parse_trait_marker            (const LSFBase* param);
    void  parse_region                  (const LSFBase* param);
    void  parse_max_ped_size            (const LSFBase* param);
    void  parse_scan_type               (const LSFBase* param);
    void  parse_distance                (const LSFBase* param);
    void  parse_output_pedigrees_option (const LSFBase* param);
    void  parse_ind_sample_table_option (const LSFBase* param);

    void  check_distance_value(double d_value);

    void  print_header() const;
    void  print_footer() const;

    /// Sends a message to the error stream when an unrecognized parameter is
    /// encountered while parsing.
    void  produce_bad_parameter_warning(const string& param_name);

    /// Returns \c true if the trait marker is not in the current
    /// AnalysisParameter's trait list, \c false otherwise.
    bool  is_trait_marker_new (size_t trait_marker_id) const;

    /// Check the parsed analysis for a few final errors/warnings for
    /// the user.  Make any final changes that can be made to fix the
    /// analysis when invalid states are encountered, if possible.
    void  validate();
    
    /// Checks the pedigree detail is not inconsistent with the scan type.
    /// If pedigree detail is too large, scan type is adjusted and a warning
    /// produced.
    void  validate_ped_detail();
    
    /// Checks to make sure there are valid traits to be analyzed and produce
    /// an invalid analysis message if there aren't.
    void  validate_trait_count();

    /// Checks to make sure the region has been specified, and produce
    /// an invalid analysis message if it hasn't.
    void  validate_region();
  
    // Data members.
    const Data&            my_data;          ///< MLOD data for checking traits, etc.
    ostream&               my_messages;      ///< Where to send parsing status messages.
                                             ///< (note that errors comes from APP::BasicParser

    long                   my_analysis_id;   ///< The counter for assigning analysis ids
    
    AnalysisParameters     my_parameters;    ///< The currently worked on parameters.
};

#include "mlod/parser.ipp"

}
}

#endif

