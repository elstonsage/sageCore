#ifndef MLOD_ANALYSIS_PARAMETERS_H
#define MLOD_ANALYSIS_PARAMETERS_H

//============================================================================
// File:      analysis_parameters.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   5/24/02  created                        - gcw
//            5/28/02  modified                       - djb
//            11/15/04 Significantly modified/updated - gcw
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <list>
#include <string>
#include <iostream>
#include <iomanip>

#include "error/internal_error.h"
#include "output/Output.h"
#include "mlocus/imodel.h"
#include "rped/genome_description.h"
#include "mlod/definitions.h"

namespace SAGE
{
namespace MLOD
{

/// \brief Storage of analysis parameters
///
/// This class stores all the parameters which define an MLOD analysis.
/// This includes analysis information like traits and regions to use, as well
/// as meta information like title, output options, and so on.
class AnalysisParameters
{
  public:
  
    /// The TraitModelPair has little to do with traits, and everything to
    /// do with inheritance models.  It stores a pair of elements, the first being
    /// the index of the model in the global inheritance model map, and the
    /// second being the global penetrance model stored there.
    ///
    /// In general, this will represent a 'trait model', which is a model we wish
    /// to place within the marker map using MLOD.  Sometimes, it might represent a
    /// marker.
    typedef std::pair<size_t, MLOCUS::penetrance_model> TraitModelPair;

    /// A TraitModelList represents an unordered group of inheritance models
    /// that must be processed for a given analysis.  Typically this is a set
    /// of inheritance models which represent a set of models for some trait or traits
    /// of interest.
    typedef std::list<TraitModelPair>                   TraitModelList;
  
    /// The RegionType is just a simplifying typedef
    ///
    typedef RPED::genome_description::region_type RegionType;
  
    /// \name Object Management
    //@{
    
    /// Constructor (default)
    ///
    /// \param analysis_id The id of the analysis.  USed to specify default options
    explicit AnalysisParameters(size_t analysis_id = 1);
    
    /// Copy Constructor
    ///
    AnalysisParameters(const AnalysisParameters&);
    
    /// Copy Operator
    ///
    AnalysisParameters& operator=(const AnalysisParameters&);
    
    /// Destructor
    ///
    ~AnalysisParameters();
    //@}
    
    /// Valididy check
    ///
    /// \returns \c true if the object contains a valid analysis specification,
    ///          \c false, otherwise.
    bool  is_valid() const;

    /// \name Output Management
    ///
    /// Output management comprises the title (as use within output files) and
    /// options which control the naming of output files.  Output files are named
    /// using a base name to which are appended various extensions indicating
    /// file type.
    //@{
    /// Returns the title
    ///
    const string& get_title() const;
    
    /// Sets the title
    ///
    void set_title(const string& s);
      
    /// Returns the base name used for output files.  Extensions (.sum, .det, etc.)
    /// should be appended to this name to create output filenames 
    /// for this analysis.
    const string& get_base_output_filename() const;

    /// Sets the base name used for output files.  Extensions (.sum, .det
    /// as needed) should be appended to this name to create 
    /// output filenames for this analysis.
    ///
    /// \param s The string to be used as a base filename
    void set_base_output_filename(const string& s);
    
    /// Ennumeration describing the detail of the per pedigree output desired
    ///
    enum PedOutputDetailEnum
    {
      PD_NONE,      ///< No per-pedigree output desired (DEFAULT)
      PD_MARKERS,   ///< Output lod scores only at the markers for each pedigree
      PD_INTERVALS, ///< Output lod scores only at the intervals for each pedigree
      PD_ALL        ///< Output all lod scores (marker and interval) for each pedigree
                    ///< (NOTE: interval information only produced when intervals are available)
    };
    /// Returns a PedOutputDetailEnum option indicating the amount of detail 
    /// desired in the per-pedigree output.
    PedOutputDetailEnum get_ped_output_detail_option() const;
    
    /// Returns a string equivalent to the PedOutputDetailEnum
    ///
    string get_ped_output_detail_string() const;

    /// Sets the pedigree output detail
    ///
    /// \param p Level of detail desired
    void set_ped_output_detail_option(PedOutputDetailEnum p);
    //@}

    /// Ennumeration describing the detail of the individual status
    /// summary table
    enum IndSampleTableEnum
    {
      IS_NONE,    ///< No table desired
      IS_REMOVED, ///< Output only removed individuals (DEFAULT)
      IS_ALL      ///< Output all individuals
    };
    /// Returns a IndSampleTableEnum option indicating the amount of detail 
    /// desired in the individual sample table output.
    IndSampleTableEnum get_ind_sample_table_option() const;
    
    /// Returns a string equivalent to the IndSamleTableEnum
    ///
    string get_ind_sample_table_string() const;

    /// Sets the pedigree output detail
    ///
    /// \param p Level of detail desired
    void set_ind_sample_table_option(IndSampleTableEnum i);
    //@}

    /// \name Analysis Control
    //@{
  
    /// Lists options specifying the level of scanning to perform
    ///
    enum ScanTypeEnum
    {
      ST_MARKER,   ///< Calculate lod scores only at markers 
      ST_INTERVAL, ///< Calculate lod scores at intervals between markers
      ST_BOTH      ///< Calculate lod scores at markers and at intervals between markers
    };
    
    /// Returns ScanTypeEnum indicating the type of scan to perform
    ///
    ScanTypeEnum get_scan_type() const;
    
    /// Returns string-equivalent of the ScanType to perform
    ///
    string       get_scan_type_string() const;
    
    /// Sets the scan type
    ///
    /// \param s Scan type to set to
    void set_scan_type(ScanTypeEnum s);

    /// Returns the distance in centimorgans to scan (if interval scanning is specified)
    ///
    double get_scan_distance() const;
    
    /// Sets the scan distance.  Only used if interval scanning is specified
    ///
    /// \param d the distance to set.
    void set_scan_distance(double d);
    
    /// Returns the maximum 2n-f a constituent pedigree can have.
    ///
    size_t get_max_ped_size() const;
    
    /// Sets the max 2n-f size for the analysis.
    ///
    /// \param sz the size to set.
    void set_max_ped_size(size_t sz);

    /// Gets access to the region to analyze
    ///
    const RegionType& get_region() const;
    
    /// Sets the region to analyze
    ///
    /// \param reg The region to analyze
    void set_region(const RegionType& reg);

    /// Gets access to the traits to analyze
    ///
    TraitModelList& get_trait_list();
    
    /// Gets access to the traits to analyze
    ///
    const TraitModelList& get_trait_list() const;
    
    //@}
    
    /// Print the analysis to a stream in detail.  Primarily for user feedback
    ///
    /// \param out     The output stream to print to.
    /// \param markers The markers for getting trait locus names
    void print_analysis_table(std::ostream& out) const;
                                    
  private:

    string my_analysis_title;          ///< Title for screen use
                                       ///<
    string my_base_filename;           ///< String to which extentions are appended to create filenames
                                       ///<
  
    ScanTypeEnum my_scan_type;         ///< Scan type option indicating the level of
                                       ///< scanning.

    PedOutputDetailEnum my_ped_option; ///< Option indicating the output of per-pedigree output
                                       ///< desired, defaulting to PD_NONE.  Cannot
                                       ///< be PD_ALL if my_scan_type is ST_MARKERS
                                       
    IndSampleTableEnum  my_ind_option; ///< Option indicating the detail of the
                                       ///< individual status table output
                                
    double my_distance;                ///< Distance for interval use (if desired)
                                       ///< Default is 2 centimorgans.
                                       
    size_t my_max_ped_size;            ///< Max pedigree size
                                       ///<
                                       
    RegionType   my_region;            ///< region from the genome_description.
                                       ///< No default behavior, ie user must specify it.
                            
    TraitModelList my_trait_loci;      ///< indices and models in the 'inheritance models' in the RefMultiPedigree
                                       ///< No default behavior, ie user must specify at least one.
};

/// \internal
///
/// This function is provided as a helper function when we want to get the size_t
/// component of the AnalysisParameters::TraitModelPair.  We often need to do
/// this in STL functors where we can't easily access member info.
size_t get_trait_model_pair_first_element(const AnalysisParameters::TraitModelPair& p);

/// \internal
///
/// This function is provided as a helper function when we want to locate a size_t
/// indexed component of the AnalysisParameters::TraitModelPair in the
/// AnalysisParameters::TraitModelList.  We often need to do
/// this in STL functors where we can't easily access member info.
bool is_same_index(const AnalysisParameters::TraitModelPair& p, size_t idx);

}
}

#include "mlod/analysis_parameters.ipp"

#endif
