#ifndef ANALYSIS_DATA_H
#define ANALYSIS_DATA_H

//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   5/24/02 created                     - gcw
//            5/28/02 modified                    - djb
//                                                                          
// Notes      Defines classes to hold analysis parameters and analysis data.
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <list>
#include <string>
#include <limits>
#include "fped/fped.h"
#include "mlod/analysis_parameters.h"
#include "mlod/PedigreeAnalysisSample.h"
#include "mlod/analysis_results.h"
#include "mlod/lod_table.h"

namespace SAGE
{
namespace MLOD
{

// Forward declaration
class AnalysisDataImpl;

/// \brief Provides public access to the data associated with a particular analysis
///
/// The AnalysisData objects provides the public interface to the data associated
/// with a particular analysis.  This includes the AnalysisParameters 
/// (the parameters), PedigreeAnalysisSample (the data set after filtering), and
/// AnalysisResults (results of performing the analysis).
///
/// \internal
///
/// Due to the large size of the data set and results, the analysis' data is
/// not actually stored in the AnalysisData object.  Instead, the AnalysisDataImpl
/// holds the data and the AnalysisData holds a shared pointer to it.  This
/// makes passing the data in and out of Analyzer objects much more efficient.
class AnalysisData
{
    friend class Analyzer;
  
  public:
  
    /// \name Object Management
    //@{
    /// Constructor
    ///
    AnalysisData();

    /// Copy Constructor
    ///
    AnalysisData(const AnalysisData& );
    
    /// Destructor
    ///
    ~AnalysisData();
    
    /// Copy Operator
    ///
    AnalysisData& operator=(const AnalysisData&);
    //@}
    
    /// Returns access to the AnalysisParameters of the analysis
    ///
    const AnalysisParameters&     get_parameters()      const;

    /// Returns access to the PedigreeAnalysisSample of the analysis
    ///
    const PedigreeAnalysisSample& get_pedigree_sample() const;

    /// Returns access to the AnalysisResults of the analysis
    ///
    const AnalysisResults&        get_results()         const;

  private:
  
    /// Shared pointer to the AnalysisDataImpl
    ///
    typedef boost::shared_ptr<const AnalysisDataImpl> AnalysisDataImplShPtr;
  
    /// Primary Constructor
    ///
    /// This Constructor is used by the Analyzer to create new AnalysisData objects
    /// which point to AnalysisDataImpl objects.  The AnalysisDataImpl objects
    /// themselves are also created by the Analyzer.
    ///
    /// \param data_store The storage object to link to.
    //lint -e{1704} Only Analysis objects can actually create these with a data store.
    AnalysisData(const AnalysisDataImplShPtr& data_store);
  
    /// Shared pointer to the current AnalysisDataImpl
    ///
    AnalysisDataImplShPtr my_data_storage;
};

/// \internal
///
/// \brief Stores all data associated with an analysis
///
/// The AnalysisDataImpl responsible for constructing and containing analysis
/// specific information.  Given parameters and data, it constructs the data
/// structures that will be directly analyzed (PedigreeDataSample and an
/// analysis specific genome containing the region and traits to be analyzed)
/// and provides access to this data to the Analyzer.
///
/// Since these structures are large and copying is impractical, AnalysisDataImpl
/// objects are hidden and created and passed around by shared pointer within
/// AnalysisData objects.  AnalysisData objects provide the public interface to
/// access this information to client code.  Construction and direct access
/// can then be limited to the Analysis classes safely despite the public
/// interface.
class AnalysisDataImpl
{
  public:
  
    AnalysisParameters       my_parameters;       ///< Storage of the analysis' parameters
                                                  ///<
    PedigreeAnalysisSample   my_pedigree_sample;  ///< Storage of the pedigree data
                                                  ///< valid for the analysis specified by
                                                  ///< my_parameters
    AnalysisResults          my_results;          ///< Storage of the results of
                                                  ///< performing the analysis
    RPED::genome_description my_genome;           ///< genome description for the
                                                  ///< analysis
  private:
  
    friend class Analyzer;
    friend class AnalysisData;

  /// \name Object Management
  //@{
    /// Constructor
    ///
    /// Given a set of parameters and the source multipedigree,
    /// the AnalysisDataImpl copies the former, and uses the latter
    /// to create the analysis' pedigree data set which is also stored.
    ///
    /// \param params     The AnalysisParameters
    /// \param source_ped The source RPED::MultiPedigree of the analysis
    /// \param errors     Stream which recieves errors (needed by my_genome)
    //lint -e{1704} Only Analysis and AnalysisData objects can actually create these
    AnalysisDataImpl(const AnalysisParameters&  params,
                     const RPED::MultiPedigree& source_rped,
                     cerrorstream&              errors);

    /// Default constructor disabled
    ///
    AnalysisDataImpl();

    /// Copy constructor disabled
    ///
    AnalysisDataImpl(const AnalysisDataImpl&);

    /// Copy operator disabled
    ///
    AnalysisDataImpl& operator=(const AnalysisDataImpl);
  //@}
  
  /// \name Analysis Specific Genome Creation
  //@{
    /// Creates the analysis specific genome.  This genome contains two
    /// regions, REGION and TRAITS.
    ///
    /// REGIONS is identical to the region stored
    /// in the (global) Data object, but differs in the scan distance as specified
    /// by my_parameters.get_scan_distance()
    ///
    /// TRAITS is a region containing all the trait_markers specified for the
    /// analysis.
    ///
    /// These regions are created using build_analysis_specific_region() and
    /// build_analysis_trait_region(), respectively.
    void setup_genome();
    
    /// Adds a region "REGION" to my_genome which is identical to the region
    /// stored in the global genome, but has a distance as specified by
    /// my_parameters.get_scan_distance()
    void build_analysis_specific_region();
    
    /// Adds a region "TRAITS" to my_genome which contains each of the
    /// traits specified in my_parameters.get_trait_list()
    void build_analysis_trait_region();
  //@}
  
};

}
}
#include "mlod/analysis_data.ipp"

#endif
