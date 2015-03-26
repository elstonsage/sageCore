#ifndef ANALYSIS_RESULTS_H
#define ANALYSIS_RESULTS_H

//============================================================================
// File:      analysis_results.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <list>
#include <string>
#include <limits>
#include "fped/fped.h"
#include "mlod/PedigreeAnalysisSample.h"
#include "mlod/lod_table.h"

namespace SAGE
{
namespace MLOD
{

/// \brief Stores pedigree specific and summary lod tables for an analysis
///
/// The AnalysisResults object stores the results of an analysis.  It keeps
/// reference to the pedigree sample that generated the data, and stores tables
/// for each subpedigree in the analysis.  A summary table is also stored.
class AnalysisResults
{
  public:

    friend class Analysis;
    friend class AnalysisData;
    
    typedef PedigreeAnalysisSample::SpedIterator SpedIterator;
    
    /// Constructor
    ///
    /// \param ped_data The PedigreeAnalysisSample analyzed
    AnalysisResults(const PedigreeAnalysisSample& ped_data);

    /// Stores the table given in the entry for the subpedigree indicated
    /// and adds the table to the summary results.
    ///
    /// Behaviors when adding tables multiple times, or of different size
    /// (number of traits, number of points in region) is undefined.  Don't
    /// do it!
    ///
    /// \param i     The subpedigree whose table we're adding
    /// \param table The table to be stored
    void set_sped_lod_table(SpedIterator i, const SpedLodTable& table);
    
    /// Returns a reference to the table stored for the specific subpedigree
    ///
    /// \param i The subpedigree whose table we would like to access
    const SpedLodTable& get_sped_lod_table(SpedIterator i) const;

    /// Returns the Summary of all added subpedigree tables
    ///
    const SummaryLodTable& get_summary_lod_table() const;
    
  private:
  
    /// Constructor (disabled)
    ///
    //lint -e{1704} Default constructor disabled
    AnalysisResults();
  
    /// Reference to the data set analyzed
    ///
    const PedigreeAnalysisSample& my_pedigree_sample;

    /// Storage for the subpedigree lod tables for each subpedigree in
    /// the data set.
    vector<SpedLodTable> my_sped_lod_scores;

    /// The summary table summarizing all the subpedigrees in the data set
    ///
    SummaryLodTable      my_summary_lod_scores;
};

}
}

#include "mlod/analysis_results.ipp"

#endif

