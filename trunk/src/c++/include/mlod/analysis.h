#ifndef MLOD_ANALYSIS_H
#define MLOD_ANALYSIS_H

//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include <list>
#include <string>
#include <limits>
#include "error/internal_error.h"
#include "app/output_streams.h"
#include "rped/rped.h"
#include "rped/genome_description.h"
#include "mlod/analysis_parameters.h"
#include "mlod/analysis_data.h"
#include "mlod/lod_score_analyzer.h"
#include "mlod/lod_table_print.h"

namespace SAGE
{
namespace MLOD
{

/// \brief Object which performs an MLOD analysis given parameters
///
/// The Analyzer performs the actual MLOD analysis which has been specified given
/// a particular set of parameters upon the basic data as provided in the constructor.
///
/// Currently, the results of this analysis do not need to be retained, as the Analyzer
/// produces output as the analysis is run.  However, future analyses may want or
/// require this output, so it's provided for future development.
class Analyzer
{
  public:

    //lint -e{1712} No Default Constructor

    /// Constructor
    ///
    /// \param rp      The data set
    /// \param genome  The genome
    /// \param out     The APP standard output streams
    Analyzer(const RPED::RefMultiPedigree& rp, const RPED::genome_description& genome, APP::Output_Streams& out);

    /// Runs an analysis given the parameters and produces output based on them.
    /// Because output is generated as the analysis is performed, no results need
    /// be returned at this time, but, just in case, the AnalysisData is returned.
    /// This could be used for future analyses (homogeneity?) that may not use
    /// this object.
    ///
    /// \param params The parameters of the analysis to be run
    /// \returns      The AnalysisData which includes all relevant analysis information.
    AnalysisData run_analysis(const AnalysisParameters& params) const;

  protected:
  
    typedef boost::shared_ptr<AnalysisDataImpl> DataShPtr;

    /// \name Analysis Initialization
    //@{
    void open_output_files      (const AnalysisParameters& params) const;
    void print_analysis_headers (const AnalysisParameters& params) const;
    //@}

    /// \name Analysis Processing
    //@{
    DataShPtr generate_data_impl(const AnalysisParameters& params) const;
    DataShPtr process_analysis(const AnalysisParameters& params) const;
    //@}

    /// \name Analysis Finalization
    //@{
    void close_output_files     ()                                 const;
    //@}

    const RPED::RefMultiPedigree&   my_peds;    ///< The unfiltered data set
    const RPED::genome_description& my_genome;  ///< The genome
    APP::Output_Streams&            my_out;     ///< The output streams
    
    mutable boost::shared_ptr<ofstream> my_sum;
    mutable boost::shared_ptr<ofstream> my_det;
};

}
}

#endif
