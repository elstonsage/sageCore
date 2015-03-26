#ifndef LOD_SCORE_ANALYZER_H
#define LOD_SCORE_ANALYZER_H

#include "boost/bind.hpp"
#include "app/output_streams.h"
#include "rped/genome_description.h"
#include "fped/fped.h"
#include "gelim/geno_eliminate.h"
#include "gelim/pedigree_region.h"
#include "lvec/lvector.h"
#include "lvec/mpoint_like.h"
#include "mlod/analysis_results.h"
#include "mlod/analysis_data.h"

namespace SAGE {
namespace MLOD {

/// \brief Computes Lod scores and information content for subpedigrees along a marker region
///
/// The LodScoreAnalyzed calculates the Lod Scores for a subpedigree for each
/// point in a region based upon the AnalysisParameters given to it.
/// 
/// For efficiency purposes, it pre-allocates and stores all the memory needed
/// to do this calculation for the largest subpedigree it finds in the
/// PedigreeAnalysisSample.  It should \b never be copied or moved.
///
class LodScoreAnalyzer
{
  public:
  
    /// Convenience typedef to give access to the SpedIterator
    ///
    typedef PedigreeAnalysisSample::SpedIterator SpedIterator;
    
  /// \name Object Management
  //@{
    
    /// Public Constructor
    ///
    /// \param data The Analysis Data to be analyzed
    /// \param out  Where output should be sent 
    LodScoreAnalyzer(const AnalysisDataImpl& data,
                     APP::Output_Streams&    out);

    /// Destructor
    ///
    ~LodScoreAnalyzer();
  //@}

    /// Returns \c true if construction was successful, \c false otherwise
    /// -- Change to exception?
    bool is_valid() const;

    /// Sets the target where results are to be placed.  Note that
    /// this is a member of the AnalysisDataImpl, in general, but
    /// this gives explicit, non-const access to the results object.
    ///
    /// \param results The target for results
    void set_result_target(AnalysisResults& results);
  
    /// Run the analysis on a specific subpedigree
    ///
    /// \param sp The subpedigree to run the analysis on
    bool analyze_subpedigree(SpedIterator sp);
  
  private:
  
  /// \name (Disabled) Object Management
  //@{
    /// Default Constructor disabled
    ///
    LodScoreAnalyzer();
    
    /// Copy Constructor disabled
    ///
    LodScoreAnalyzer(const LodScoreAnalyzer&);
    
    /// Copy Operator disabled
    ///
    LodScoreAnalyzer operator=(const LodScoreAnalyzer&);
  //@}
  
    typedef Likelihood_Vector LikelihoodVector;
  
//    void build_trait_genome();
    void request_marker_resources();
    void build_likelihood_data_structures();
    void build_temporary_likelihood_vectors();
  
    void initialize_meiosis_map(const FPED::Subpedigree& s);
    void initialize_marker_data(const FPED::Subpedigree& s);
    void initialize_trait_data (const FPED::Subpedigree& s);
  
    void compute_lod_scores(SpedLodTable& sptable,
                            const Likelihood_Vector& mkr_vect,
                            size_t sptable_index);
  
    /// Computes the multipoint likelihood vector at point point_idx in the
    /// interval between marker_idx and marker_idx+1
    ///
    /// The resulting likelihood vector is stored in temp1.
    ///
    /// \param marker_idx The index of the marker to the left
    /// \param point_idx  The index of the point being calculated
    /// \returns Nothing.  The results are stored in the temp1 vector.
    void compute_interval_multipoint_lvec(size_t marker_idx, size_t point_idx);
    
    size_t num_loci()   const;
    size_t num_points() const;
    size_t num_traits() const;
    size_t num_bits()   const;
    
    bool using_intervals() const;
  
    const AnalysisDataImpl&       my_data;
    const AnalysisParameters&     my_parameters;
    const PedigreeAnalysisSample& my_peds;
    
    RPED::genome_description::region_type my_traits;
    
    AnalysisResults* my_result_target;
    
    APP::Output_Streams&  my_output;
    
    bool my_valid;
  
    // subpedigree specific data structures.
    
    meiosis_map my_meiosis_map;
  
    // likelihood vector data structures
    
    LikelihoodVector temp1, temp2;
    
    mpoint_likelihood_data         my_marker_data;
    mpoint_likelihood_data         my_trait_data;
};

} // end namespace MLOD
} // end namespace SAGE

#include "mlod/lod_score_analyzer.ipp"

#endif
