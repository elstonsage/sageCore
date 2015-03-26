#ifndef FREQ_ESTIMATOR
#define FREQ_ESTIMATOR


#include <string>
#include <set>
#include <map>
#include "globals/SAGEConstants.h"
#include "rped/rped.h"
#include "mped/mp_utilities.h"
#include "error/errormanip.h"
#include "error/errorstream.h"
#include "maxfunapi/maxfunapi.h"
#include "freq/Sample.h"
#include "freq/Peeler.h" 
#include "freq/Results.h"
#include "freq/Configuration.h"
#include "freq/LikelihoodCalculator.h"

namespace SAGE {
namespace FREQ {

/// \brief Carries out a single analysis.
///
/// Given a data source (a multipedigree) and a Configuration, the Estimator will 
/// perform two steps for every marker id in the configuration:
///
/// (1) If the marker is codominant, calculates the initial allele frequencies.
/// (2) Unless the skip_mle is true, calculates maximium likelihood estimates for the allele frequencies.
///
class Estimator
{
public:

    ///
    /// Executes the given analysis on the given data.
    ///
    /// \param config The configuration object describing the analysis to run
    /// \param mp The multipedigree storing the source data
    /// \param errors The errorstream to which messages should be directed
    static Results runAnalysis(const Configuration & config, const FPED::FilteredMultipedigree & mp, cerrorstream & errors = SAGE::sage_cerr);

private:

    ///
    /// Calculates the allele frequencies in the given sample.
    ///
    /// NOTE: Assumes that the marker in question is codominant (one genotype per phenotype!).
    static AlleleFrequencies calculateAlleleFrequencies(const Configuration & config, const Sample & sample);

    ///
    /// Maximizes the likelihood.
    /// NOTE: Assumes my_marker and my_gmodel have been set, and my_results[my_marker] exists.
    static void estimateLikelihood(const Configuration & config, const Sample & sample, MarkerAnalysis & analysis);

    static void setupComponents(
      const Configuration        & config, 
      const Sample               & sample,
      const AlleleFrequencies    & init_freqs,
            MAXFUN::ParameterMgr & mgr, 
            MAXFUN::Function     & func, 
            LikelihoodCalculator & lhcalc);

    static void runMaximize(MAXFUN::Function & func, const MAXFUN::SequenceCfg & seq, const MAXFUN::DebugCfg & dbg, MAXFUN::Results & results);

};

} // End namespace FREQ
} // End namespace SAGE

#endif
