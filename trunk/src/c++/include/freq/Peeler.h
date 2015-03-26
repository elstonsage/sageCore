#ifndef FREQ_PEELER
#define FREQ_PEELER

#include <vector>
#include <list>
#include "error/errormanip.h"
#include "error/errorstream.h"
#include "containers/bitfield.h"
#include "numerics/log_double.h"
#include "peeling/peeler3.h"
#include "util/AutoTrace.h"
#include "rped/loop.h"   
#include "rped/rped.h"
#include "maxfunapi/maxfunapi.h"
#include "freq/Sample.h"
#include "freq/Cache.h"
#include "freq/AlleleRemapping.h"
#include "freq/TransmissionProbCalc.h"
#include "freq/GenotypeProbCalc.h"

namespace SAGE {
namespace FREQ {

typedef peeling::peeler<MLOCUS::unphased_genotype, log_double, peeling::individual_cache<MLOCUS::unphased_genotype, log_double> > PeelerBase;

/// \brief Calculates the freq likelihood
///
/// \par How the likelihood is calculated
///
/// The multipedigree likelihood, expressed on a log scale, is the sum of the
/// pedigree likelihoods. A pedigree likelihood, in turn, is the sum of its
/// subpedigree likelihoods and unconnecteds' likelihoods.
///
/// \par Unconnecteds
///
/// The likelihood for a singleton (unconnected individual) is the sum of the
/// probabilities of the possible genotypes for that individual's phenotype.
///
/// \par Subpedigrees
///
/// The likelihood for a subpedigree is calculated via the anterior/posterior
/// system. (If you're not familiar with this system, read up on Elston-Stewart
/// algorithm).
///
/// The anterior for a founder is the same as that for an unconnected.
///
/// The anterior for a non-founder is 
class Peeler : PeelerBase
{
  public:

  /// @name Constructor and likelihood
  //@{
  
    ///
    /// Constructor
    /// \param mgr The MAXFUN::ParameterMgr that stores the parameter estimates.
    /// \param mped The multipedigree source data
    /// \param phenosets The data structure describing the phenosets corresponding to each individual
    /// \param err The errorstream to use
    Peeler(
      const MAXFUN::ParameterMgr   & mgr,
      const Sample::SpedInfo       & sped_info,
      const Sample                 & sample,
            bool                     use_inbreeding,
            cerrorstream           & err       = sage_cerr);
            
    Peeler(const Peeler & other);
    
    virtual ~Peeler();

    ///
    /// Calculates the likelihood given the constructor data and the
    /// current estimates in the ParameterMgr.
    log_double calculateLikelihood();
    
  //@}

  /// @name Debug
  //@{

    void dump();
  
  //@}

private:

  Peeler& operator= (const Peeler & other); // Disallowed
  
  /// @name Required virtual interface
  //@{

    const result_type & internal_anterior              (const member_type& ind,                          const data_type & g, result_type & ia);
    const result_type & internal_anterior_terminal     (const member_type& ind,                          const data_type & g, result_type & iatg);
    const result_type & internal_posterior             (const member_type& ind,                          const data_type & g, result_type & ip);
    const result_type & internal_posterior_with_mate   (const member_type& ind, const member_type& mate, const data_type & g, result_type & ipwm);
    const result_type & internal_posterior_except_mate (const member_type& ind, const member_type& mate, const data_type & g, result_type & ipem);
    const result_type & internal_posterior_terminal    (const member_type& ind,                          const data_type & g, result_type & ipt);

  //@}

  /// @name Helper functions
  //@{

    log_double calculateAnteriorParent(const member_type & ind, const member_type & spouse_to_exclude, const GenotypeInfo & g);

    log_double calculateAnteriorChild(const member_type & ind, const data_type & g_m, const data_type & g_f);
  
  //@}

  /// Error stream
  cerrorstream my_errors;

  /// Phenoset data source
  const Sample & my_sample;

  /// The ParameterMgr with parameter estimates
  const MAXFUN::ParameterMgr & my_mgr;

  /// The calculator for getting transmission probs
  TransmissionProbCalc my_transm_prob_calc;
  
  /// The calculator for getting genotype probabilities
  GenotypeProbCalc my_geno_prob_calc;

  /// The allele remapping
  AlleleRemapping my_allele_remapping;

  /// The storage system for remapped allele freqs
  RemappedAlleleFreqs my_remapped_allele_freqs;
  
  /// The likelihood for this peeler
  log_double my_lh;
  
  /// Number of times this peeler's lh has been evaluated
  size_t my_lh_count;
  
  /// Bit field indicating which alleles are present (where each bit corresponds to the i'th un-remapped allele id)
  bit_field my_relevant_alleles;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
