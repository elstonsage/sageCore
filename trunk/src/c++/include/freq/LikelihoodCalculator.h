#ifndef FREQ_LHCALC
#define FREQ_LHCALC

#include "maxfunapi/maxfunapi.h"
#include "freq/Peeler.h"
#include "freq/Sample.h"
#include "freq/TransmissionProbCalc.h"
#include "freq/GenotypeProbCalc.h"
#include "freq/SingletonCalc.h"

namespace SAGE {
namespace FREQ {

class LikelihoodCalculator
{
  public:

  /// @name Misc.
  //@{
  
    ///
    /// Constructor
    ///
    /// NOTE: Assumes that the sample was built on a multipedigree with NO loops.
    ///
    /// \param sample The sample to use
    /// \param err The errorstream to use
    LikelihoodCalculator(const Sample & sample, MAXFUN::ParameterMgr & mgr, bool use_inbreeding, const cerrorstream & err = sage_cerr);

    LikelihoodCalculator(const LikelihoodCalculator & other);

    ///
    /// Destructor.
    ~LikelihoodCalculator();

    ///
    /// Calculates dependent parameters (in particular, the 'last' allele frequency).
    int calculateDependents(MAXFUN::ParameterMgr & mgr) const;

    ///
    /// Tests that genotype frequencies as given by the allele frequencies and
    /// inbreeding parameter are all > 0.0
    ///
    /// Returns 0 if all are ok, 1 if the first negative value is homozygous and 2 if 
    /// the first is heterozygous.
    int testInbredGenotypeFrequencies(MAXFUN::ParameterMgr& mgr) const;
    
    ///
    /// Calculates the likelihood given the constructor data and the
    /// current estimates in the ParameterMgr.
    double calculateLikelihood(MAXFUN::ParameterMgr & mgr) const;
    
    ///
    /// Indicates whether or not inbreeding coeff. is to be estimated.
    bool getEstimateInbreeding() const { return my_use_inbreeding; }

  //@}

private:

  LikelihoodCalculator& operator=(const LikelihoodCalculator &); // Disallowed

  /// Error stream
  cerrorstream my_errors;

  /// Phenoset data source
  const Sample & my_sample;
  
  typedef std::map<FPED::SubpedigreeConstPointer, boost::shared_ptr<Peeler> > PeelerMap;

  /// Peelers for each subpedigree:
  mutable PeelerMap my_peelers;
  
  /// Use inbreeding?
  bool my_use_inbreeding;
  
  /// Singleton calculator:
  mutable SingletonCalc my_singleton_calc;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
