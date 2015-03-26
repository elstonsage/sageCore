#ifndef FREQ_SINGLETON_CALC
#define FREQ_SINGLETON_CALC

#include <vector>
#include <list>
#include "error/errormanip.h"
#include "error/errorstream.h"
#include "numerics/log_double.h"
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

class SingletonCalc
{
  public:

    SingletonCalc(
      const MAXFUN::ParameterMgr   & mgr,
      const Sample                 & sample,
            bool                     use_inbreeding);
            
    SingletonCalc(const SingletonCalc & other);
    
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

  /// Phenoset data source
  const Sample & my_sample;

  /// The ParameterMgr with parameter estimates
  const MAXFUN::ParameterMgr & my_mgr;

  /// The calculator for getting genotype probabilities
  GenotypeProbCalc my_geno_prob_calc;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
