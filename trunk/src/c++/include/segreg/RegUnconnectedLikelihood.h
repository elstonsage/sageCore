#ifndef REG_UNCONNECTED_LIKELIHOOD
#define REG_UNCONNECTED_LIKELIHOOD

#include <mped/sp.h>
#include <mped/mp.h>
#include "peeling/peeler3.h"
#include "numerics/log_double.h"
#include "segreg/model.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/freq_sub_model.h"
#include "segreg/transmission_sub_model.h"
#include "segreg/RegPenetranceCalculator.h"
#include "segreg/polygenic_penetrance_calculator.h"
#include "segreg/peeling_caches.h"
#include "segreg/polygenic_transition_calculator.h"

namespace SAGE {
namespace SEGREG {

class RegUnconnectedLikelihood
{
  public:
    typedef RegPenetranceCalculator::penetrance_info   penetrance_info;

    RegUnconnectedLikelihood(const RegPenetranceCalculator& pen,
                             const model& modp);

    // -------------------------------------------------
    // unconnected likelihood for unconnected individual
    // -------------------------------------------------

    log_double unconnected_likelihood(int                      genotype, 
				      FPED::MemberConstPointer memit)
    {
        log_double p(1.0);

        double psi = mod.freq_sub_model.prob(genotype);

        p = psi * pen.get_penetrance(penetrance_info(*memit,genotype));

        return p;
    }

  protected:  

    const RegPenetranceCalculator& pen;
    const model &                  mod;
};

}
}

#include "segreg/RegUnconnectedLikelihood.ipp"

#endif

