#ifndef REG_UNCONNECTED_LIKELIHOOD
#include "segreg/RegUnconnectedLikelihood.h"
#endif

namespace SAGE {
namespace SEGREG {

inline
RegUnconnectedLikelihood::RegUnconnectedLikelihood
  (const RegPenetranceCalculator& rpc,
   const model& modp)
  : pen(rpc),
    mod(modp)
{}

}
}
