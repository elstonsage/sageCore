#ifndef AO_EXTRA_OUTPUT_H
#define AO_EXTRA_OUTPUT_H

#include "ageon/MemberCovariateCalculator.h"

namespace SAGE {
namespace AO   {

/// Represents the conditional probabilities and survival residuals for a single model.
struct ModelTraits
{
  std::vector<double>      cond_probs;
  std::vector<double>      surv_resid;
};

/// Vector where the index number is the model type and the value is a ModelTraits.
typedef std::vector<ModelTraits> ModelTraitsVector;

class ExtraOutput
{
  public:
  
    static void generateFile(
      const SAMPLING::PartitionedMemberDataSample & sample,
      const ModelTraitsVector                     & model_traits, 
      const std::string                           & analysis_name,
            std::ostream                          & info);

    static void populateModelTraits(const MemberCovariateCalculator & mcc, ModelTraits & model_traits);
};

} // End namespace AO
} // End namespace SAGE


#endif
