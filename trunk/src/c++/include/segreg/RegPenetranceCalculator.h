#ifndef REG_PENETRANCE_CALCULATOR_H
#define REG_PENETRANCE_CALCULATOR_H

#include <math.h>
#include "segreg/member_calculator.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/model.h"
#include "segreg/ContinuousPenetranceFunction.h"
#include "segreg/RegPenetranceCommonAdjustments.h"
#include "segreg/RegPenetranceContext.h"
#include "mped/mp_utilities.h"

namespace SAGE {
namespace SEGREG {

/// The RegPenetranceCalculator calculates the penetrance probabilities
/// under the regressive models for individuals (See eqs. 48, 48a, 49a, 63,
/// 64 and 66.)
///
/// The RegPenetranceCalculator calculated the penetrance probabilities for
/// individuals.  An optional RegPenetranceContext can be passed to the
/// penetrance function to allow for the computation of penetrances within a
/// particular context.
///
/// Contexts are controlled by the get_context() method.  This creates a new
/// context object of type RegPenetranceContext.  The context object is
/// manipulated by adding individuals specific to the context desired.  Once
/// complete, the context can be given as a function parameter to the
/// get_penetrance() function.
///
/// Keeping up to date is important.  The update() function updates internal
/// states based upon the current state of the model (which was supplied at
/// construction time).
class RegPenetranceCalculator
{
  public:
    typedef continuous_member_calculator::member_class member_class;

    typedef FPED::Member             member_type;
    typedef FPED::MemberConstPointer member_pointer;
    typedef FPED::Family             family_type;

    //==============================================================
    // Constructor/destructor/operators:
    //==============================================================

    /// The only supplied constructor.
    ///
    /// \param mped     The pedigree data that will be operated on.  Is used
    ///                 primarily by internal members.
    /// \param modp     The model used to perform calculations.
    /// \param use_asc  If the model includes ascertainment, it is indicated here.
    RegPenetranceCalculator(const FPED::Multipedigree & mped,
                            const model               & modp,
                            bool                        use_asc = false);

    //==============================================================
    // Structs:
    //==============================================================

    struct penetrance_info 
    { 
      penetrance_info(const member_type& ind,
                      genotype_index     genotype_param     = index_AA)
      { member = &ind; genotype = genotype_param; }
      penetrance_info()
      { member = NULL; genotype = index_AA; }

      penetrance_info(const penetrance_info& other)
      {
        member       = other.member;
        genotype     = other.genotype;
      }

      member_pointer member; 
      genotype_index genotype; 
    };

    int update();

    /// Returns a new, initially empty, context.  See RegPenetranceContext.
    PenetranceContext get_context() const;

    //================================================================
    // Penetrance functions:
    //================================================================

    /// Given individual, genotype, and context, calculate penetrance.
    ///
    /// \param i       penetrance information for individual we're calculating.
    /// \param context context under which the penetrance is calculated.
    double get_penetrance (penetrance_info i, const PenetranceContext& context)    const;

    /// Given individual and genotype, calculate penetrance.  Individual is
    /// independent of any specific context
    double get_penetrance (penetrance_info i) const;

    //================================================================
    // Private data members:
    //================================================================

  private:

    /// Calculate the number of prior sibs.
    ///
    /// \todo Come up with a better way to determine prior sib count. 
    ///       Currently this is done by counting them, which is inefficient.
    size_t get_prior_sib_count(const member_type& i) const;

    /// The class of the penetrance model.  Calculates only for class A and
    /// D, which differ only by the use of the siblings in family contexts
    model_class                           my_model_class;
    
    /// The residual sub model we use to get alphas, deltas, and
    /// correlations.
    const residual_correlation_sub_model& my_residuals;
    
    /// The member calculator calculates many of the values used by the
    /// penetrance calculator.  Since no other entity requires its services,
    /// it is hidden internally to the penetrance calculator.
    continuous_member_calculator          my_member_calc;

    /// Stores adjustments that are common to many contexts and can
    /// efficiently be precalculated.  A hidden member of the penetrance
    /// model.
    RegPenetranceCommonAdjustments              my_corr_adjs;
};

}
}

#include "segreg/RegPenetranceCalculator.ipp"

#endif
