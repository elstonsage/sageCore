#ifndef PARENT_OF_ORIGIN_H
#define PARENT_OF_ORIGIN_H

#include "rped/rped.h"
#include "func/FunctionParser.h"

namespace SAGE {
namespace FUNC {

/// \brief Class that calculates the (maternal) parent of origin statistic.
/// 
/// From the paper: ... we assign a covariate value of 1 to an individual whose mother
/// is a non-founder and whose father is a founder (i.e., married into the family);
/// and we assign a covariate value of -1 to an individual whose father is a
/// non-founder and whose mother is a founder. We assign a code value of 0 if the
/// individual is a founder or if both parents are non-founders. The value 0, which
/// is halfway between -1 and 1, reflects a lack of knowledge about the source of a
/// deleterious allele that is transmitted in the pedigree

class ParentOfOriginCalculator
{
public:

  /// @name Creating the POO status variable
  //@{

    /// Creates a POO trait within the indicated RefMultiPedigree.
    ///
    /// \param mp The RefMultiPedigree instance in which to create the POO status variable
    /// \param data Parser data indicating the trait name and such for creation of the trait
    static void createPooStatusTrait(RPED::RefMultiPedigree & mp, const FunctionParser::TraitData & data);
    
  //@}

private:

  static double getPooStatus(const RPED::RefMember & mem);

};    
    
} // End namespace FUNC
} // End namespace SAGE

#endif
