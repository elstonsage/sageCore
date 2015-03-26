#include "func/ParentOfOrigin.h"

namespace SAGE {
namespace FUNC {

/// Specifically, we assign a covariate value of 1 to an individual whose mother
/// is a non-founder and whose father is a founder (i.e., married into the family);
/// and we assign a covariate value of -1 to an individual whose father is a
/// non-founder and whose mother is a founder. We assign a code value of 0 if the
/// individual is a founder or if both parents are non-founders. The value 0, which
/// is halfway between -1 and 1, reflects a lack of knowledge about the source of a
/// deleterious allele that is transmitted in the pedigree
double
ParentOfOriginCalculator::getPooStatus(const RPED::RefMember & mem)
{
  // Check basic conditions - founder, parents unsexed.
  if(mem.is_founder() == true)
    return 0;
  
  if(!mem.get_mother() || !mem.get_father())
    return 0;
  
  bool mother_is_founder = mem.get_mother()->is_founder();
  bool father_is_founder = mem.get_father()->is_founder();
  
  if(mother_is_founder == father_is_founder)
    return 0;
  
  return (mother_is_founder) ? -1 : 1;
}

void
ParentOfOriginCalculator::createPooStatusTrait(RPED::RefMultiPedigree & mp, const FunctionParser::TraitData& pd)
{
  // 1. Create the trait entry in both the RefMPedInfo:

  RPED::RefMPedInfo & mp_info   = mp.info();
  size_t              trait_num = mp_info.add_continuous_trait(pd.trait_name, pd.usage);

  mp_info.trait_info(trait_num).set_string_missing_code("");
  mp_info.trait_info(trait_num).set_numeric_missing_code(-999);

  // 2. Populate the RefMultiPedigree with values:

  for(RPED::RefMultiPedigree::pedigree_iterator ped = mp.pedigree_begin(); ped != mp.pedigree_end(); ++ped)
  {
    ped->info().resize_traits(ped->info().trait_count() + 1);

    for(RPED::RefMultiPedigree::member_const_iterator mem = ped->member_begin(); mem != ped->member_end(); ++mem)
      ped->info().set_trait(mem->index(), trait_num, getPooStatus(*mem));
  }
}

} // End namespace FUNC
} // End namespace SAGE
