#include "ageon/Validator.h"

namespace SAGE {
namespace AO   {

bool
Validator::isValid(size_t i, const SAMPLING::IndividualTraitData& trait_data) const
{
  SAMPLING::GroupInfoConstIterator group_info_itr  = trait_data.begin();
  for( ; group_info_itr != trait_data.end   (); ++group_info_itr )
  {
    if( group_info_itr->name == "CORE_TRAITS" )
    {
      bool AO_present   = group_info_itr->find("age of onset") ->isValid(),
           AE_present   = group_info_itr->find("age at exam")  ->isValid(),
           /* affectedness = aff_present ? group_info_itr->find("affectedness")->value : false, */
           aff_present  = group_info_itr->find("affectedness") ->isValid();

      if( !aff_present )
        return false;

      if(     !MPED::mp_utilities::is_founder(my_mp->member_index(i))
          && (!(AO_present || AE_present)) )
        return false;
    }
    else
    {
      if( !isValidGroup(trait_data, group_info_itr->name) )
        return false;
    }
  }

  //if( MPED::mp_utilities::is_founder(my_mp->member_index(i)) )
  //  return false;
                                             
  return true;
}

} // End namespace AO
} // End namespace SAGE
