#ifndef INDIVIDUAL_VALIDATOR_H
#define INDIVIDUAL_VALIDATOR_H

#include "sampling/TraitSet.h"

namespace SAGE {
namespace SAMPLING {

/** \internal
  * \brief Pure virtual validator parent class.
  *
  * Invalidates any individual with \b any missing values.
  */
class IndividualValidator
{
public:

  // Required to make compiler happy.
  virtual ~IndividualValidator() { }

  ///
  /// This is the core virtual interface for any object derived from the IndividualValidator.
  /// Given a Member index (mpindex()), and a structure representing that member's trait data, this 
  /// function evaluates if the member should be considered valid. It returns \c true if the
  /// member \b is valid, \c false otherwise.
  ///
  /// See detailed class description for more information.
  /// \param i The individual in question
  /// \param trait_data The trait values associated with this individual
  virtual bool isValid (size_t i, const IndividualTraitData& trait_data) const = 0;

protected:
  ///
  /// Given a field group name, returns \c true if all the group's contituent fields are valid
  /// (for the given individual). If any single trait value is not valid, this function returns false.
  /// \param i The individual in question
  /// \param trait_data The trait valeus assocaited with this individual
  /// \param group_name The name of the field group
  bool isValidGroup(const IndividualTraitData& trait_data, const string& group_name) const;
};

class DefaultValidator : public IndividualValidator
{
public:
  virtual inline bool isValid (size_t i, const IndividualTraitData& trait_data) const;
};

//=====================
//  INLINE FUNCTIONS
//  IndividualValidator
//=====================

inline bool 
IndividualValidator::isValidGroup(
	const IndividualTraitData& trait_data, 
	const string&              group_name) const
{
  GroupInfoConstIterator group = trait_data.find(group_name);
  
  for(TraitInfoConstIterator trait_info_itr  = group->begin ();
                             trait_info_itr != group->end   (); ++trait_info_itr)
  {
   if(trait_info_itr->isValid() == false)
       return false;
  }
  
  return true;
}
  
//=====================
//  INLINE FUNCTIONS
//  DefaultValidator
//=====================

inline bool
DefaultValidator::isValid(size_t, const IndividualTraitData & trait_data) const
{
  for(GroupInfoConstIterator group_info_itr  = trait_data.begin ();
                             group_info_itr != trait_data.end   (); ++group_info_itr)
  {
    if(!isValidGroup(trait_data, group_info_itr->name))
      return false;
  }
  
  return true;
}


} // End namespace SAMPLING
} // End namespace SAGE

#endif
