#ifndef TRAITSET_H
#define TRAITSET_H

namespace SAGE {
namespace SAMPLING {

/** \brief Stores a trait value and whether or not this trait, if missing, will be replaced with the Field's mean
  *
  * This structure is very simple: it stores the name of a trait, an individual's trait value (as a double), and a boolean indicating
  * whether or not the trait value, if missing (QNAN, that is), will be replaced with the Field's mean.
  */
struct TraitInfo
{
public:

  ///
  /// Returns \c true if getValue() is not QNAN, or if getValue() is QNAN and allowAveraging() is \c true, or
  /// if the trait was user defined.
  ///
  /// Returns \c false if getValue() is QNAN and allowAveraging() is \c false.
  bool isValid() const;

  ///
  /// The name of this trait (as it was named in the import, not in the originial RefMultiPedigree name!).
  string name;

  ///
  /// Indicates whether or not allow_averaging is enabled on this trait (see Field::ALLOW_AVERAGING).
  bool allow_averaging;

  ///
  /// Individual's value for this trait (missing indicated by \c QNAN).
  double value;

  ///
  /// Boolean value indicating whether this trait corresponds to one in a RefMultiPedigree (user_defined
  /// = false) or one created via the MemberDataSample's createField() function (user_defined = true).
  /// Please see the documentation on creating new fields in the MemberDataSample reference for more
  /// information.
  bool user_defined;
};

typedef vector<TraitInfo>::iterator       TraitInfoIterator;
typedef vector<TraitInfo>::const_iterator TraitInfoConstIterator;

/** \brief Stores information about a group of traits
  *
  * This object is simple: it stores a group's name, as well as a set of TraitInfo's comprising that group.
  */
struct GroupInfo
{
public:

  friend class MemberDataSample;
  friend class IndividualTraitData;
  
  ///
  /// The name of this trait group.
  string name;

  ///
  /// The set of TraitInfo's comprising this group.
  vector<TraitInfo> data;

  ///
  /// \internal
  /// Returns a non-const iterator to the named trait.
  /// \param trait_name The name of the trait in question
  TraitInfoConstIterator find(const string & trait_name) const;

  ///
  /// \internal
  /// Returns a non-const iterator to the beginning of this group's TraitInfo list.
  TraitInfoConstIterator begin() const;

  ///
  /// \internal
  /// Returns a non-const iterator to the end of this group's TraitInfo list.
  TraitInfoConstIterator end() const;

private:

  ///
  /// \internal
  /// Returns a non-const iterator to the named trait.
  /// \param trait_name The name of the trait in question
  TraitInfoIterator find(const string & trait_name);

  ///
  /// \internal
  /// Returns a non-const iterator to the beginning of this group's TraitInfo list.
  TraitInfoIterator begin();

  ///
  /// \internal
  /// Returns a non-const iterator to the end of this group's TraitInfo list.
  TraitInfoIterator end();
};

typedef vector<GroupInfo>::iterator       GroupInfoIterator;
typedef vector<GroupInfo>::const_iterator GroupInfoConstIterator;

/** \brief Stores information about an individual's trait values
  *
  * This object is relatively simple: it simply provides functions for iterating across the trait groups
  * that comprise the individual's trait values.
  */
struct IndividualTraitData
{
public:

  friend class MemberDataSample;

  ///
  /// Returns a non-const iterator to the named group.
  /// \param group_name The name of the group in question
  GroupInfoConstIterator find(const string & group_name) const;

  ///
  /// Returns a non-const iterator to the beginning of this object's group list.
  GroupInfoConstIterator begin() const;

  ///
  /// Returns a non-const iterator to the end of this object's group list.
  GroupInfoConstIterator end() const;

private:

  // ===== INTERNAL STUFF =====

  ///
  /// \internal
  /// Data storage of the GroupInfo's.
  vector<GroupInfo> data;

  ///
  /// \internal
  /// Returns a non-const iterator to the named group.
  /// \param group_name The name of the group in question
  GroupInfoIterator find(const string & group_name);

  ///
  /// \internal
  /// Returns a non-const iterator to the beginning of this object's group list.
  GroupInfoIterator begin();

  ///
  /// \internal
  /// Returns a non-const iterator to the end of this object's group list.
  GroupInfoIterator end();
};

//==============================
//  INLINE FUNCTIONS
//  TraitInfo
//==============================

inline bool TraitInfo::isValid() const
{
  return  (SAGE::isnan(value) == false) || 
          (SAGE::isnan(value) == true && allow_averaging == true) || 
          (user_defined == true);
}

//==============================
//  INLINE FUNCTIONS
//  GroupInfo
//==============================

inline TraitInfoConstIterator
GroupInfo::find(const string & trait_name) const
{
  for(TraitInfoConstIterator i = data.begin(); i != data.end(); ++i)
    if(i->name == trait_name)
      return i;
      
  return data.end();
}

inline TraitInfoIterator
GroupInfo::find(const string & trait_name)
{
  for(TraitInfoIterator i = data.begin(); i != data.end(); ++i)
    if(i->name == trait_name)
      return i;
      
  return data.end();
}

inline TraitInfoIterator      GroupInfo::begin ()       { return data.begin (); }
inline TraitInfoIterator      GroupInfo::end   ()       { return data.end   (); }
inline TraitInfoConstIterator GroupInfo::begin () const { return data.begin (); }
inline TraitInfoConstIterator GroupInfo::end   () const { return data.end   (); }
                    
//=============================
//  INLINE FUNCTIONS
//  IndividualTraitData
//=============================

inline GroupInfoConstIterator 
IndividualTraitData::find(const string & group_name) const
{
  for(GroupInfoConstIterator g = data.begin(); g != data.end(); ++g)
    if(g->name == group_name)
      return g;
      
  return data.end();
}

inline GroupInfoIterator 
IndividualTraitData::find(const string & group_name)
{
  for(GroupInfoIterator g = data.begin(); g != data.end(); ++g)
    if(g->name == group_name)
      return g;
      
  return data.end();
}        

inline GroupInfoIterator      IndividualTraitData::begin ()       { return data.begin (); }
inline GroupInfoIterator      IndividualTraitData::end   ()       { return data.end   (); }
inline GroupInfoConstIterator IndividualTraitData::begin () const { return data.begin (); }
inline GroupInfoConstIterator IndividualTraitData::end   () const { return data.end   (); }

} // End namespace SAMPLING
} // End namespace SAGE

#endif
