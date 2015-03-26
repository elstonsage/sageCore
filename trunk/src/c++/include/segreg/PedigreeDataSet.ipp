#ifndef PEDIGREE_DATA_SET_H
#include <segreg/PedigreeDataSet.h>
#endif

namespace SAGE {
namespace SEGREG {

/// Default Constructor
///
inline
  PedigreeDataSet::SubpedigreeIterator::SubpedigreeIterator() 
    : iterator_adaptor_((SubpedigreeList::const_iterator)0)
{ }

/// Copy Constructor
///
inline
  PedigreeDataSet::SubpedigreeIterator::SubpedigreeIterator
    (const SubpedigreeList::const_iterator& i)
    : iterator_adaptor_(i)
{ }
          
/// Dereference operation, used by boost::iterator_adaptor to dereference
/// the underlying iterator.
inline
const FPED::Subpedigree&
  PedigreeDataSet::SubpedigreeIterator::dereference() const
{
  return **this->base_reference();
}

/// Default Constructor
///
inline
  PedigreeDataSet::MemberIterator::MemberIterator() 
  : iterator_adaptor_((MemberList::const_iterator)0) 
{ }

/// Copy Constructor
///
inline
  PedigreeDataSet::MemberIterator::MemberIterator
    (const MemberList::const_iterator& i)
    : iterator_adaptor_(i)
{ }

/// Dereference operation, used by boost::iterator_adaptor to dereference
/// the underlying iterator.
inline
const FPED::Member&
  PedigreeDataSet::MemberIterator::dereference() const
{
  return **this->base_reference();
}
  
/// Default Constructor
///
inline
    PedigreeDataSet::PedigreeDataSet()
    : my_data(),
      my_subpedigrees(),
      my_members(),
      my_unconnecteds()
{ }

/// Copy Constructor
///
inline
    PedigreeDataSet::PedigreeDataSet(const PedigreeDataSet& other)
    : my_data(other.my_data),
      my_subpedigrees(other.my_subpedigrees),
      my_members(other.my_members),
      my_unconnecteds(other.my_unconnecteds)

{ }

/// Copy Operator
///
inline
PedigreeDataSet&
    PedigreeDataSet::operator=(const PedigreeDataSet& other)
{
  if(this != &other)
  {
    my_data         = other.my_data;
    my_subpedigrees = other.my_subpedigrees;
    my_members      = other.my_members;
    my_unconnecteds = other.my_unconnecteds;
  }

  return *this;
}

/// Destructor
///
inline
    PedigreeDataSet::~PedigreeDataSet()
{ }

/// Returns access to the pedigree data
///
inline
boost::shared_ptr<const FPED::Multipedigree>
    PedigreeDataSet::get_raw_data() const
{
  return my_data;
}

/// Returns \c true if there is valid data in the dataset, \c false otherwise.
///
inline
bool
    PedigreeDataSet::is_valid() const
{
  return !is_empty();
}

/// Returns \c true if there is no data in the dataset, \c false otherwise.
inline
bool
    PedigreeDataSet::is_empty() const
{
  return my_data.get() == NULL ||
         my_data->pedigree_count() == 0;
}

/// Returns an iterator to the start of the subpedigrees in the dataset
///
inline
PedigreeDataSet::SubpedigreeIterator
    PedigreeDataSet::get_subpedigree_begin() const
{
  return SubpedigreeIterator(my_subpedigrees.begin());
}

/// Returns an iterator to the end of the subpedigrees in the dataset
///
inline
PedigreeDataSet::SubpedigreeIterator
    PedigreeDataSet::get_subpedigree_end() const
{
  return SubpedigreeIterator(my_subpedigrees.end());
}

/// Returns pair of iterators to begin and end of the subpedigrees in the dataset
///
inline
PedigreeDataSet::SubpedigreeCursor
    PedigreeDataSet::get_subpedigrees() const
{
  return std::make_pair(get_subpedigree_begin(), get_subpedigree_end());
}

/// Returns number of subpedigrees in the dataset
///
inline
size_t
    PedigreeDataSet::get_subpedigree_count() const
{
  return my_subpedigrees.size();
}

/// Returns an iterator to the start of the members in the dataset
///
inline
PedigreeDataSet::MemberIterator
    PedigreeDataSet::get_member_begin() const
{
  return MemberIterator(my_members.begin());
}

/// Returns an iterator to the end of the members in the dataset
///
inline
PedigreeDataSet::MemberIterator
    PedigreeDataSet::get_member_end() const
{
  return MemberIterator(my_members.end());
}

/// Returns pair of iterators to begin and end of the members
/// in the dataset
inline
PedigreeDataSet::MemberCursor
    PedigreeDataSet::get_members() const
{
  return std::make_pair(get_member_begin(), get_member_end());
}

/// Returns number of members in the dataset
///
inline
size_t
    PedigreeDataSet::get_member_count() const
{
  return my_members.size();
}

/// Returns an iterator to the start of the unconnected members in the dataset
///
inline
PedigreeDataSet::MemberIterator
    PedigreeDataSet::get_unconnected_begin() const
{
  return MemberIterator(my_unconnecteds.begin());
}

/// Returns an iterator to the end of the unconnected members in the dataset
///
inline
PedigreeDataSet::MemberIterator
    PedigreeDataSet::get_unconnected_end() const
{
  return MemberIterator(my_unconnecteds.end());
}

/// Returns pair of iterators to begin and end of the unconnected members
/// in the dataset
inline
PedigreeDataSet::MemberCursor
    PedigreeDataSet::get_unconnecteds() const
{
  return std::make_pair(get_unconnected_begin(), get_unconnected_end());
}

/// Returns number of unconnected members in the dataset
///
inline
size_t
    PedigreeDataSet::get_unconnected_count() const
{
  return my_unconnecteds.size();
}

inline std::list<FPED::MemberConstPointer> PedigreeDataSet::getMemberList() const
{
  return my_members;
}
  

}
}

