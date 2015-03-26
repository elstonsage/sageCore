#ifndef FPED_H
#include "fped/fped.h"
#endif

namespace SAGE {
namespace FPED {

//lint --e{613}  Lots of spurious 'null pointer' errors with my_ped_info, but 
//               we check it everywhere before calling functions on it.

//------------------------------------------------------------------
//  Inline implementation of filtered_mmultipedigree
//------------------------------------------------------------------
inline
FilteredMultipedigree::FilteredMultipedigree(const RPED::RefMultiPedigree& source)
  : MPED::multipedigree<FilteredMemberInfo, MPED::no_info, MPED::no_info,
                        FilteredPedigreeInfo, FilteredMultipedigreeInfo>(),
    my_source(&source)
{
  // Make our info a copy of the info provided by source.

  info() = source.info();
  
  info().set_source_rped(&source);
}

//lint -e{1738} Not really a copy constructor
//lint -e{1554} Pointer copy of my_source is ok
inline
FilteredMultipedigree::FilteredMultipedigree(const FilteredMultipedigree& source)
  : MPED::multipedigree<FilteredMemberInfo, MPED::no_info, MPED::no_info,
                        FilteredPedigreeInfo, FilteredMultipedigreeInfo>(),
    my_source(source.my_source)
{
  // Make our info a copy of the info provided by source.

  info() = source.info();

  info().set_source_rped(source.info().get_source_rped());
}
inline
FilteredMultipedigree::~FilteredMultipedigree()
{
  my_source = NULL;
}

inline
FilteredMultipedigree::FilteredMultipedigree()
  : my_source(NULL)
{ }

// Copying is disabled (made private).  This minimal version is here to make 
// compilers and linkers not complain.
//lint -e{1529,1745}  This isn't a real copy constructor, so we have to disable lint
inline
FilteredMultipedigree& FilteredMultipedigree::operator=(const FilteredMultipedigree&)
{ 
  return *this;
}

//------------------------------------------------------------------
//  Inline implementation of FilteredMemberInfo
//------------------------------------------------------------------

inline
FilteredMemberInfo::FilteredMemberInfo()
  : my_member     (NULL),
    my_ref_pedinfo(NULL),
    my_ref_index  ((size_t) -1)
{ }

inline
FilteredMemberInfo::FilteredMemberInfo(const RPED::RefMember& member)
  : my_member     (&member),
    my_ref_pedinfo(&member.pedigree()->info()),
    my_ref_index  (member.index())
{ }

inline
FilteredMemberInfo::FilteredMemberInfo(const Member& member)
  : my_member     (member.info().my_member),
    my_ref_pedinfo(member.info().my_ref_pedinfo),
    my_ref_index  (member.info().my_ref_index)
{ }

//lint -e{1554} Direct copy of pointer ok.
inline
FilteredMemberInfo::FilteredMemberInfo(const FilteredMemberInfo& info)
  : my_member     (info.my_member),
    my_ref_pedinfo(info.my_ref_pedinfo),
    my_ref_index  (info.my_ref_index)
{ }

//lint -e{1555} Direct copy of pointer ok.
inline
FilteredMemberInfo& FilteredMemberInfo::operator=
    (const FilteredMemberInfo& info)
{
  if(this == &info) return *this;
  
  my_member      = info.my_member;
  my_ref_pedinfo = info.my_ref_pedinfo;
  my_ref_index   = info.my_ref_index;
  
  return *this;
}

inline
FilteredMemberInfo::~FilteredMemberInfo()
{
  my_member      = NULL;
  my_ref_pedinfo = NULL;
}

inline
void FilteredMemberInfo::set_source_member(const RPED::RefMember& member)
{
  my_member      = &member;
  my_ref_pedinfo = &member.pedigree()->info();
  my_ref_index   = member.index();
}

inline
void FilteredMemberInfo::set_source_member(const Member& member)
{
  *this = member.info();
}

inline
const RPED::RefMember* FilteredMemberInfo::get_source_member() const
{
  return my_member;
}

inline size_t
FilteredMemberInfo::trait_count() const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  return my_ref_pedinfo->trait_count();
}

inline size_t
FilteredMemberInfo::string_count() const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  return my_ref_pedinfo->string_count();
}

inline size_t
FilteredMemberInfo::marker_count() const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  return my_ref_pedinfo->marker_count();
}

inline double
FilteredMemberInfo::trait(size_t t) const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  if( my_ref_index >= my_ref_pedinfo->member_count() )
    return std::numeric_limits<double>::quiet_NaN();

  return my_ref_pedinfo->trait(my_ref_index, t);
}

inline bool
FilteredMemberInfo::trait_missing(size_t t) const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  if( my_ref_index >= my_ref_pedinfo->member_count() )
    return true;

  return my_ref_pedinfo->trait_missing(my_ref_index, t);
}

inline string
FilteredMemberInfo::get_string(size_t s) const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  if( my_ref_index >= my_ref_pedinfo->member_count() )
    return "";

  return my_ref_pedinfo->get_string(my_ref_index, s);
}

inline uint
FilteredMemberInfo::phenotype(size_t m) const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  if( my_ref_index >= my_ref_pedinfo->member_count() )
    return MLOCUS::NPOS;

  return my_ref_pedinfo->phenotype(my_ref_index, m);
}

inline bool
FilteredMemberInfo::phenotype_missing(size_t m, const RPED::RefMarkerInfo& mi) const
{
  // This function should never be called if the FilteredMemberInfo is
  // empty, so this indicates something seriously wrong.  We call an internal_error
  if(!my_ref_pedinfo)
    SAGE_internal_error();

  if( my_ref_index >= my_ref_pedinfo->member_count() )
    return true;

  return my_ref_pedinfo->phenotype_missing(my_ref_index, m, mi);
}

//------------------------------------------------------------------
//  Inline implementation of FilteredPedigreeInfo
//------------------------------------------------------------------
inline
FilteredPedigreeInfo::FilteredPedigreeInfo()
  : my_ped_source (NULL),
    my_ref_pedinfo(NULL),
    my_ref_indices()
{ }

//lint -e{1554} Direct copy of pointer ok.
inline
FilteredPedigreeInfo::FilteredPedigreeInfo(const FilteredPedigreeInfo& info)
  : my_ped_source (info.my_ped_source),
    my_ref_pedinfo(info.my_ref_pedinfo),
    my_ref_indices(info.my_ref_indices)
{ }

//lint -e{1555} Direct copy of pointer ok.
inline
FilteredPedigreeInfo& FilteredPedigreeInfo::operator=(const FilteredPedigreeInfo& info)
{
  if(this == &info) return *this;
  
  my_ped_source   = info.my_ped_source;
  my_ref_pedinfo  = info.my_ref_pedinfo;
  my_ref_indices  = info.my_ref_indices;

  return *this;                  
}

inline
FilteredPedigreeInfo::~FilteredPedigreeInfo()
{
  my_ped_source  = NULL;
  my_ref_pedinfo = NULL;
}

inline size_t
FilteredPedigreeInfo::trait_count() const
{
  return my_ref_pedinfo->trait_count();
}

inline size_t
FilteredPedigreeInfo::string_count() const
{
  return my_ref_pedinfo->string_count();
}

inline size_t
FilteredPedigreeInfo::marker_count() const
{
  return my_ref_pedinfo->marker_count();
}

inline size_t
FilteredPedigreeInfo::member_count() const
{
  return my_ref_indices.size();
}

inline double
FilteredPedigreeInfo::trait(size_t i, size_t t) const
{
  if( i >= my_ref_indices.size() )
    return std::numeric_limits<double>::quiet_NaN();
 
  if( my_ref_indices[i] >= my_ref_pedinfo->member_count() )
    return std::numeric_limits<double>::quiet_NaN();

  return my_ref_pedinfo->trait(my_ref_indices[i], t);
}

inline bool
FilteredPedigreeInfo::trait_missing(size_t i, size_t t) const
{
  if( i >= my_ref_indices.size() )
    return true;

  if( my_ref_indices[i] >= my_ref_pedinfo->member_count() )
    return true;

  return my_ref_pedinfo->trait_missing(my_ref_indices[i], t);
}

inline string
FilteredPedigreeInfo::get_string(size_t i, size_t s) const
{
  if( i >= my_ref_indices.size() )
    return "";

  if( my_ref_indices[i] >= my_ref_pedinfo->member_count() )
    return "";

  return my_ref_pedinfo->get_string(my_ref_indices[i], s);
}

inline uint
FilteredPedigreeInfo::phenotype(size_t i, size_t m) const
{
  if( i >= my_ref_indices.size() )
    return MLOCUS::NPOS;

  if( my_ref_indices[i] >= my_ref_pedinfo->member_count() )
    return MLOCUS::NPOS;

  return my_ref_pedinfo->phenotype(my_ref_indices[i], m);
}

inline bool
FilteredPedigreeInfo::phenotype_missing(size_t i, size_t m, const RPED::RefMarkerInfo& mi) const
{
  if( i >= my_ref_indices.size() )
    return true;

  if( my_ref_indices[i] >= my_ref_pedinfo->member_count() )
    return true;

  return my_ref_pedinfo->phenotype_missing(my_ref_indices[i], m, mi);
}

inline const RPED::Pedigree* 
FilteredPedigreeInfo::get_source_pedigree() const
{
  return my_ped_source;
}

inline
FilteredMultipedigreeInfo::FilteredMultipedigreeInfo()
  : my_rped_source(NULL)
{ }
    
inline const RPED::MultiPedigree* 
FilteredMultipedigreeInfo::get_source_rped() const
{
  return my_rped_source;
}
inline
FilteredMultipedigreeInfo::FilteredMultipedigreeInfo(const RPED::RefMPedInfo& rped)
  : RPED::RefMPedInfo(rped),
    my_rped_source(NULL)
{ }
  
inline void 
FilteredMultipedigreeInfo::set_source_rped(const RPED::MultiPedigree* rp)
{
  my_rped_source = rp;
}

} // End namespace RPED
} // End namespace SAGE
