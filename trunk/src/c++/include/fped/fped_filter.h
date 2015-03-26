#ifndef FPED_FILTER_H
#define FPED_FILTER_H

//==========================================================================
//  File:       fped_func.h
//                                                                          
//  Author:     Geoff Wedig
//                                                                          
//  History:    Initial implementation.                          gcw Oct 04
//                                                                          
//  Notes:      filtered_pedigree filtration functions
//                                                                          
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include <boost/call_traits.hpp>
#include "error/internal_error.h"
#include "fped/fped_obj.h"
#include "output/Output.h"

namespace SAGE
{
namespace FPED
{

/// \defgroup FPEDFilterObjects Multipedigree Filtration Objects
///
/// These objects perform basic filtering functions.  Included are:
///  - MPFilterer,
///  - FilterResults, and
///  - \ref FPEDFunctors
///
//@{

// forware declaration
class MPFilterer;
  
/// \brief Provides lists of included and excluded members after filtering a multipedigree
///      
/// The FilterResults object stores and provides access to two lists of members.
/// The first contains all members which were included during filtration of a
/// MPED data element into a FilteredMultipedigree by the MPFilterer.
/// The second contains all members excluded during the filtration. 
///
/// \internal
///
/// The FilterResults object uses a boost::shared_ptr to contain each of the
/// two lists of member pointers it contains.  This aids in copying and passing
/// the FilterResults around at an external level.  This requires internal
/// uniquifying functions, and care needs to be taken with some functions, especially
/// splicing-type operations, since the list splice requires that the lists being
/// spliced are different (for the splices we use).  See splice() for more information.
class FilterResults
{
  public:
    /// The MemberPtrList is a list of members, stored as const MPED::member_base
    /// pointers.
    typedef list<const MPED::member_base*> MemberPtrList;
    
    /// The MemberPtrIterator iterates over the member list.
    ///
    typedef MemberPtrList::const_iterator  MemberPtrIterator;
  
    /// \name Object Management
    //@{

    /// Constructor
    ///
    FilterResults();

    /// Copy Constructor
    ///
    /// \param r The FilterResults we're copying from.
    FilterResults(const FilterResults& r);

    /// Destructor
    ///
    ~FilterResults();

    /// Copy Operator
    ///
    FilterResults& operator=(const FilterResults&);
    //@}
    
    /// \name Result Information
    //@{

    /// Get the number of members included
    ///
    /// \returns Number of members included
    size_t get_included_member_count() const;

    /// Get the list of members included
    ///
    /// \returns List of members included
    const MemberPtrList& get_included_members() const;
    
    /// Get the beginning of the list of members included
    ///
    /// \returns Begin Iterator for the list of members included
    MemberPtrIterator get_included_member_begin() const;

    /// Get the ending of the list of members included
    ///
    /// \returns End Iterator for the list of members included
    MemberPtrIterator get_included_member_end  () const;
    
    /// Get the number of members excluded
    ///
    /// \returns Number of members excluded
    size_t get_excluded_member_count() const;

    /// Get the list of members excluded
    ///
    /// \returns List of members excluded
    const MemberPtrList& get_excluded_members() const;

    /// Get the beginning of the list of members excluded
    ///
    /// \returns Begin Iterator for the list of members excluded
    MemberPtrIterator get_excluded_member_begin() const;

    /// Get the ending of the list of members excluded
    ///
    /// \returns End Iterator for the list of members excluded
    MemberPtrIterator get_excluded_member_end  () const;
    
    /// Returns the excluded member list as an OUTPUT::Table.
    ///
    /// \returns The excluded member list as an OUTPUT::Table.
    OUTPUT::Table get_excluded_member_table() const;

    //@}
    
  private:
    friend class MPFilterer;
    
    typedef boost::shared_ptr<MemberPtrList> MemberPtrListShPtr;
    
    /// \name Internal Maintenance
    //@{
    
    /// Check each of our member lists, which are potentially shared, and
    /// if they are, make them unique.
    void uniquify();
      
    /// Add a member to the included list
    ///
    /// \param m The member to be added
    void add_member_to_included(const MPED::member_base& m);

    /// Add a member to the excluded list
    ///
    /// \param m The member to be added
    void add_member_to_excluded(const MPED::member_base& m);

    /// Take all the elements of the lists in \c r and splice them into our
    /// own list.  Leaves \c r empty.  If the lists in \c r are the same as
    /// those in this object (shared), then an internal error is called.  This
    /// is because we cannot splice the elements from a list into itself.
    /// This should not be a problem, as the splice shouldn't be used with
    /// lists that contain the same members.
    ///
    /// \param r The results to splice in.
    void splice(FilterResults& r);
    //@}
    
    /// The list of included members
    ///
    MemberPtrListShPtr my_included_members;
    
    /// The list of excluded members
    ///
    MemberPtrListShPtr my_excluded_members;
};

/// \brief Basic filtering of RPED::RefMultiPedigree and FilteredMultipedigree
///
/// The MPFilterer filters multipedigrees based upon some filtration
/// function.  This function is a unary function which classifies a multipedigree element as
/// included (true) or excluded (false).  Elements which cause the function to
/// return true into the user-provided FilteredPedigree.  Elements which return
/// false are excluded.
///
/// Lineage information is preserved whenever both parents and the child are all 
/// included in the new pedigree.
///
/// This is a static class that may not be instantiated.
class MPFilterer
{
  public:

    /// \name Adding Multipedigrees to FilteredMultipedigree
    ///
    /// These functions add the elements of a multipedigree to a FilteredMultipedigree.
    /// The multipedigree being added can be filtered by a test function, or each
    /// element (pedigree, subpedigree, unconnecteds or members) may be filtered individually.
    /// A functor of the appropriate type (taking an argument of the
    /// type on which filtering is performed) must be passed to the function
    /// doing the filtering.
    //@{
      
    /// Given a MPED::multipedigree-derived object, add all members to the 
    /// FilteredMultipedigree and all lineage information relating to those
    /// members
    ///
    /// \param fped The FilteredMultipedigree to copy data into
    /// \param mped The source multipedigree
    template<class MPTYPE>
    static FilterResults add_multipedigree
            (FilteredMultipedigree& fped,
             const MPTYPE&          mped);
  
    /// Given a MPED::multipedigree-derived object, add it conditional on
    /// the FILTER function returning \c true.  The FILTER must follow the interface of
    /// unary_function<MPTYPE, bool> (taking as an argument  a base class of MPTYPE is acceptable).
    ///
    /// Adds lineage information for all parent-child triples if the multipedigree is added.
    ///
    /// \param fped The FilteredMultipedigree to copy data into
    /// \param mped The source multipedigree
    /// \param f    The filter function (unary_function<MPTYPE, bool>).
    template<class MPTYPE, class FILTER>
    static FilterResults add_multipedigree_filtered
            (FilteredMultipedigree& fped,
             const MPTYPE&          mped, 
             FILTER                 f);

    /// Given a MPED::multipedigree-derived object, add each pedigree it contains
    /// conditional on the FILTER function returning \c true.  The FILTER must 
    /// follow the interface of unary_function<MPTYPE::pedigree_type, bool>
    /// (taking as an argument  a base class of MPTYPE::pedigree_type is acceptable).
    ///
    /// Adds lineage information for all parent-child triples in any pedigree which
    /// is added.
    ///
    /// \param fped The FilteredMultipedigree to copy data into
    /// \param mped The source multipedigree
    /// \param f    The filter function (unary_function<MPTYPE::pedigree_type, bool).
    template<class MPTYPE, class FILTER>
    static FilterResults add_multipedigree_filtered_by_pedigrees
            (FilteredMultipedigree& fped,
             const MPTYPE&          mped,
             FILTER                 f);

    /// Given a MPED::multipedigree-derived object, add each subpedigree it contains
    /// conditional on the FILTER function returning \c true.  The FILTER must 
    /// follow the interface of unary_function<MPTYPE::subpedigree_type, bool>
    /// (taking as an argument  a base class of MPTYPE::subpedigree_type is acceptable).
    ///
    /// Adds lineage information for all parent-child triples in any subpedigree which
    /// is added.
    ///
    /// \param fped The FilteredMultipedigree to copy data into
    /// \param mped The source multipedigree
    /// \param f    The filter function (unary_function<MPTYPE::subpedigree_type, bool).
    template<class MPTYPE, class FILTER>
    static FilterResults add_multipedigree_filtered_by_subpedigrees
            (FilteredMultipedigree& fped,
             const MPTYPE&          mped,
             FILTER                 f);

    /// Given a MPED::multipedigree-derived object, add all unconnected members for which the
    /// filter \c f returns true to the FilteredMultipedigree.  The FILTER must 
    /// follow the interface of unary_function<MPTYPE::member_type, bool>
    /// (taking as an argument  a base class of MPTYPE::member_type is acceptable).
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param mped The source multipedigree
    /// \param f    The filter function (unary_function<MPTYPE::member_type, bool).
    template<class MPTYPE, class FILTER>
    static FilterResults add_multipedigree_filtered_by_unconnecteds
            (FilteredMultipedigree& fped, 
             const MPTYPE&          mped,
             FILTER                 f);

    /// Given a MPED::multipedigree-derived object, add all members for which the
    /// filter \c f returns true to the FilteredMultipedigree.  The FILTER must 
    /// follow the interface of unary_function<MPTYPE::member_type, bool>
    /// (taking as an argument  a base class of MPTYPE::member_type is acceptable).
    ///
    /// Adds lineage information for all parent-child triples for which all three 
    /// have been added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param mped The source multipedigree
    /// \param f    The filter function (unary_function<MPTYPE::member_type, bool).
    template<class MPTYPE, class FILTER>
    static FilterResults add_multipedigree_filtered_by_members
            (FilteredMultipedigree& fped,
             const MPTYPE&          mped,
             FILTER                 f);
    
    //@}

    /// \name Adding Pedigrees to FilteredMultipedigree
    ///
    /// These functions add the elements of a pedigree to a FilteredMultipedigree.
    /// The pedigree being added can be filtered by a test function, or each
    /// element (pedigree, subpedigree, unconnecteds or members) may be filtered individually.
    /// A functor of the appropriate type (taking an argument of the
    /// type on which filtering is performed) must be passed to the function
    /// doing the filtering.
    //@{
      
    /// Given a MPED::pedigree-derived object, add all members to the 
    /// FilteredMultipedigree and all lineage information relating to those
    /// members
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param ped  The source pedigree
    template<class PTYPE>
    static FilterResults add_pedigree
            (FilteredMultipedigree& fped,
             const PTYPE&           ped);
  
    /// Given a MPED::pedigree-derived object, add it conditional on
    /// the FILTER function returning \c true.  The FILTER must follow the interface of
    /// unary_function<PTYPE, bool>
    /// (taking as an argument a base class of PTYPE is acceptable).
    ///
    /// Adds lineage information for all parent-child triples if the pedigree is added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param ped  The source pedigree
    /// \param f    The filter function (unary_function<PTYPE,bool)
    template<class PTYPE, class FILTER>
    static FilterResults add_pedigree_filtered
            (FilteredMultipedigree& fped, 
             const PTYPE&           ped,
             FILTER                 f);

    /// Given a MPED::pedigree-derived object, add each subpedigree it contains
    /// conditional on the FILTER function returning \c true.  The FILTER must 
    /// follow the interface of unary_function<PTYPE::subpedigree_type, bool>
    /// (taking as an argument a base class of PTYPE::subpedigree_type is acceptable).
    ///
    /// Adds lineage information for all parent-child triples in any subpedigree which
    /// is added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param ped  The source pedigree
    /// \param f    The filter function (unary_function<PTYPE::subpedigree_type, bool>).
    template<class PTYPE, class FILTER>
    static FilterResults add_pedigree_filtered_by_subpedigrees
            (FilteredMultipedigree& fped,
             const PTYPE&           ped, 
             FILTER                 f);

    /// Given a MPED::pedigree-derived object, add all unconnected members for which the
    /// filter \c f returns true to the FilteredMultipedigree.  The FILTER must 
    /// follow the interface of unary_function<PTYPE::member_type, bool>
    /// (taking as an argument a base class of PTYPE::member_type is acceptable).
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param ped  The source pedigree
    /// \param f    The filter function (unary_function<MPTYPE::member_type, bool>).
    template<class PTYPE, class FILTER>
    static FilterResults add_pedigree_filtered_by_unconnecteds
            (FilteredMultipedigree& fped,
             const PTYPE&           ped, 
             FILTER                 f);

    /// Given a MPED::pedigree-derived object, add all members for which the
    /// filter \c f returns true to the FilteredMultipedigree.  The FILTER must 
    /// follow the interface of unary_function<MPTYPE::member_type, bool>
    /// (taking as an argument a base class of PTYPE::member_type is acceptable).
    ///
    /// Adds lineage information for all parent-child triples for which all three 
    /// have been added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param ped  The source pedigree
    /// \param f    The filter function (unary_function<MPTYPE::member_type, bool>).
    template<class PTYPE, class FILTER>
    static FilterResults add_pedigree_filtered_by_members
            (FilteredMultipedigree& fped,
             const PTYPE&           ped,
             FILTER                 f);

    //@}

    /// \name Adding Subpedigrees to FilteredMultipedigree
    ///
    /// These functions add the elements of a subpedigree to a FilteredMultipedigree.
    /// The subpedigree being added can be filtered by a test function, or each member
    /// may be filtered individually.
    /// A functor of the appropriate type (taking an argument of the
    /// type on which filtering is performed) must be passed to the function
    /// doing the filtering.
    //@{

    /// Given a MPED::subpedigree-derived object, add all members to the 
    /// FilteredMultipedigree and all lineage information relating to those
    /// members
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param sped The source subpedigree
    template<class SPTYPE>
    static FilterResults add_subpedigree
            (FilteredMultipedigree& fped,
             const SPTYPE&          sped);
  
    /// Given a MPED::subpedigree-derived object, add it conditional on
    /// the FILTER function returning \c true.  The FILTER must follow the interface of
    /// unary_function<SPTYPE, bool>
    /// (taking as an argument a base class of SPTYPE is acceptable).
    ///
    /// Adds lineage information for all parent-child triples if the pedigree is added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param sped The source subpedigree
    /// \param f    The filter function (unary_function<SPTYPE, bool>).
    template<class SPTYPE, class FILTER>
    static FilterResults add_subpedigree_filtered
            (FilteredMultipedigree& fped, 
             const SPTYPE&          sped,
             FILTER                 f);
  
    /// Given a MPED::subpedigree-derived object, add all members for which the
    /// filter \c f returns true to the FilteredMultipedigree. The FILTER must follow the interface of
    /// unary_function<SPTYPE::member_type, bool>
    /// (taking as an argument a base class of SPTYPE::member_type is acceptable).
    ///
    /// Add lineage information
    /// for all parent-child triples for which all three have been added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param sped The source subpedigree
    /// \param f    The filter function (unary_function<SPTYPE, bool>).
    template<class SPTYPE, class FILTER>
    static FilterResults add_subpedigree_filtered_by_members
            (FilteredMultipedigree& fped,
             const SPTYPE&          sped,
             FILTER                 f);
  
    //@}

    /// \name Adding Families to FilteredMultipedigree
    ///
    /// These functions add the elements of a family to a FilteredMultipedigree.
    /// The family being added can be filtered by a test function, or each member
    /// may be filtered individually.
    /// A functor of the appropriate type (taking an argument of the
    /// type on which filtering is performed) must be passed to the function
    /// doing the filtering.
    //@{

    /// Given a MPED::family-derived object, add all members to the 
    /// FilteredMultipedigree and all lineage information relating to those
    /// members
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param fam  The source family
    template<class FTYPE>
    static FilterResults add_family
            (FilteredMultipedigree& fped,
             const FTYPE&           fam);
  
    /// Given a MPED::family-derived object, add the family conditional on the
    /// FILTER returning \c true. The FILTER must follow the interface of
    /// unary_function<FTYPE, bool>
    /// (taking as an argument a base class of FTYPE is acceptable).
    ///
    /// Lineage is added for all children.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param fam  The source family
    /// \param f    The filter function (unary_function<FTYPE, bool>).
    template<class FTYPE, class FILTER>
    static FilterResults add_family_filtered
            (FilteredMultipedigree& fped,
             const FTYPE&           fam,
             FILTER                 f);

    /// Given a MPED::family-derived object, add all members for which the
    /// filter \c f returns true to the FilteredMultipedigree. The FILTER must 
    /// follow the interface of unary_function<FTYPE::member_type, bool>
    /// (taking as an argument a base class of FTYPE::member_type is acceptable).
    ///
    /// If both parents
    /// are added, lineage is added for all children which were also added.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param fam  The source family
    /// \param f    The filter function (unary_function<FTYPE::member_type, bool>).
    template<class FTYPE, class FILTER>
    static FilterResults add_family_filtered_by_members
            (FilteredMultipedigree& fped,
             const FTYPE&           fam,
             FILTER                 f);

    //@}

    /// \name Adding Members to FilteredMultipedigree
    ///
    /// These functions add a member to a FilteredMultipedigree.
    /// The member being added can be filtered by a test function.
    /// A functor of the appropriate type (taking an argument of the
    /// type on which filtering is performed) must be passed to the function
    /// doing the filtering.
    //@{

    /// Given a MPED::member-derived object, add the member to the 
    /// FilteredMultipedigree
    ///
    /// \param fped The FilteredMultipedigree to copy the member into
    /// \param mem  The source member
    template<class MTYPE>
    static FilterResults add_member(FilteredMultipedigree& fped, const MTYPE& mem);

    /// Given a MPED::member-derived object, add the member to the 
    /// FilteredMultipedigree if the filter \c f returns true.  The FILTER must 
    /// follow the interface of unary_function<MTYPE, bool>
    /// (taking as an argument a base class of MTYPE is acceptable).
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param mem  The source member
    /// \param f    The filter function (unary_function<MTYPE, bool>).
    template<class MTYPE, class FILTER>
    static FilterResults add_member_filtered(FilteredMultipedigree& fped, const MTYPE& mem, FILTER f);
    
    //@}
    
  private:

    /// \name Object Management (restricted)
    ///@{

    /// Constructor (restricted static class)
    ///
    MPFilterer();

    /// Copy Constructor (restricted static class)
    ///
    MPFilterer(const MPFilterer&);
    ///@}
    
    
    /// \name Internal add functions.
    //@{
  
    /// Iterate over a group of members, and add the members for which the filter 
    /// \c f returns true.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param b    The begin point
    /// \param e    The end point
    /// \param f    The filter function.
    template<class Iterator>
    static FilterResults add_members(FilteredMultipedigree& fped, Iterator b, Iterator e);

    /// Iterate over a group of members, and add the members for which the filter 
    /// \c f returns true.  FILTER follows interface of unary_function<MTYPE, bool>.
    ///
    /// \param fped The FilteredMultipedigree to copy members into
    /// \param b    The begin point
    /// \param e    The end point
    /// \param f    The filter function.
    template<class Iterator, class FILTER>
    static FilterResults add_members_filtered(FilteredMultipedigree& fped, Iterator b, Iterator e, FILTER& f);
    //@}
};
  
//@}

} // End FPED namespace
} // End SAGE namespace

#include "fped/fped_filter.ipp"

#endif
