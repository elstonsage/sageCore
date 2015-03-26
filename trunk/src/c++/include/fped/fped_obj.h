#ifndef FPED_OBJ_H
#define FPED_OBJ_H

//==========================================================================
//  File:       fped.h
//                                                                          
//  Author:     Geoff Wedig & Yeunjoo Song
//                                                                          
//  History:    Initial implementation.                          yjs Nov. 03
//              Major Revision to make general, documentation    gcw Oct. 04 
//                                                                          
//  Notes:      filtered_pedigree objects.
//                                                                          
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "error/internal_error.h"
#include "mped/mp.h"
#include "rped/rped.h"

namespace SAGE {
namespace FPED {

// Forward declarations.  Descriptions to follow with classes.
  
class FilteredMemberInfo;
class FilteredPedigreeInfo;
class FilteredMultipedigreeInfo;

/// \defgroup FPEDPedObjects Filtered Multipedigree Objects
///
/// The Filtered Multipedigree objects are the data structures which provide
/// the access to the data once it has been filtered.  The primary class here
/// is the FilteredMultipedigree, but also included are the info classes which are
/// needed to make it function.
///
//@{

/// \brief Stores multipedigree info for the FilteredMultipedigree
///
/// The FilteredMultipedigreeInfo stores multipedigree info for the FilteredMultipedigree.
/// Primarily, this is the same information storied in the RPED::RefMPedInfo.  The
/// only additional information is the source RPED::RefMultipedigree of the 
/// FilteredMultipedigree.
class FilteredMultipedigreeInfo : public RPED::RefMPedInfo
{
  public:
  
    friend class FilteredMultipedigree;
  
    /// Default Constructor.
    FilteredMultipedigreeInfo();
    
    /// Returns the RPED::Multipedigree source of the FPED's data (members, pedigrees,
    /// and so on).
    const RPED::MultiPedigree* get_source_rped() const;
  
  private:
    
    FilteredMultipedigreeInfo(const RPED::RefMPedInfo& rped);
  
    /// Sets the RPED::Multipedigree source of the FPED's data (members, pedigrees,
    /// and so on).
    void set_source_rped(const RPED::MultiPedigree* rp);
    
    const RPED::MultiPedigree* my_rped_source; ///< Storage of the RPED::Multipedigree source
};    

/// \brief Class for storing filtered pedigree data
///
/// The FilteredMultipedigree is the primary class of the FPED library.  It
/// is a SAGE::MPED::multipedigree derived type which allows for the construction
/// of multipedigrees of arbitrary structure based upon an underlying
/// SAGE::RPED::RefMultiPedigree.
class FilteredMultipedigree : 
    public MPED::multipedigree<FilteredMemberInfo, MPED::no_info, MPED::no_info,
                               FilteredPedigreeInfo, FilteredMultipedigreeInfo>
{
  public:
  
    /// Constructor.  Given a source RPED::RefMultiPedigree, creates the 
    /// FilteredMultipedigree based upon that source.  Note that the new 
    /// FilteredMultipedigree is left clean, no individual data is copied, 
    /// just the RPED::RefMPedInfo of the source as well as a pointer to it for 
    /// future reference.
    ///
    /// \param source The source RPED::RefMultiPedigree.
    inline explicit FilteredMultipedigree(const RPED::RefMultiPedigree& source);

    /// Constructor.  Given a source FilteredMultipedigree, creates the 
    /// FilteredMultipedigree based upon that source.  Note that the new 
    /// FilteredMultipedigree is left clean, no individual data is copied, 
    /// just the RPED::RefMPedInfo of the source as well as a pointer to it for 
    /// future reference, so this is \em not a copy constructor.
    ///
    /// \param source The source FilteredMultipedigree.
    inline explicit FilteredMultipedigree(const FilteredMultipedigree& source);

    /// Destructor
    ///
    inline ~FilteredMultipedigree();

    /// This function is the last thing that should be done to a
    /// FilteredMultipedigree before it is ready to be used.  It
    /// builds the multipedigree, populates the FilteredPedigreeInfo
    /// and fixes the FilteredMemberInfo's to the right pedigree infos
    void construct();

  private:

    /// Does a partial sort on the filtered_multipedigree, such that no
    /// child has a smaller pedigree or subpedigree index than their parents
    void sort_into_descent_order();

    /// Does a partial sort on the pedigree, such that no
    /// child has a smaller pedigree index than their parents
    ///
    /// \param p The Pedigree to sort
    void sort_ped_into_descent_order (pedigree_type&     p);

    /// Does a partial sort on the subpedigree, such that no
    /// child has a smaller subpedigree index than their parents
    ///
    /// \param sp The Subpedigree to sort
    void sort_sped_into_descent_order(subpedigree_type& sp);
    
    /// Default Constructor is private to prevent bad instantiation
    ///
    //lint -e{1704}
    inline FilteredMultipedigree();

    /// Copy operator is private to prevent bad instantiation
    ///
    FilteredMultipedigree& operator=(const FilteredMultipedigree&);

    /// Storage of the original source.  This is useful for making sure
    /// any additions came from the correct place.
    const RPED::RefMultiPedigree* my_source;
};


/// \name Public Typedefs
///
/// These typedefs define useful bits of the FilteredMultipedigree.  They
/// are equivalent to similar typedefs for the RPED::RefMultiPedigree.
//@{
typedef FilteredMultipedigree                             Multipedigree;
typedef FilteredMultipedigree::pedigree_type              Pedigree;
typedef FilteredMultipedigree::subpedigree_type           Subpedigree;
typedef FilteredMultipedigree::family_type                Family;
typedef FilteredMultipedigree::member_type                Member;

typedef FilteredMultipedigree::pedigree_iterator          PedigreeIterator;
typedef FilteredMultipedigree::subpedigree_iterator       SubpedigreeIterator;
typedef FilteredMultipedigree::family_iterator            FamilyIterator;
typedef FilteredMultipedigree::member_iterator            MemberIterator;
typedef FilteredMultipedigree::offspring_iterator         OffspringIterator;
typedef FilteredMultipedigree::progeny_iterator           ProgenyIterator;
typedef FilteredMultipedigree::sibling_iterator           SiblingIterator;
typedef FilteredMultipedigree::mate_iterator              MateIterator;

typedef FilteredMultipedigree::pedigree_const_iterator    PedigreeConstIterator;
typedef FilteredMultipedigree::subpedigree_const_iterator SubpedigreeConstIterator;
typedef FilteredMultipedigree::family_const_iterator      FamilyConstIterator;
typedef FilteredMultipedigree::member_const_iterator      MemberConstIterator;
typedef FilteredMultipedigree::offspring_const_iterator   OffspringConstIterator;
typedef FilteredMultipedigree::progeny_const_iterator     ProgenyConstIterator;
typedef FilteredMultipedigree::sibling_const_iterator     SiblingConstIterator;
typedef FilteredMultipedigree::mate_const_iterator        MateConstIterator;

typedef FilteredMultipedigree*                            MultipedigreePointer;
typedef FilteredMultipedigree::pedigree_pointer           PedigreePointer;
typedef FilteredMultipedigree::subpedigree_pointer        SubpedigreePointer;
typedef FilteredMultipedigree::family_pointer             FamilyPointer;
typedef FilteredMultipedigree::member_pointer             MemberPointer;

typedef const FilteredMultipedigree*                      MultipedigreeConstPointer;
typedef FilteredMultipedigree::pedigree_const_pointer     PedigreeConstPointer;
typedef FilteredMultipedigree::subpedigree_const_pointer  SubpedigreeConstPointer;
typedef FilteredMultipedigree::family_const_pointer       FamilyConstPointer;
typedef FilteredMultipedigree::member_const_pointer       MemberConstPointer;
//@}

/// \brief Stores a reference to member data for the FilteredMultipedigree
///
/// The FilteredMemberInfo provides access to member data (traits, phenotypes, etc)
/// of a member in an underlying RPED::RefMultiPedigree.  The interface is designed to
/// mimic that of the RPED::RefPedInfo, with the exclusion of the member's index, 
/// which isn't required. Only the const interface is replicated.

class FilteredMemberInfo
{
    friend class FilteredMultipedigree;
    friend class FilteredPedigreeInfo;

  public:

    /// \name Object Management (Constructors, destructors, copying)
    //@{
    /// Default Constructor
    ///
    FilteredMemberInfo();
    
    /// Create a FilteredMemberInfo based upon the member given
    ///
    /// \param member The member to which we link
    FilteredMemberInfo(const RPED::RefMember& member);

    /// Create a FilteredMemberInfo based upon the member given
    ///
    /// \param member The member which contains the member to which we link
    FilteredMemberInfo(const Member& member);

    /// Copy Constructor
    ///
    /// \param info The FilteredMemberInfo we're copying
    FilteredMemberInfo(const FilteredMemberInfo& info);

    /// Copy Operator
    ///
    /// \param info The FilteredMemberInfo we're copying
    FilteredMemberInfo& operator=(const FilteredMemberInfo& info);
    
    /// Destructor
    ///
    ~FilteredMemberInfo();
    //@}

    /// \name Source Member Access
    //@{
      
    /// Sets the source member to the member given
    ///
    /// \param member The member we want set to
    void set_source_member(const RPED::RefMember& member);
    
    /// Sets the source member to the source member in the filtered_member given
    ///
    /// \param member The filtered_member which contains our new source member
    void set_source_member(const Member& member);
    
    /// Returns the RPED::RefMember to which we currently refer.  Returns null if we
    /// don't point to any valid RPED::RefMember.
    const RPED::RefMember* get_source_member() const;
    //@}
    
    /// \name Data Access
    //@{
    
    /// Return the number of traits
    ///
    size_t trait_count()     const;

    /// Return the number of strings
    ///
    size_t string_count()    const;
    
    /// Return the number of markers
    size_t marker_count()    const;

    /// Return a particular trait value
    ///
    /// \param t The trait we want
    double  trait             (size_t t) const;
    
    /// Returns \c true if the trait is missing or otherwise unavailable,
    /// \c false otherwise
    ///
    /// \param t The trait we want to check
    bool    trait_missing     (size_t t) const;

    /// Return a particular string value
    ///
    /// \param s The string we want
    string  get_string        (size_t s) const;

    /// Return a particular phenotype value
    ///
    /// \param m The marker id of the marker whose phenotype we want
    uint    phenotype         (size_t m) const;

    /// Returns \c true if the phenotype is missing or otherwise unavailable,
    /// \c false otherwise
    ///
    /// \param m  The marker id of the marker whose phenotype whose we want to check
    /// \param mi The marker information regarding marker \c m.
    bool    phenotype_missing (size_t m, const RPED::RefMarkerInfo& mi) const;
    
    //@}

  private:
  
    const RPED::RefMember*    my_member;       ///< The member we're based upon 
    const RPED::RefPedInfo*   my_ref_pedinfo;  ///< The RPED::RefPedInfo where our data is stored
    size_t              my_ref_index;    ///< The index into the RPED::RefPedInfor where our data is stored
};

/// \brief Stores a reference to pedigree data for the FilteredMultipedigree
///
/// The FilteredPedigreeInfo provides access to the pedigree info
/// of a pedigree in an underlying RPED::RefMultiPedigree.  The interface is designed to
/// mimic that of the RPED::RefPedInfo. Only the const interface is replicated.
class FilteredPedigreeInfo
{
  public:

    /// \name Object Management (Constructors, destructors, copying)
    //@{
    /// Default Constructor
    ///
    FilteredPedigreeInfo();
    
    /// Copy Constructor
    ///
    /// \param info The FilteredMemberInfo we're copying
    FilteredPedigreeInfo(const FilteredPedigreeInfo& info);

    /// Copy Operator
    ///
    /// \param info The FilteredMemberInfo we're copying
    FilteredPedigreeInfo& operator=(const FilteredPedigreeInfo& info);
    
    /// Destructor
    ///
    ~FilteredPedigreeInfo();
    //@}
    
    /// \name Data Access
    //@{
    
    /// Return the number of traits
    ///
    size_t trait_count()     const;
    
    /// Return the number of strings
    ///
    size_t string_count()    const;
    
    /// Return the number of markers
    ///
    size_t marker_count()    const;
    
    /// Return the number of members
    ///
    size_t member_count()    const;

    /// Return a particular trait value
    ///
    /// \param i The index of the member whose value we want
    /// \param t The trait we want
    double  trait             (size_t i, size_t t) const;

    /// Returns \c true if the trait is missing or otherwise unavailable,
    /// \c false otherwise
    ///
    /// \param i The index of the member whose value we want to check
    /// \param t The trait we want to check
    bool    trait_missing     (size_t i, size_t t) const;

    /// Return a particular string value
    ///
    /// \param i The index of the member whose value we want
    /// \param s The string we want
    string  get_string        (size_t i, size_t s) const;

    /// Return a particular phenotype value
    ///
    /// \param i The index of the member whose value we want
    /// \param m The marker id of the marker whose phenotype we want
    uint    phenotype         (size_t i, size_t m) const;

    /// Returns \c true if the phenotype is missing or otherwise unavailable,
    /// \c false otherwise
    ///
    /// \param i  The index of the member whose value we want to check
    /// \param m  The marker id of the marker whose phenotype whose we want to check
    /// \param mi The marker information regarding marker \c m.
    bool    phenotype_missing (size_t i, size_t m, const RPED::RefMarkerInfo& mi) const;

    /// Switches the indices of m1 and m2.  This needs doing any time members
    /// are swapped in the filtered_pedigree to which this belongs.
    void    swap_members(size_t m1, size_t m2);

    /// Dumps the current indices of all members to cout.  Should only be used
    /// for debugging purposes.
    void    dump_ref_indices() const;
    
    /// Returns the source pedigree from the RPED
    ///
    const RPED::Pedigree* get_source_pedigree() const;

  private:
    
    friend class FilteredMultipedigree;
  
    /// Constructs the info class based upon the fped (this is the fped to which
    /// this filter belongs) and the RPED::RefPedigree which is its source
    ///
    /// \param fped    The filtered_pedigree to which this belongs.
    /// \param ref_ped The source RPED::RefPedigree
    void construct(Pedigree& fped, const RPED::RefPedigree* ref_ped);

    const RPED::RefPedigree*  my_ped_source;  ///< The pedigree which is our source
    const RPED::RefPedInfo*   my_ref_pedinfo; ///< The info of our source pedigree
    vector<size_t>            my_ref_indices; ///< The indices of our members in the source pedigree's info object
};
//@}

} // end of namespace RPED
} // end of namespace SAGE

#include "fped/fped_obj.ipp"

#endif
