#ifndef FPED_FUNC_H
#define FPED_FUNC_H

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

namespace SAGE
{
namespace FPED
{
  
/// \defgroup FPEDFunctors Multipedigree Filtration Functors
/// 
/// These functor objects are used by the MPFilterer to provide common 
/// options for filtering multipedigrees and their parts.
///
/// Filtering functions should be of the type:
///
/// \code
///
/// template<class TYPE>
/// class filter_function : public unary_function<TYPE, bool>
/// {
///   public:
///     // This function should return true if m is to be included in the
///     // filtered subset, false otherwise 
///     bool operator()(const TYPE& m) const;
/// };
///
/// \endcode
///
/// This code can be copied as a starting point for new filter functors.
///
/// \par Notes
///
/// -# These functors use the functional naming convention rather than the
///    class naming convention since they are used as functions, rather than
///    classes.  One minor modification of this method is used for \ref METAFILTERS "meta-filters"
///    is the \c _t appended to the name to allow actual build functions to
///    use that name.  See the meta-filters section for more information. 
/// -# Filtration functions are templatized to allow them to work on various
///    multipedigree types.  The boolean operator is assumed to get a
///    specific type (multipedigree, pedigree, subpedigree, family or member) reference 
///    which can be tested and whose members (if not itself a member reference)
///    have access to a RefMember which will be the source of the FilteredMember. 
///    Currently, the only two types which provide this functionality are the
///    RefMember and the FilteredMember, but more can be added in the future.
///    If filtration can only be done from a specific type (some custom multipedigree
///    derived type), there is no need for templatization, but this should be avoided.
/// -# The unary function derivation adds typedefs which allow certain STL algorithms
///    to work.  In many cases, it can be dropped (see always_keep), but
///    that should be avoided where possible.  See STL documentation for more information.
/// -# Functors are passed by value (copied) to the MPFilterer.  This is
///    consistent with the way functors work in the STL.  For this reason, the
///    following restrictions should be observed:
///    -# Functors should be small (minimal internal state)
///    -# Functors should be monomorphic (no virtual functions)
///    -# Copy constructors and operators should be provided if a functor has
///       internal state.
///
/// There are two categories of filter functors, basic filters and meta-filters.
/// Basic filters do a simple, fixed test of the type given to them.  Examples
/// include always_keep, and has_informative_loci. Meta-filters use
/// other filters as their test, but extend that filter's functionality in some
/// pre-defined way.  Examples include is_inf_within_sped_t and is_family_sib_pair_inf_t.
///
/// \anchor METAFILTERS
///
/// \par Meta-Filters
///
/// A few notes about meta-filters.  Since their types can get quite unweildy due
/// to template instantiations, they use a slightly different convention.  First,
/// they have an _t appended to their name, indicating the type of them.  The
/// name without the _t is then used for a template function which creates the
/// appropriate type.  This is similar to the std::make_pair() and other
/// such templatized functions in the STL, and therefore follows that convention.
///
/// \ingroup FPEDFilterObjects
//@{
  
/// \brief Functor which given any type, always returns true
///
/// This function can be used to add specific parts of the source multipedigree
/// unconditionally.
///
/// Note, that since it is templatized at function evaluation time, this is
/// not an adaptable STL function, and won't be useable in certain contexts (but
/// seems to be fine within FPED).
class always_keep
{
  public:
  
    template<class TYPE>
    bool operator()(const TYPE&) const;
};

/// \brief Functor which given a member, returns true if any of the loci we're checking are informative
///
/// The has_informative_loci function stores a set of loci to check.  Upon
/// being given a member, it checks that member for informativity at those
/// loci.  If even one of those loci is informative, it returns true.  Otherwise,
/// it returns false.
template<class MTYPE>
class has_informative_loci : public unary_function<MTYPE, bool>
{
  public:
  
    //typedef typename MTYPE::multipedigree_type multipedigree_type;
    //typename MTYPE::subpedigree_type   subpedigree_type;
  
    /// Constructor.
    ///
    /// The primary constructor requires the multipedigree to which the members
    /// we're filtering belong.  This is to get access to the loci it has stored.
    ///
    /// \param mped    The multipedigree we'll be filtering from
    /// \param use_all Default status of the check_status flag for all markers.
    ///                Can be set to true if all markers are to be checked (default)
    ///                or false if individual markers need to be set.
    has_informative_loci(const typename MTYPE::multipedigree_type& mped, bool use_all = true);
    
    /// Copy Constructor
    ///
    has_informative_loci(const has_informative_loci&);

    /// Destructor
    ///
    //lint -e{1510} unary_function has no destructor, but that's no big thing
    ~has_informative_loci();
    
    /// Copy Operator
    ///
    has_informative_loci& operator=(const has_informative_loci&);
    
    /// Sets the check status for a particular locus
    ///
    /// \param locus_index The index of the locus we're setting
    /// \param check       Boolean stating whether we should check this locus.
    void set_check_status_for_locus(size_t locus_index, bool check);
    
    /// Gets the check status for a particular locus
    ///
    /// \param locus_index The index of the locus we're setting
    /// \returns \c true if we are checking that locus, \c false otherwise
    bool get_check_status_for_locus(size_t locus_index) const;
    
    /// Our test function.
    ///
    /// \internal
    /// This function is in the ipp file because it is templatized and therefore
    /// cannot be in the cpp, but it is not inline due to the loop it must run to
    /// check all loci.
    ///
    /// \param member The member to check
    /// \returns \c true if the member has at least one phenotype which is 
    ///          not missing or unavailable, \c false otherwise.
    bool operator()(const MTYPE& member) const;
    
  private:

    /// Default constructor is private to prevent default creation.
    ///
    //lint -e{1704}
    has_informative_loci();

    /// A vector keeping track of which loci we are to check upon being queried.
    ///
    /// While there may be many loci, boolean vectors are very space efficient, so
    /// this is ok for pass-by-value operations.
    vector<bool> my_loci_check_states;
};

/// \brief Functor which given a member, returns true if any of the traits we're checking are informative
///
/// The has_informative_traits function stores a set of traits to check.  Upon
/// being given a member, it checks that member for informativity at those
/// traits.  If even one of those loci is informative, it returns true.  Otherwise,
/// it returns false.
template<class MTYPE>
class has_informative_traits : public unary_function<MTYPE, bool>
{
  public:
  
    /// Constructor.
    ///
    /// The primary constructor requires the multipedigree to which the members
    /// we're filtering belong.  This is to get access to the traits it has stored.
    ///
    /// \param mped    The multipedigree we'll be filtering from
    /// \param use_all Default status of the check_status flag for all markers.
    ///                Can be set to true if all markers are to be checked (default)
    ///                or false if individual markers need to be set.
    has_informative_traits(const typename MTYPE::multipedigree_type& mped, bool use_all = true);
    
    /// Copy Constructor
    ///
    has_informative_traits(const has_informative_traits&);

    /// Destructor
    ///
    //lint -e{1510} unary_function has no destructor, but that's no big thing
    ~has_informative_traits();
    
    /// Copy Operator
    ///
    has_informative_traits& operator=(const has_informative_traits&);
    
    /// Sets the check status for a particular trait
    ///
    /// \param trait_index The index of the trait we're setting
    /// \param check       Boolean stating whether we should check this trait.
    void set_check_status_for_trait(size_t trait_index, bool check);
    
    /// Gets the check status for a particular trait
    ///
    /// \param trait_index The index of the trait we're setting
    /// \returns \c true if we are checking that trait, \c false otherwise
    bool get_check_status_for_trait(size_t trait_index) const;
    
    /// Our test function.
    ///
    /// \internal
    /// This function is in the ipp file because it is templatized and therefore
    /// cannot be in the cpp, but it is not inline due to the loop it must run to
    /// check all loci.
    ///
    /// \param member The member to check
    /// \returns \c true if the member has at least one trait which is 
    ///          not missing or unavailable, \c false otherwise.
    bool operator()(const MTYPE& member) const;
    
  private:

    /// Default constructor is private to prevent default creation.
    ///
    //lint -e{1704}
    has_informative_traits();

    /// A vector keeping track of which traits we are to check upon being queried.
    ///
    /// While there may be many traits, boolean vectors are very space efficient, so
    /// this is ok for pass-by-value operations.
    vector<bool> my_trait_check_states;
};

/// \brief Functor which, given a member, returns true if they're informative within the context of their subpedigree.
///
/// This filter is a meta-filter.  It takes another filter as an argument
/// at construction, and uses that as its base informativity.
///
/// It calculates, for each individual, whether they are informative within
/// their subpedigree (unconnected individuals are assumed to be in a subpedigree of
/// size 1).  An individual is informative within the subpedigree if:
///
/// -# It is informative by the filter function provided, or
/// -# It helps connect two members who are informative by the provided function
///    through a common ancestor
///
/// A note about connection:  A member \em i helps connect two other members \em j 
/// and \em k through a common ancestor if and only if \em j and \em k have a common
/// ancestor, and \em i is necessary to preserve the relationship between that
/// ancestor and \em j and \em k.  This includes the members on the direct path between
/// the ancestor and \em j and \em k as well as those spouses needed to maintain
/// lineage (since a member may only have 0 or 2 parents)
///
/// Like other meta-filters, this has the _t appended to the end to allow for the
/// is_inf_within_sped() function to use that name.  This simplifies construction
/// of this kind of object.  See section on \ref METAFILTERS "meta-filters" for more
/// information.
template<class MTYPE, class FILTER>
class is_inf_within_sped_t : public unary_function<MTYPE, bool>
{
  public:
  
    /// Constructor.
    ///
    /// \param f The initial function we use for determining if a member
    ///          is informative.
    explicit is_inf_within_sped_t(const FILTER& f);
    
    /// Copy Constructor
    ///
    is_inf_within_sped_t(const is_inf_within_sped_t&);

    /// Destructor
    ///
    //lint -e{1510} unary_function has no destructor, but that's no big thing
    ~is_inf_within_sped_t();
    
    /// Copy Operator
    ///
    is_inf_within_sped_t& operator=(const is_inf_within_sped_t&);
    
    /// Our test function
    ///
    /// \param member The member to test
    bool operator()(const MTYPE& member) const;
    
  private:

    /// Default constructor is private to prevent default creation.
    ///
    //lint -e{1704}
    is_inf_within_sped_t();

    enum InfChildStateEnum
    {
      NONE = 0,
      ONE,
      TWO_OR_MORE
    };
  
    /// Classifies the number of informative children into one of three states.
    /// Informativity of children is determined by the current state of
    /// \c my_informative_states.
    /// 
    /// \param mem         The member whose informative chilren we wish to count
    /// \returns A pair, containing the informative children status and
    ///          one of the informative children (generally the last found).
    pair<InfChildStateEnum, const MPED::member_base*>
        classify_inf_children(const MPED::member_base&  mem) const;
    
    /// Calculates the informativity of every member of my_current_subpedigree
    /// and stores it within my_informative_states
    void calc_subped_inf_vector() const;
    
    FILTER my_filter;

    mutable const typename MTYPE::subpedigree_type* my_current_subpedigree;
    
    mutable vector<bool> my_informative_states;
};

/// \brief Function which creates an is_inf_within_sped_t object
///
/// This function allows a much easier creation of an is_inf_within_sped_t object
/// without resorting to explicit templatization.  This is similar to what the
/// std::make_pair() function does for pairs in the STL, as well as various others.
/// See section on \ref METAFILTERS "meta-filters" for more
/// information.
template<class FILTER>
is_inf_within_sped_t<typename FILTER::argument_type, FILTER> is_inf_within_sped(const FILTER& f);

/// \brief Functor which, given a member, returns true if they're informative under the constraints (see detailed description).
///
/// This filter is a meta-filter.  It takes another filter as an argument
/// at construction, and uses that to determine informativity of members within
/// a family.
///
/// An individual returns true for this filter if:
///
/// -# They belong to the family specified at construction time,
/// -# There are at least two children informative by the filter function in this family,
/// -# They are either an informative child or a parent.
///
/// ie, if there aren't two informative children, no one in the family is informative.
/// If there are at least two, then all informative children and the parents are
/// retained, and uninformative children are removed.
///
/// Note that parental informativity \b doesn't matter in this case.  Only the
/// children's state matters.  Primarily, this functor is provided for GENIBD's use,
/// and parent-child interactions are not relevant there, only sibling pairs.  If
/// parent-child connections must be preserved, consider using the is_inf_within_subped
/// functor.
///
/// Like other meta-filters, this has the _t appended to the end to allow for the
/// is_inf_within_sped() function to use that name.  This simplifies construction
/// of this kind of object.  See section on \ref METAFILTERS "meta-filters" for more
/// information.
template<class MTYPE, class FILTER>
class is_family_sib_pair_inf_t : public unary_function<MTYPE, bool>
{
  public:
    /// \name Basic Maintenance Operations
    //@{
  
    /// Constructor.
    ///
    /// \param fam The family we're interested in.
    /// \param f   The initial function we use for determining if a member
    ///            is informative.
    explicit is_family_sib_pair_inf_t(const typename MTYPE::family_type& fam, const FILTER& f);
    
    /// Copy Constructor
    ///
    is_family_sib_pair_inf_t(const is_family_sib_pair_inf_t&);

    /// Destructor
    ///
    //lint -e{1510} unary_function has no destructor, but that's no big thing
    ~is_family_sib_pair_inf_t();
    
    /// Copy Operator
    ///
    is_family_sib_pair_inf_t& operator=(const is_family_sib_pair_inf_t&);
    
    //@}
    
    /// \name Family Get/Set operations
    //@{
    
    /// Set the current family.
    ///
    /// \param fam The family to set to
    void set_family(const typename MTYPE::family_type& fam);

    /// Get the current family.
    ///
    /// \returns The current family.
    const typename MTYPE::family_type& get_family() const;

    //@}

    /// Our test function
    ///
    /// \param member The member to test
    bool operator()(const MTYPE& member) const;
    
  private:
  
    /// Returns \c true if the member is a parent of the current family,
    /// \c false otherwise.
    bool is_member_a_parent(const MTYPE& mem) const;

    /// Returns \c true if the member is a parent of the current family,
    /// \c false otherwise.
    bool is_member_a_child(const MTYPE& mem) const;

    /// Counts the number of informative sibs.  Returns \c true if
    /// the number of informative sibs is >= 2, \c false otherwise.
    bool does_fam_have_enough_inf_sibs() const;
    
    /// Default constructor is private to prevent default creation.
    ///
    //lint -e{1704}
    is_family_sib_pair_inf_t();

    FILTER my_filter;

    const typename MTYPE::family_type* my_current_family;
    
    bool my_family_has_informative_spairs;
};

/// \brief Function which creates an is_family_sib_pair_inf_t object
///
/// This function allows a much easier creation of an is_family_sib_pair_inf_t object
/// without resorting to explicit templatization.  This is similar to what the
/// std::make_pair() function does for pairs in the STL, as well as various others.
/// See section on \ref METAFILTERS "meta-filters" for more
/// information.
template<class FTYPE, class FILTER>
is_family_sib_pair_inf_t<typename FTYPE::member_type, FILTER>
    is_family_sib_pair_inf(const FTYPE& fam, const FILTER& f);

//@}
  
/// \brief Creates a new pedigree with a single nuclear family with the two individuals as sibs and dummy parents.
///    
/// Primarily for RELTEST, filter_to_sib_pair takes a pair of individuals and
/// adds a nuclear family in the target FilteredMultipedigree with the pair
/// as siblings.  Parents are assigned names "~dummy1" and "~dummy2" within this
/// family, and are given no source member (they are considered to have missing data).
/// The pedigree name is assigned based upon ped_name as follows:
///
/// When the ped_name field is empty, the function creates a custom name based upon
/// the names of the two members and their pedigrees. The structure of this is
/// "P1:I1 x P2:I2", where P1 and P2 are the pedigrees of the first and second member,
/// respectively, and I1 and I2 are their names.  This should mostly assure uniqueness
/// of pedigree names to allow us to add many of these pairs to a single
/// FilteredMultipedigree.
///
/// This is a simple case of building a new structure during filtering, which is
/// necessary given RELTEST's analysis (which uses an assumed sibling relationship
/// to test for allele sharing percentages).  It could be moved to that library, but
/// is left in fped as an example of building new structures independent of the
/// RPED::RefMultiPedigree's structure.
///
/// \param fped     The FilteredMultipedigree target
/// \param ind1     The first sibling
/// \param ind2     The second sibling
/// \param ped_name Pedigree name to assign in the fped.  This is provided so
///                 that multiple new sibling groups can be created with
///                 unique pedigree identifiers.  
///
/// \ingroup FPEDFilterObjects
template <class MTYPE>
void filter_to_sib_pair(FilteredMultipedigree& fped,
                        const MTYPE& ind1,
                        const MTYPE& ind2,
                        string ped_name = "");



} // End FPED namespace
} // End SAGE namespace

#include "fped/fped_func.ipp"

#endif
