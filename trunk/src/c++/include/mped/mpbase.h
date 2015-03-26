#ifndef _MPBASE_HPP
#define _MPBASE_HPP

//============================================================================
//  File:       mpbase.h
//
//  Purpose:    This header file defined the new SAGE multi-pedigree data
//              type.
//
//  Author:     Bob Steagall
//
//  History:    Version 0.90
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mped/spbase.h"

namespace SAGE {
/// MPED Namespace
namespace MPED {

/** \ingroup BaseStorageClasses
  * \brief Stores pedigree structure only; trait data is stored via the templatized, derived version multipedigree
  *
  * \par Introduction
  * 
  * The multipedigree_base stores information about pedigree structure. It organizes pedigree data into
  * pedigree_base s (group of individuals sharing the same family id in the data source),
  * subpedigree_base s (group of individuals related by marriage or blood within a pedigree_base),
  * and family_base s (that is, nuclear families).
  *
  * For individuals, the multipedigree_base also tracks sibships, parents, and offspring.
  *
  * The multipedigree_base does \b not, however, store individual trait data. It is designed to be inherited
  * and, through the inherited class, implement a method for storing individual trait data. This is accomplished
  * through the multipedigree, as well as the filtered_multipedigree.
  *
  * \par Getting started
  *
  * There are two stages to using the multipedgree_base. First, you add entries for all individuals, including
  * information about sibships, marriages, parent/offspring lineages, etc. Second, you build the multipedigree_base,
  * which effectively finalizes the internal data storage.
  *
  * \par Adding entries
  *
  * You can add information about individuals using four functions: add_member(), add_marriage(), add_lineage(),
  * and add_sibship(). There are a few variants for each function, allowing you to use either pedigree name
  * (from the source data) or pedigree_iterator uniquely identifying the pedigree.
  *
  * To add a listing for a male individual '2' in pedigree '1', for instance, you would enter the following code:
  *
  * \code
  * SAGE::MPED::multipedigree_base my_multipedigree_base();
  * my_multipedigree_base.add_member("1", "2", SAGE::MPED::male);
  * \endcode
  * 
  * \par Building the multipedigree_base
  *
  * Once you have added all the necessary entries for individuals, marriages, lineages, and sibships, you must
  * invoke build(). 
  *
  * build() takes care of finalizing the internal data storage and preparing it for use after data importation.
  *
  * Example:
  *
  * \code
  * SAGE::MPED::multipedigree_base my_multipedigree_base();
  * my_multipedigree_base.add_member("1", "2", SAGE::MPED::male);
  * my_multipedigree_base.build();
  * \endcode
  *
  * \par Indexing & iteration
  *
  * Functions are available for iterating across the stored pedigrees (both const and non-const).
  *
  * Also, there are functions available for looking up members, pedigrees, and querying the number of pedigrees.
  *
  * \par Additional notes
  *
  */
class multipedigree_base
{
  public:

  /** @name Typedefs
    * Typedef-ed local to this class.
    */
  //@{

    typedef SAGE::MPED::member_base                         member_type;
    typedef SAGE::MPED::family_base                         family_type;
    typedef SAGE::MPED::subpedigree_base                    subpedigree_type;
    typedef SAGE::MPED::pedigree_base                       pedigree_type;
    typedef SAGE::MPED::multipedigree_base                  multipedigree_type;

    typedef SAGE::MPED::member_base*                        member_pointer;
    typedef SAGE::MPED::family_base*                        family_pointer;
    typedef SAGE::MPED::subpedigree_base*                   subpedigree_pointer;
    typedef SAGE::MPED::pedigree_base*                      pedigree_pointer;
    typedef SAGE::MPED::multipedigree_base*                 multipedigree_pointer;
    typedef SAGE::MPED::family_base_iterator                family_iterator;
    typedef SAGE::MPED::mate_base_iterator                  mate_iterator;
    typedef SAGE::MPED::member_base_iterator                member_iterator;
    typedef SAGE::MPED::offspring_base_iterator             offspring_iterator;
    typedef SAGE::MPED::parent_base_iterator                parent_iterator;
    typedef SAGE::MPED::pedigree_base_iterator              pedigree_iterator;
    typedef SAGE::MPED::progeny_base_iterator               progeny_iterator;
    typedef SAGE::MPED::sibling_base_iterator               sibling_iterator;
    typedef SAGE::MPED::subpedigree_base_iterator           subpedigree_iterator;

    typedef const SAGE::MPED::member_base*                  member_const_pointer;
    typedef const SAGE::MPED::family_base*                  family_const_pointer;
    typedef const SAGE::MPED::subpedigree_base*             subpedigree_const_pointer;
    typedef const SAGE::MPED::pedigree_base*                pedigree_const_pointer;
    typedef const SAGE::MPED::multipedigree_base*           multipedigree_const_pointer;
    typedef SAGE::MPED::family_base_const_iterator          family_const_iterator;
    typedef SAGE::MPED::mate_base_const_iterator            mate_const_iterator;
    typedef SAGE::MPED::member_base_const_iterator          member_const_iterator;
    typedef SAGE::MPED::offspring_base_const_iterator       offspring_const_iterator;
    typedef SAGE::MPED::parent_base_const_iterator          parent_const_iterator;
    typedef SAGE::MPED::pedigree_base_const_iterator        pedigree_const_iterator;
    typedef SAGE::MPED::progeny_base_const_iterator         progeny_const_iterator;
    typedef SAGE::MPED::sibling_base_const_iterator         sibling_const_iterator;
    typedef SAGE::MPED::subpedigree_base_const_iterator     subpedigree_const_iterator;

  //@}

  /// @name Constructor & destructor
  //@{

    ///
    /// Constructor.
    multipedigree_base();

    ///
    /// Destructor.
    virtual ~multipedigree_base();

  //@}

  /// @name Adding members & relationships
  //@{

    ///
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param ped The name of the pedigree in which the member is indicated
    /// \param name The name of the member to be added
    /// \param s The sex of the individual (default is SEX_MISSING)
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const string& ped, const string& name, SexCode s=SEX_MISSING);

    ///
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the member is indicated
    /// \param name The name of the member to be added
    /// \param s The sex of the individual (default is SEX_MISSING)
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const pedigree_iterator& ped, const string& name, SexCode s=SEX_MISSING);

    ///
    /// Amends the internal pedigree storage to reflect a spousal relationship.
    /// \param ped The name of the pedigree in which the marriage is indicated
    /// \param spouse1 The name of the first spouse in the relationship
    /// \param spouse2 The name of the second spouse in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_marriage(const string& ped, const string& spouse1, const string& spouse2);

    ///
    /// Amends the internal pedigree storage to reflect a spousal relationship.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the marriage is indicated
    /// \param spouse1 The name of the first spouse in the relationship
    /// \param spouse2 The name of the second spouse in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_marriage(const pedigree_iterator& ped, const string& spouse1, const string& spouse2);

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param ped The name of the pedigree in which the lineage relationship is indicated
    /// \param child The name of the child in the relationship
    /// \param parent The name of the parent in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_lineage(const string& ped, const string& child, const string& parent);

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the lineage relationship is indicated
    /// \param child The name of the child in the relationship
    /// \param parent The name of the parent in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_lineage(const pedigree_iterator& ped, const string& child, const string& parent);

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param ped The name of the pedigree in which the lineage relationship is indicated
    /// \param child The name of the child in the relationship
    /// \param parent1 The name of the first parent in the relationship
    /// \param parent2 The name of the second parent in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_lineage(const string& ped, const string& child, const string& parent1, const string& parent2);

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the lineage relationship is indicated
    /// \param child The name of the child in the relationship
    /// \param parent1 The name of the first parent in the relationship
    /// \param parent2 The name of the second parent in the relationship
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_lineage(const pedigree_iterator& ped, const string& child, const string& parent1, const string& parent2);

    ///
    /// Amends the internal pedigree storage to reflect a sibling relationship.
    /// \param ped The name of the pedigree in which the sibship is indicated
    /// \param sib1 The name of the one of the sibs
    /// \param sib2 The name of another sib
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_sibship(const string& ped, const string& sib1, const string& sib2);

    ///
    /// Amends the internal pedigree storage to reflect a sibling relationship.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the sibship is indicated
    /// \param sib1 The name of the one of the sibs
    /// \param sib2 The name of another sib
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_sibship(const pedigree_iterator& ped, const string& sib1, const string& sib2);

  //@}

  /// @name Building functions
  //@{

    ///
    /// Having added entries for all members as well as structural relationships, this function finalizes
    /// the internal storage for all pedigrees.
    void build();

    ///
    /// Having invoked build(), you can then use freeze() to disallow any additional changes to the
    /// data. Any subsequent attempt to alter the data will cause an error.
    void freeze();

    ///
    /// Having invoked build(), you can then use freeze() to disallow any additional changes to the
    /// data in the indicated pedigree. Any subsequent attempt to alter the data will cause an error.
    /// \param ped The name of the pedigree
    /// \return An iterator for the indicated pedigree
    pedigree_iterator   freeze(const string& ped);

    ///
    /// Having invoked build(), you can then use freeze() to disallow any additional changes to the
    /// data in the indicated pedigree. Any subsequent attempt to alter the data will cause an error.
    /// \param ped A const_iterator for the indicated pedigree
    /// \return An iterator for the indicated pedigree
    pedigree_iterator   freeze(const pedigree_iterator& ped);

  //@}

  /// @name Descriptive statistics
  //@{

    ///
    /// Returns the number of pedigrees.
    uint pedigree_count() const;

    ///
    /// Returns the number of members.
    uint member_count() const;

  //@}

  /// @name Iterating across pedigrees
  //@{

    ///
    /// Returns a non-const begin iterator for the pedigrees.
    pedigree_iterator pedigree_begin();

    ///
    /// Returns a non-const end iterator for the pedigrees.
    pedigree_iterator pedigree_end();

    ///
    /// Returns a non-const iterator for the last pedigree.
    pedigree_iterator pedigree_last();

    ///
    /// Returns a const begin iterator for the pedigrees.
    pedigree_const_iterator pedigree_begin() const;

    ///
    /// Returns a const end iterator for the pedigrees.
    pedigree_const_iterator pedigree_end() const;

    ///
    /// Returns a const iterator for the last pedigree.
    pedigree_const_iterator pedigree_last() const;

  //@}

  /// @name Pedigree/member lookup functions
  //@{

    ///
    /// Given a pedigree id number (coming from the multipedigree_base), returns the corresponding
    /// non-const pedigree_base.
    /// \param i The pedigree id number
    /// \return The corresponding pedigree_base
    pedigree_base & pedigree_index(uint i);

    ///
    /// Given a pedigree id number (coming from the multipedigree_base), returns the corresponding
    /// const pedigree_base.
    /// \param i The pedigree id number
    /// \return The corresponding pedigree_base
    const pedigree_base & pedigree_index(uint i) const;

    ///
    /// Swaps the internally stored pedigree id numbers.
    /// \param i The id number of the first pedigree
    /// \param j The id number of the second pedigree
    void pedigree_index_swap(uint i, uint j);

    ///
    /// Given a pedigree name, returns the corresponding non-const pedigree pointer.
    /// \param P The name of the pedigree (as indicated in the pedigree data file)
    pedigree_pointer pedigree_find(const string& P);

    ///
    /// Given a pedigree name, returns the corresponding const pedigree pointer.
    /// \param P The name of the pedigree (as indicated in the pedigree data file)
    pedigree_const_pointer pedigree_find(const string& P) const;

    ///
    /// Given a member id number (coming from the multipedigree_base), returns the corresponding
    /// non-const member_base.
    /// \param i The member id number
    /// \return The corresponding member_base
    member_base & member_index(uint i);

    ///
    /// Given a member id number (coming from the multipedigree_base), returns the corresponding
    /// const member_base.
    /// \param i The member id number
    /// \return The corresponding member_base
    const member_base & member_index(uint i) const;

    ///
    /// Swaps the internally stored member id numbers.
    /// \param i The id number of the first member
    /// \param j The id number of the second member
    void member_index_swap(uint i, uint j);

    ///
    /// Given a pedigree name and a member name, ("3" and "2", for instance), returns the
    /// corresponding non-const pointer to the member.
    /// \param ped The name of the pedigree to which the member belong (as indicated in the pedigree data file)
    /// \param name The name of the member (as indicated in the pedigree data file)
    member_pointer member_find(const string & ped, const string & name);

    ///
    /// Given a pedigree name and a member name, ("3" and "2", for instance), returns the
    /// corresponding const pointer to the member.
    /// \param ped The name of the pedigree to which the member belong (as indicated in the pedigree data file)
    /// \param name The name of the member (as indicated in the pedigree data file)
    member_const_pointer member_find(const string & ped, const string & name) const;

  //@}

  protected:

  /// @name Additional helper functions
  //@{

    ///
    /// Returns the name of the most recently added pedigree.
    const string & last_name() const;

    ///
    /// Although the multipedigree_base doesn't technically store trait data information for
    /// individuals, it *does* in fact store a void pointer for each individual. Then, on the derived
    /// level, the void pointer is cast as the correct info type.
    ///
    /// Consequently, you can use the add_member() function within the multipedigree_base to add
    /// an individual with an entry for that individual's info object.
    /// \param ped A pedigree_iterator for the pedigree to which the individual belongs
    /// \param name The name of the individual
    /// \param s The sex of the individual
    /// \param info The const pointer to the individual's info
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const pedigree_iterator& ped, const string& name, SexCode s, const void* info);

    ///
    /// Although the multipedigree_base doesn't technically store trait data information for
    /// individuals, it *does* in fact store a void pointer for each individual. Then, on the derived
    /// level, the void pointer is cast as the correct info type.
    ///
    /// Consequently, you can use the add_member() function within the multipedigree_base to add
    /// an individual with an entry for that individual's info object.
    /// \param ped The name of the pedigree to which the individual belongs
    /// \param name The name of the individual
    /// \param s The sex of the individual
    /// \param info The const pointer to the individual's info
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const string& ped, const string& name, SexCode s, const void* info);

  //@}

  /// @name Pedigree allocation / deallocation
  //@{

    ///
    /// Allocates an entry for the indicated pedigree.
    ///
    /// Please note that this function is declared virtual in case the derived type needs to
    /// override it in any way.
    /// \param name The name of the pedigree
    /// \return The id number of the pedigree
    virtual pedigree_id allocate_pedigree(const string & name);

    ///
    /// De-allocates an entry for the indicated pedigree.
    ///
    /// Please note that this function is declared virtual in case the derived type needs to
    /// override it in any way.
    /// \param P The pedigree_id of the pedigree
    virtual void deallocate_pedigree(pedigree_id P);

  //@}

    ///
    /// Clears all entries from all pedigrees.
    void clear();

    ///
    /// Having added entries for all members as well as structural relationships, this function finalizes
    /// the internal storage for a specific pedigree.
    /// \param ped The name of the pedigree to build
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator build(const string& ped);

    ///
    /// Having added entries for all members as well as structural relationships, this function finalizes
    /// the internal storage for a specific pedigree.
    /// \param ped A const pedigree_iterator corresponding to the pedigree to build
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator build(const pedigree_iterator& ped);


  private:

    pedigree_iterator   flush(const string& ped);
    pedigree_iterator   flush(const pedigree_iterator& ped);
    void        flush();

    /// \internal
    /// Clears all entries from the indicated pedigree.
    /// \param ped The name of the pedigree to clear
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator   clear(const string& ped);

    /// \internal
    /// Clears all entries from the indicated pedigree.
    /// \param ped A const pedigree_iterator corresponding to the pedigree to clear
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator   clear(const pedigree_iterator& ped);



    typedef std::pair<pedigree_map::iterator,bool>  insert_pair;
    typedef pedigree_map::value_type                value_pair;

    mutable pedigree_map        my_pedmap;
    mutable pidx_vector         my_ped_index;
    mutable midx_vector         my_mem_mpindex;
    mutable string              my_last_name;
    mutable pedigree_iterator   my_last_iter;

    multipedigree_base(const multipedigree_base&);
    void    operator=(const multipedigree_base&);

    pedigree_iterator   update_cache(const string& ped);
    pedigree_iterator   update_cache(const pedigree_iterator& ped);
    void                update_indices(pedigree_id P);

    void                assign_member_indices();
};


} // End namespace MPED
} // End namespace SAGE

#include "mped/mpbase.ipp"

#endif  //- _MPBASE_HPP
