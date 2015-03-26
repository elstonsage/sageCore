#ifndef _MP_HPP
#define _MP_HPP

//============================================================================
//  File:       mp.h
//
//  Purpose:    This header file defined the new SAGE multi-pedigree data
//              type.
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mped/sp.h"

namespace SAGE {
namespace MPED {

/** \ingroup DerivedStorageClasses
  * \brief This class is a templatized, derived form of the multipedigree_base
  *
  * \par Template parameters
  *
  * - \c GI Info object associated with individuals (SAGE::MPED::member)
  * - \c FI Info object associated with families (SAGE::MPED::family)
  * - \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * - \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * - \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  *
  * \par Introduction
  *
  * First of all, you should review the detailed description of the multipedigree_base. In brief,
  * the multipedigree_base takes care of organized pedigree structure, but that's it! Any data you might
  * want associated with the structural components (multipedigree, pedigree, subpedigree, family, or individual)
  * is done through templatization at the multipedigree level.
  *
  * \par Template parameters
  *
  * All of the derived storage classes in the mped library are templatized on all five template parameters,
  * even though each derived storage classes only makes use of the template parameter relevant to its storage type.
  *
  * The multipedigree makes use of the template parameter MI (MI stands for multipedigree info). If you are planning
  * on building your own specialization of the multipedigree class, you can design your own info object to be 
  * associated with the multipedigree.
  *
  * \par Getting started
  *
  * Basically, the multipedigree acts identically to the multipedigree_base, except that it is templatized
  * on the various info parameters. First, you add structural information about the individuals; then, you invoke
  * build() to fully construct the internal data storage.
  *
  * In addition to the structural relationship functions standard in the multipedigree_base, the multipedigree
  * also has two more versions of add_member() available. Both versions allow you to specify the info object
  * associated with the member that is being added. 
  *
  */
template <class GI, class FI=no_info, class SI=no_info, class PI=no_info, class MI=no_info>
class multipedigree : public multipedigree_base
{
  public:

  /** @name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef GI                                                        geninfo_type; // individual genetic information
    typedef FI                                                        faminfo_type; // family information
    typedef SI                                                        subinfo_type; // subpedigree information
    typedef PI                                                        pedinfo_type; // pedigree information
    typedef MI                                                        mpinfo_type;  // multipedigree information
    typedef multipedigree_base                                        base;
    typedef SAGE::MPED::member<GI,FI,SI,PI,MI>                        member_type;
    typedef SAGE::MPED::family<GI,FI,SI,PI,MI>                        family_type;
    typedef SAGE::MPED::subpedigree<GI,FI,SI,PI,MI>                   subpedigree_type;
    typedef SAGE::MPED::pedigree<GI,FI,SI,PI,MI>                      pedigree_type;
    typedef SAGE::MPED::multipedigree<GI,FI,SI,PI,MI>                 multipedigree_type;

    typedef member_type*                                              member_pointer;
    typedef family_type*                                              family_pointer;
    typedef subpedigree_type*                                         subpedigree_pointer;
    typedef pedigree_type*                                            pedigree_pointer;
    typedef multipedigree_type*                                       multipedigree_pointer;
    typedef SAGE::MPED::family_iterator<GI,FI,SI,PI,MI>               family_iterator;
    typedef SAGE::MPED::mate_iterator<GI,FI,SI,PI,MI>                 mate_iterator;
    typedef SAGE::MPED::member_iterator<GI,FI,SI,PI,MI>               member_iterator;
    typedef SAGE::MPED::offspring_iterator<GI,FI,SI,PI,MI>            offspring_iterator;
    typedef SAGE::MPED::parent_iterator<GI,FI,SI,PI,MI>               parent_iterator;
    typedef SAGE::MPED::pedigree_iterator<GI,FI,SI,PI,MI>             pedigree_iterator;
    typedef SAGE::MPED::progeny_iterator<GI,FI,SI,PI,MI>              progeny_iterator;
    typedef SAGE::MPED::sibling_iterator<GI,FI,SI,PI,MI>              sibling_iterator;
    typedef SAGE::MPED::subpedigree_iterator<GI,FI,SI,PI,MI>          subpedigree_iterator;

    typedef const member_type*                                        member_const_pointer;
    typedef const family_type*                                        family_const_pointer;
    typedef const subpedigree_type*                                   subpedigree_const_pointer;
    typedef const pedigree_type*                                      pedigree_const_pointer;
    typedef const multipedigree_type*                                 multipedigree_const_pointer;
    typedef SAGE::MPED::family_const_iterator<GI,FI,SI,PI,MI>         family_const_iterator;
    typedef SAGE::MPED::mate_const_iterator<GI,FI,SI,PI,MI>           mate_const_iterator;
    typedef SAGE::MPED::member_const_iterator<GI,FI,SI,PI,MI>         member_const_iterator;
    typedef SAGE::MPED::offspring_const_iterator<GI,FI,SI,PI,MI>      offspring_const_iterator;
    typedef SAGE::MPED::parent_const_iterator<GI,FI,SI,PI,MI>         parent_const_iterator;
    typedef SAGE::MPED::pedigree_const_iterator<GI,FI,SI,PI,MI>       pedigree_const_iterator;
    typedef SAGE::MPED::progeny_const_iterator<GI,FI,SI,PI,MI>        progeny_const_iterator;
    typedef SAGE::MPED::sibling_const_iterator<GI,FI,SI,PI,MI>        sibling_const_iterator;
    typedef SAGE::MPED::subpedigree_const_iterator<GI,FI,SI,PI,MI>    subpedigree_const_iterator;

  //@}

  /// @name Constructor & destructor
  //@{

    ///
    /// Constructor.
    multipedigree();

    ///
    /// Destructor.
    virtual ~multipedigree();

  //@}

  /// @name Multipedigree info object
  //@{ 

    ///
    /// Returns the multipedigree info object associated with the template-specific instantiation
    /// of this object.
    const mpinfo_type & info() const;

    ///
    /// Returns the multipedigree info object associated with the template-specific instantiation
    /// of this object.
    mpinfo_type & info();
    
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
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param ped The name of the pedigree in which the member is indicated
    /// \param name The name of the member to be added
    /// \param s The sex of the individual
    /// \param info The geninfo_type object associated with this member
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const string& ped, const string& name, SexCode s, const geninfo_type& info);

    ///
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param ped A const pedigree_iterator corresponding to the pedigree in which the member is indicated
    /// \param name The name of the member to be added
    /// \param s The sex of the individual
    /// \param info The geninfo_type object associated with this member
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator add_member(const pedigree_iterator& ped, const string& name, SexCode s, const geninfo_type& info);

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
    pedigree_iterator freeze(const string& ped);

    ///
    /// Having invoked build(), you can then use freeze() to disallow any additional changes to the
    /// data in the indicated pedigree. Any subsequent attempt to alter the data will cause an error.
    /// \param ped A const_iterator for the indicated pedigree
    /// \return An iterator for the indicated pedigree
    pedigree_iterator freeze(const pedigree_iterator& ped);

  //@}

  /// @name Descriptive statistics
  //@{

    ///
    /// Returns the number of pedigrees.
    uint pedigree_count() const;

    ///
    /// Returns the number of pedigrees.
    uint member_count() const;

  //@}

  /// @name Iterating across pedigrees
  //@{

    ///
    /// Returns a non-const begin iterator for this object's list of pedigrees.
    pedigree_iterator pedigree_begin();

    ///
    /// Returns a non-const end iterator for this object's list of pedigrees.
    pedigree_iterator pedigree_end();

    ///
    /// Returns a non-const last iterator for this object's list of pedigrees.
    pedigree_iterator pedigree_last();

    ///
    /// Returns a const begin iterator for this object's list of pedigrees.
    pedigree_const_iterator pedigree_begin() const;

    ///
    /// Returns a const end iterator for this object's list of pedigrees.
    pedigree_const_iterator pedigree_end() const;

    ///
    /// Returns a const last iterator for this object's list of pedigrees.
    pedigree_const_iterator pedigree_last() const;

  //@}

  /// @name Pedigree/member lookup functions
  //@{

    ///
    /// Returns the non-const pedigree object corresponding to the indicated index number.
    /// \param i The index number of the pedigree
    pedigree_type & pedigree_index(uint i);

    ///
    /// Returns the const pedigree object corresponding to the indicated index number.
    /// \param i The index number of the pedigree
    const pedigree_type & pedigree_index(uint i) const;

    ///
    /// Swaps the index numbers for the indicated pedigrees.
    /// \param i The first pedigree to swap
    /// \param j The second pedigree to swap
    void pedigree_index_swap(uint i, uint j);

    ///
    /// Returns a non-const pointer to the indicated pedigree.
    /// \param P The name of the pedigree
    pedigree_pointer pedigree_find(const string& P);

    ///
    /// Returns a const pointer to the indicated pedigree.
    /// \param P The name of the pedigree
    pedigree_const_pointer pedigree_find(const string& P) const;

    ///
    /// Returns the non-const member object corresponding to the indicated index number.
    /// \param i The index number of the member
    member_type & member_index(uint i);

    ///
    /// Returns the const member object corresponding to the indicated index number.
    /// \param i The index number of the member
    const member_type & member_index(uint i) const;

    ///
    /// Swaps the index numbers for the indicated members.
    /// \param i The first member to swap
    /// \param j The second member to swap
    void member_index_swap(uint i, uint j);

    ///
    /// Returns a non-const pointer to the indicated member.
    /// \param ped The name of the pedigree to which the member belongs
    /// \param name The name of the member
    member_pointer member_find(const string& ped, const string& name);

    ///
    /// Returns a const pointer to the indicated member.
    /// \param ped The name of the pedigree to which the member belongs
    /// \param name The name of the member
    member_const_pointer member_find(const string& ped, const string& name) const;

  //@}

    void        flush();
    pedigree_iterator   flush(const string& ped);
    pedigree_iterator   flush(const pedigree_iterator& ped);

  protected:

  /// @name Pedigree allocation / deallocation
  //@{

    /// 
    /// Allocates an entry for the indicated pedigree.
    ///
    /// \param name The name of the pedigree
    /// \return The id number of the pedigree
    virtual pedigree_id allocate_pedigree(const string& name);

    ///
    /// De-allocates an entry for the indicated pedigree.
    ///
    /// \param P The pedigree_id of the pedigree
    virtual void deallocate_pedigree(pedigree_id P);

  //@}

  private:
    multipedigree(const multipedigree&);
    void    operator=(const multipedigree&);

    /// 
    /// Clears all entries from all pedigrees.
    void clear();

    ///
    /// Clears all entries from the indicated pedigree.
    /// \param ped The name of the pedigree to clear
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator clear(const string& ped);

    ///
    /// Clears all entries from the indicated pedigree.
    /// \param ped A const pedigree_iterator corresponding to the pedigree to clear
    /// \return A pedigree_iterator for the indicated pedigree
    pedigree_iterator clear(const pedigree_iterator& ped);

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

    mpinfo_type    my_info;
};

} // End namespace MPED
} // End namespace SAGE

#include "mped/mp.ipp"

#endif  //- _MP_HPP
