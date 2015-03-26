#ifndef FAMILY_BASE_H
#define FAMILY_BASE_H

//============================================================================
//  File:       spbase.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Notes:      This header defines the containers used to manipulate
//              single pedigree data and its components.
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mped/spbaseiter.h"

namespace SAGE {
namespace MPED {

/** \ingroup BaseStorageClasses
  * \brief This class provides the basic representation of a nuclear family in a multi-pedigree object
  *
  * This is also fairly simple class, providing only accessor 
  * functions as its public members.
  *
  * However, this class has the additional responsibility of 
  * building the pedigree graph when called upon to do so by its
  * friend class 'pedigree_base'. 
  *
  * It should also be noted that the names of the parents are
  * always kept in lexicographical order.  That is, 'name1()' 
  * will always return a name string that is lexicographically
  * less than 'name2()'.  This simplifies the task of looking up
  * families by the names of the parents.
  */
class family_base
{
  public:
    friend  class member_base;
    friend  class subpedigree_base;
    friend  class pedigree_base;
    friend  class PedigreeBuilder;
    friend  class cursor_proxy;

  /** \name Typdefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base                         member_type;
    typedef family_base                         family_type;
    typedef subpedigree_base                    subpedigree_type;
    typedef pedigree_base                       pedigree_type;
    typedef multipedigree_base                  multipedigree_type;

    typedef member_base*                        member_pointer;
    typedef family_base*                        family_pointer;
    typedef subpedigree_base*                   subpedigree_pointer;
    typedef pedigree_base*                      pedigree_pointer;
    typedef multipedigree_base*                 multipedigree_pointer;
    typedef family_base_iterator                family_iterator;
    typedef mate_base_iterator                  mate_iterator;
    typedef member_base_iterator                member_iterator;
    typedef offspring_base_iterator             offspring_iterator;
    typedef parent_base_iterator                parent_iterator;
    typedef pedigree_base_iterator              pedigree_iterator;
    typedef progeny_base_iterator               progeny_iterator;
    typedef sibling_base_iterator               sibling_iterator;
    typedef subpedigree_base_iterator           subpedigree_iterator;

    typedef const member_base*                  member_const_pointer;
    typedef const family_base*                  family_const_pointer;
    typedef const subpedigree_base*             subpedigree_const_pointer;
    typedef const pedigree_base*                pedigree_const_pointer;
    typedef const multipedigree_base*           multipedigree_const_pointer;
    typedef family_base_const_iterator          family_const_iterator;
    typedef mate_base_const_iterator            mate_const_iterator;
    typedef member_base_const_iterator          member_const_iterator;
    typedef offspring_base_const_iterator       offspring_const_iterator;
    typedef parent_base_const_iterator          parent_const_iterator;
    typedef pedigree_base_const_iterator        pedigree_const_iterator;
    typedef progeny_base_const_iterator         progeny_const_iterator;
    typedef sibling_base_const_iterator         sibling_const_iterator;
    typedef subpedigree_base_const_iterator     subpedigree_const_iterator;

  //@}

  public:

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor
    family_base();

    ///
    /// Destructor
    ~family_base() {}

  //@}

  /// @name Attributes
  //@{

    ///
    /// Returns the index number of this family within the pedigree.
    uint index() const;

    ///
    /// Returns the index number of this family within the subpedigree.
    uint subindex() const;

    ///
    /// Returns the number of offspring in this family.
    uint offspring_count() const;

    ///
    /// Returns a string formed by concatenating the names of both parents in the family, with a colon separator character.
    ///
    /// Equivalent to \code name1() + ":" + name2() \endcode
    string name() const;

    ///
    /// Returns the name of the first parent in this family.
    const string & name1() const;

    ///
    /// Returns the name of the second parent in this family.
    const string & name2() const;

  //@}

  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to this family's multipedigree.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to this family's multipedigree.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to this family's pedigree.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to this family's pedigree.
    pedigree_const_pointer pedigree() const;

    ///
    /// Returns a non-const pointer to this family's subpedigree.
    subpedigree_pointer subpedigree();

    ///
    /// Returns a const pointer to this family's subpedigree.
    subpedigree_const_pointer subpedigree() const;

    ///
    /// Returns a non-const pointer to the first parent in this family.
    member_pointer parent1();

    ///
    /// Returns a const pointer to the first parent in this family.
    member_const_pointer parent1() const;

    ///
    /// Returns a non-const pointer to the second parent in this family.
    member_pointer parent2();

    ///
    /// Returns a const pointer to the second parent in this family.
    member_const_pointer parent2() const;

    member_pointer       get_mother();
    member_const_pointer get_mother() const;
    member_pointer       get_father();
    member_const_pointer get_father() const;
    
  //@}

  /// @name Pedigree structure iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this family's offspring list.
    offspring_iterator offspring_begin();

    ///
    /// Returns a const begin iterator for this family's offspring list.
    offspring_const_iterator offspring_begin() const;

    ///
    /// Returns a non-const end iterator for this family's offspring list.
    offspring_iterator offspring_end();

    ///
    /// Returns a const end iterator for this family's offspring list.
    offspring_const_iterator offspring_end() const;

    ///
    /// Returns a non-const begin iterator for this family's parent list.
    parent_iterator parent_begin();

    ///
    /// Returns a const begin iterator for this family's parent list.
    parent_const_iterator parent_begin() const;

    ///
    /// Returns a non-const end iterator for this family's parent list.
    parent_iterator parent_end();

    ///
    /// Returns a const end iterator for this family's parent list.
    parent_const_iterator parent_end() const;

  //@}

  private:

    family_base(const family_base&);
    family_base(member_id P1, member_id P2, const no_info&);
    family_base(subped_id SP, member_id P1, member_id P2, member_id C);
    void    operator=(const family_base&);

    member_id   offspring() const;
    bool        in_offspring(member_id P) const;
    void        add_offspring(member_id P);
    void        reset_links(subped_id SP, member_id P1, member_id P2, member_id C);
    void        set_index(uint i);
    void        set_subped(subped_id SP);
    void        set_subindex(uint i);

    uint        my_index;
    uint        my_subindex;
    subped_id   my_subped;
    member_id   my_parent1;
    member_id   my_parent2;
    member_id   my_mother;
    member_id   my_father;
    member_id   my_offspring;
    uint        my_offspring_count;
};

}
}

#endif

