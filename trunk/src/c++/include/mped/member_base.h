#ifndef MEMBER_BASE_H
#define MEMBER_BASE_H

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
  * \brief This class provides the basic representation of a person in a multi-pedigree object.
  *
  * This is a very simple class, providing only accessor
  * functions as its public members.  This is so the derived
  * class template \code 'member<GI>' \endcode can be used by clients of the
  * multi-pedigree type.
  *
  * Its private methods are used by its friends 'family_base' 
  * and 'pedigree_base' to build the graph of the pedigree.  It 
  * performs no building in and of itself -- all the expertise 
  * is with 'family_base' and 'pedigree_base'.
  */
class member_base
{
  public:
    friend  class family_base;
    friend  class subpedigree_base;
    friend  class pedigree_base;
    friend  class PedigreeBuilder;
    friend  class multipedigree_base;
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

  ~member_base() {}

  /// @name Attributes
  //@{

    ///
    /// Returns the sex of the individual, reduced to basic values
    /// (SEX_MALE, SEX_FEMALE or SEX_MISSING).
    SexCode get_effective_sex() const;

    ///
    /// Returns the sex of the individual, completely detailed (including
    /// the inference bit)
    SexCode get_detailed_sex() const;

    ///
    /// Assigns this member's sex to \c s.
    /// \param s The sex of this member
    void set_sex(SexCode s);

    ///
    /// Returns \c true if male, \c false if female or missing
    bool is_male() const;

    ///
    /// Returns \c true if female, \c false if male or missing
    bool is_female() const;

    ///
    /// Returns \c true if sex is missing or arbitrary, \c false otherwise
    bool is_sex_unknown() const;

    ///
    /// Returns the index number of the individual within the multipedigree.  
    uint mpindex() const;

    ///
    /// Returns the index number of the individual within the pedigree.
    uint index() const;

    ///
    /// Returns the index number of the individual within the subpedigree.
    uint subindex() const;

    ///
    /// Returns the number of mates of this individual.
    uint mate_count() const;

    ///
    /// Returns the number of offspring of this individual.
    uint offspring_count() const;

    ///
    /// Returns the number of siblings of this individual.
    uint sibling_count() const;

    ///
    /// Returns the name of this individual (generally, the number of the individual in the original source data)
    const string & name() const;

  //@}

  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to this individual's multipedigree.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to this individual's multipedigree.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to this individual's pedigree.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to this individual's pedigree.
    pedigree_const_pointer pedigree() const;

    ///
    /// Returns a non-const pointer to this individual's subpedigree.
    subpedigree_pointer subpedigree();

    ///
    /// Returns a const pointer to this individual's subpedigree.
    subpedigree_const_pointer subpedigree() const;

    ///
    /// Returns a non-const pointer to this individual's family.
    family_pointer family();

    ///
    /// Returns a const pointer to this individual's family.
    family_const_pointer family() const;

    ///
    /// Returns a non-const pointer to this individual's first parent.
    member_pointer parent1();

    ///
    /// Returns a const pointer to this individual's first parent.
    member_const_pointer parent1() const;

    ///
    /// Returns a non-const pointer to this individual's second parent.
    member_pointer parent2();

    ///
    /// Returns a const pointer to this individual's second parent.
    member_const_pointer parent2() const;

    member_pointer       get_mother();
    member_const_pointer get_mother() const;
    member_pointer       get_father();
    member_const_pointer get_father() const;

  //@}

  /// @name Mate iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this individual's mate list.
    mate_iterator mate_begin();

    ///
    /// Returns a const begin iterator for this individual's mate list.
    mate_const_iterator mate_begin() const;

    ///
    /// Returns a non-const end iterator for this individual's mate list.
    mate_iterator mate_end();

    ///
    /// Returns a const end iterator for this individual's mate list.
    mate_const_iterator mate_end() const;

  //@}

  /// @name Offspring iterators
  //@{

    ///
    /// Given a mate, returns a non-const begin iterator for this individual's offspring list.
    /// \param i The mate of this individual
    offspring_iterator offspring_begin(const mate_const_iterator& i);

    ///
    /// Given a mate, returns a const begin iterator for this individual's offspring list.
    /// \param i The mate of this individual
    offspring_const_iterator offspring_begin(const mate_const_iterator& i) const;

    ///
    /// Given a mate, returns a non-const end iterator for this individual's offspring list.
    /// \param i The mate of this individual
    offspring_iterator offspring_begin(const member_const_iterator& i);

    ///
    /// Given a mate, returns a const begin iterator for this individual's offspring list.
    /// \param i The mate of this individual
    offspring_const_iterator offspring_begin(const member_const_iterator& i) const;

    ///
    /// Given a mate, returns a non-const begin iterator for this individual's offspring list.
    /// \param m The mate of this individual
    offspring_iterator offspring_begin(const member_base& m);

    ///
    /// Given a mate, returns a non-const begin iterator for this individual's offspring list.
    /// \param m The mate of this individual
    offspring_const_iterator offspring_begin(const member_base& m) const;

    ///
    /// Given a mate, returns a non-const begin iterator for this individual's offspring list.
    /// \param m The mate of this individual
    offspring_iterator offspring_begin(const string& m);

    ///
    /// Given a mate, returns a const begin iterator for this individual's offspring list.
    /// \param m The mate of this individual
    offspring_const_iterator offspring_begin(const string& m) const;

    ///
    /// Returns a non-const end iterator for this individual's offspring list.
    offspring_iterator offspring_end();

    ///
    /// Returns a const end iterator for this individual's offspring list.
    offspring_const_iterator offspring_end() const;

  //@}

  /// @name Parent iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this individual's parent list.
    parent_iterator parent_begin();

    ///
    /// Returns a const begin iterator for this individual's parent list.
    parent_const_iterator parent_begin() const;

    ///
    /// Returns a non-const end iterator for this individual's parent list.
    parent_iterator parent_end();

    ///
    /// Returns a const end iterator for this individual's parent list.
    parent_const_iterator parent_end() const;

  //@}

  /// @name Progeny iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this individual's progeny list.
    progeny_iterator progeny_begin();

    ///
    /// Returns a const begin iterator for this individual's progeny list.
    progeny_const_iterator progeny_begin() const;

    ///
    /// Returns a non-const end iterator for this individual's progeny list.
    progeny_iterator progeny_end();

    ///
    /// Returns a const end iterator for this individual's progeny list.
    progeny_const_iterator progeny_end() const;

  //@}

  /// @name Sibling iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this individual's sibling list.
    sibling_iterator sibling_begin();

    ///
    /// Returns a const begin iterator for this individual's sibling list.
    sibling_const_iterator sibling_begin() const;

    ///
    /// Returns a non-const end iterator for this individual's sibling list.
    sibling_iterator sibling_end();

    ///
    /// Returns a const end iterator for this individual's sibling list.
    sibling_const_iterator sibling_end() const;

  //@}

  /// \name Member Relationship Functions
  ///
  /// These functions provide answers to common questions about the kind of relationships
  /// a member has with other members
  /// 
  //@{

    /// Returns \c true if the individual is a founder, \c false otherwise.
    /// Unconnected individuals always return \c true
    bool is_founder() const;

    /// Returns \c true if the individual is a nonfounder, \c false otherwise.
    /// Unconnected individuals always return \c false
    bool is_nonfounder() const;

    /// Returns \c true if the individual belongs to a valid subpedigree, 
    /// \c false otherwise.
    bool is_connected() const;

    /// Returns \c true if the individual does not belong to a valid subpedigree, 
    /// \c false otherwise.
    bool is_unconnected() const;

    /// Returns \c true if the individual is a parent, 
    /// \c false otherwise.  Unconnects always return \c false
    bool is_parent() const;

    /// Returns \c true if the individual is not a parent, 
    /// \c false otherwise.  Unconnects always return \c true
    bool is_non_parent() const;
  //@}

  protected:

    member_base(const string& name, SexCode s=SEX_MISSING);
    member_base(const string& name, const no_info&);

  private:

    member_base();
    member_base(const member_base& P);
    void    operator=(const member_base& P);

    member_id   siblings() const;
    uint        add_child();
    uint        add_mate();
    void        set_mpindex(uint i);
    void        set_index(uint i);
    void        set_origin(family_id F);
    void        set_sibling(member_id P);
    void        set_subped(subped_id SP);
    void        set_subindex(uint i);

    uint        my_mpindex;
    uint        my_index;
    uint        my_subindex;
    string      my_name;
    subped_id   my_subped;
    family_id   my_origin;
    member_id   my_siblings;
    SexCode     my_gender;
    uint        my_offspring_count;
    uint        my_mate_count;

};

}
}

#endif
