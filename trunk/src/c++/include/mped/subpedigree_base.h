#ifndef SUBPEDIGREE_BASE_H
#define SUBPEDIGREE_BASE_H

#include "mped/spbaseiter.h"

namespace SAGE {
namespace MPED {

/** \ingroup BaseStorageClasses
  * \brief Provides a basic representation of a connected sub-pedigree within a pedigree in a multi-pedigree object
  *
  * This is also fairly simple class, providing only accessor 
  * functions as its members.  All of its members are private,
  * and are used primarily by its friends to establish ownership
  * or containment relationships.
  */
class subpedigree_base
{
  public:
    friend  class family_base;
    friend  class member_base;
    friend  class PedigreeBuilder;
    friend  class pedigree_base;
    friend  class cursor_proxy;

  /** \name Typedefs
    * Typdef-ed local to this class
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

  /// @name Destructor
  //@{

    ///
    /// Destructor
    ~subpedigree_base() {}

  //@}

  /// @name Attributes
  //@{

    ///
    /// Returns the name of the subpedigree in the form 'pedigree_name:subpedigree_index'
    const string & name() const;

    ///
    /// Returns the index number of this subpedigree within the pedigree.
    uint index() const;

    ///
    /// Returns the number of families in this subpedigree.
    uint family_count() const;

    ///
    /// Returns the number of individuals in this subpedigree.
    uint member_count() const;

  //@}

  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to this subpedigree's multipedigree.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to this subpedigree's multipedigree.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to this subpedigree's pedigree.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to this subpedigree's pedigree.
    pedigree_const_pointer pedigree() const;

  //@}

  /// @name Pedigree structure iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this pedigree's family list.
    family_iterator family_begin();

    ///
    /// Returns a const begin iterator for this pedigree's family list.
    family_const_iterator family_begin() const;

    ///
    /// Returns a non-const end iterator for this pedigree's family list.
    family_iterator family_end();

    ///
    /// Returns a const end iterator for this pedigree's family list.
    family_const_iterator family_end() const;

    ///
    /// Returns a non-const begin iterator for this pedigree's member list.
    member_iterator member_begin();

    ///
    /// Returns a const begin iterator for this pedigree's member list.
    member_const_iterator member_begin() const;

    ///
    /// Returns a non-const end iterator for this pedigree's member list.
    member_iterator member_end();

    ///
    /// Returns a const end iterator for this pedigree's member list.
    member_const_iterator member_end() const;

  //@}

  /// @name Indices
  //@{

    ///
    /// Given a family index number, returns the corresponding non-const family_base reference.
    /// \param i The index number of the family in question
    family_base & family_index(uint i);

    ///
    /// Given a family index number, returns the corresponding const family_base reference.
    /// \param i The index number of the family in question
    const family_base & family_index(uint i) const;

    ///
    /// Given a member index number, returns the corresponding non-const member_base reference.
    /// \param i The index number of the member in question
    member_base & member_index(uint i);

    ///
    /// Given a member index number, returns the corresponding const member_base reference.
    /// \param i The index number of the member in question
    const member_base & member_index(uint i) const;

    ///
    /// Swaps the index numbers of the indicated families.
    /// \param i The index number of the first family
    /// \param j The index number of the second family
    void family_index_swap(uint i, uint j);

    ///
    /// Swaps the index numbers of the indicated members.
    /// \param i The index number of the first member
    /// \param j The index number of the second member
    void member_index_swap(uint i, uint j);

  //@}

  protected:
    subpedigree_base(const string& name, pedigree_id pped);
    subpedigree_base(const string& name, const no_info&);


  private:
    uint        my_index;
    string      my_name;
    pedigree_id my_pedigree;

    //- Sets for lookup of families and members
    member_set      my_members;
    family_set      my_families;

    //- Vectors for random indexing of families and members.
    //
    midx_vector my_member_index;
    fidx_vector my_family_index;

  private:
    subpedigree_base();
    subpedigree_base(const subpedigree_base&);
    void    operator =(const subpedigree_base&);

    void        add_family(family_id f);
    void        add_member(member_id m);
    void        cleanup();
    void        clear_index_arrays();
    void        reset(pedigree_id pid);
    void        reset(const string& name, pedigree_id pid);
    void        set_index(uint i);
    void        build_indices();
};

}
}

#endif

