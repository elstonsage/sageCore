#ifndef SUBPEDIGREE_H
#define SUBPEDIGREE_H

namespace SAGE {
namespace MPED {
  
/** \ingroup DerivedStorageClasses
  * \brief This is the templatized, derived form of the subpedigree_base
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  */
template <class GI, class FI=no_info, class SI=no_info, class PI=no_info, class MI=no_info>
class subpedigree : public subpedigree_base
{
  public:

  /** \name Typdefs
    * Typdef-ed local to this class
    */
  //@{

    typedef GI                                  geninfo_type;
    typedef FI                                  faminfo_type;
    typedef SI                                  subinfo_type;
    typedef PI                                  pedinfo_type;
    typedef MI                                  mpinfo_type;
    typedef subpedigree_base                    base;
    typedef SAGE::MPED::member<GI,FI,SI,PI,MI>        member_type;
    typedef SAGE::MPED::family<GI,FI,SI,PI,MI>        family_type;
    typedef SAGE::MPED::subpedigree<GI,FI,SI,PI,MI>   subpedigree_type;
    typedef SAGE::MPED::pedigree<GI,FI,SI,PI,MI>      pedigree_type;
    typedef SAGE::MPED::multipedigree<GI,FI,SI,PI,MI> multipedigree_type;

    typedef member_type*                                member_pointer;
    typedef family_type*                                family_pointer;
    typedef subpedigree_type*                           subpedigree_pointer;
    typedef pedigree_type*                              pedigree_pointer;
    typedef multipedigree_type*                         multipedigree_pointer;
    typedef SAGE::MPED::family_iterator<GI,FI,SI,PI,MI>       family_iterator;
    typedef SAGE::MPED::mate_iterator<GI,FI,SI,PI,MI>         mate_iterator;
    typedef SAGE::MPED::member_iterator<GI,FI,SI,PI,MI>       member_iterator;
    typedef SAGE::MPED::offspring_iterator<GI,FI,SI,PI,MI>    offspring_iterator;
    typedef SAGE::MPED::parent_iterator<GI,FI,SI,PI,MI>       parent_iterator;
    typedef SAGE::MPED::pedigree_iterator<GI,FI,SI,PI,MI>     pedigree_iterator;
    typedef SAGE::MPED::progeny_iterator<GI,FI,SI,PI,MI>      progeny_iterator;
    typedef SAGE::MPED::sibling_iterator<GI,FI,SI,PI,MI>      sibling_iterator;
    typedef SAGE::MPED::subpedigree_iterator<GI,FI,SI,PI,MI>  subpedigree_iterator;

    typedef const member_type*                                  member_const_pointer;
    typedef const family_type*                                  family_const_pointer;
    typedef const subpedigree_type*                             subpedigree_const_pointer;
    typedef const pedigree_type*                                pedigree_const_pointer;
    typedef const multipedigree_type*                           multipedigree_const_pointer;
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

  public:

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor.
    /// \param name The name of the subpedigree
    /// \param pped The pedigree_base of which this subpedigree is a member
    subpedigree(const string& name, pedigree_base* pped);

    ///
    /// Destructor.
    ~subpedigree(){}

  //@}

  /// @name Attributes
  //@{

    ///
    /// Returns the index number of this subpedigree within the pedigree.
    uint index() const;

    ///
    /// Returns the number of families in this subpedigree.
    uint family_count() const;

    ///
    /// Returns the number of members in this subpedigree.
    uint member_count() const;

    ///
    /// Returns the non-const subinfo_type instance associated with this subpedigree.
    subinfo_type & info();

    ///
    /// Returns the const subinfo_type instance associated with this subpedigree.
    const subinfo_type & info() const;

    ///
    /// Returns this subpedigree's name (see subpedigree_base::name() ).
    const string & name() const;

  //@}

  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to the multipedigree to which this subpedigree belongs.
    multipedigree_pointer   multipedigree();

    ///
    /// Returns a const pointer to the multipedigree to which this subpedigree belongs.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to the pedigree to which this subpedigree belongs.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to the pedigree to which this subpedigree belongs.
    pedigree_const_pointer pedigree() const;

  //@}

  /// @name Index-based lookup functions
  //@{

    ///
    /// Given a family's id number, return the corresponding non-const family reference.
    /// \param i The family's id number
    family_type & family_index(uint i);

    ///
    /// Given a family's id number, return the corresponding const family reference.
    /// \param i The family's id number
    const family_type & family_index(uint i) const;

    ///
    /// Given a member's id number, return the corresponding non-const member reference.
    /// \param i The member's id number
    member_type & member_index(uint i);

    ///
    /// Given a member's id number, return the corresponding const member reference.
    /// \param i The member's id number
    const member_type & member_index(uint i) const;

  //@}

  /// @name Index swapping functions
  //@{

    ///
    /// Swaps the index numbers of two families.
    /// \param i The index number of the first family
    /// \param j The index number of the second family
    void family_index_swap(uint i, uint j);

    ///
    /// Swaps the index numbers of two members.
    /// \param i The index number of the first member
    /// \param j The index number of the second member
    void member_index_swap(uint i, uint j);

  //@}

  /// @name Family iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this subpedigree's family list.
    family_iterator family_begin();

    ///
    /// Returns a non-const end iterator for this subpedigree's family list.
    family_iterator family_end();

    ///
    /// Returns a const begin iterator for this subpedigree's family list.
    family_const_iterator family_begin() const;

    ///
    /// Returns a const end iterator for this subpedigree's family list.
    family_const_iterator family_end() const;

  //@}

  /// @name Member iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this subpedigree's member list.
    member_iterator member_begin();

    ///
    /// Returns a non-const end iterator for this subpedigree's member list.
    member_iterator member_end();

    ///
    /// Returns a const begin iterator for this subpedigree's member list.
    member_const_iterator member_begin() const;

    ///
    /// Returns a const end iterator for this subpedigree's member list.
    member_const_iterator member_end() const;

  //@}

  private:
    subinfo_type    my_info;

    subpedigree(const subpedigree&);
    void    operator =(const subpedigree&);
};

}
}

#include "mped/subpedigree.ipp"

#endif

