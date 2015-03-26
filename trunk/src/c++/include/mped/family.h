#ifndef FAMILY_H
#define FAMILY_H

namespace SAGE {
namespace MPED {
  
  /** \ingroup DerivedStorageClasses
  * \brief This is the templatized, derived form of the family_base
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
class family : public family_base
{
  public:

  /** \name Typedefs
    * Typdef-ed local to this class.
    */
  //@{

    typedef GI                                  geninfo_type;
    typedef FI                                  faminfo_type;
    typedef SI                                  subinfo_type;
    typedef PI                                  pedinfo_type;
    typedef MI                                  mpinfo_type;
    typedef family_base                         base;
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
    family() {}

    ///
    /// Destructor.
    ~family() {}

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
    /// Returns the non-const faminfo_type instance associated with this family.
    faminfo_type & info();

    ///
    /// Returns the const faminfo_type instance associated with this family.
    const faminfo_type & info() const;

    ///
    /// Returns the name of this family (see family_base::name() )
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
    /// Returns a non-const pointer to the multipedigree to which this family belongs.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to the multipedigree to which this family belongs.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to the pedigree to which this family belongs.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to the pedigree to which this family belongs.
    pedigree_const_pointer pedigree() const;

    ///
    /// Returns a non-const pointer to the subpedigree to which this family belongs.
    subpedigree_pointer subpedigree();

    ///
    /// Returns a const pointer to the subpedigree to which this family belongs.
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

  /// @name Offspring iterators
  //@{

    ///
    /// Returns a non-const begin iterator to this family's offspring list.
    offspring_iterator offspring_begin();

    ///
    /// Returns a non-const end iterator to this family's offspring list.
    offspring_iterator offspring_end();

    ///
    /// Returns a const begin iterator to this family's offspring list.
    offspring_const_iterator offspring_begin() const;

    ///
    /// Returns a const end iterator to this family's offspring list.
    offspring_const_iterator offspring_end() const;

  //@}

  /// @name Parent iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this family's parent list.
    parent_iterator parent_begin();

    ///
    /// Returns a non-const end iterator for this family's parent list.
    parent_iterator parent_end();

    ///
    /// Returns a const begin iterator for this family's parent list.
    parent_const_iterator parent_begin() const;

    ///
    /// Returns a const end iterator for this family's parent list.
    parent_const_iterator parent_end() const;

  //@}

  private:
    faminfo_type    my_info;

    family(const family&);
    void    operator =(const family&);
};

}
}

#include "mped/family.ipp"

#endif

