#ifndef MEMBER_H
#define MEMBER_H

namespace SAGE {
namespace MPED {

/** \ingroup DerivedStorageClasses
  * \brief This is the templatized, derived form of the member_base
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
class member : public member_base
{
  public:

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef GI                                  geninfo_type;
    typedef FI                                  faminfo_type;
    typedef SI                                  subinfo_type;
    typedef PI                                  pedinfo_type;
    typedef MI                                  mpinfo_type;
    typedef member_base                         base;
    typedef member<GI,FI,SI,PI,MI>        member_type;
    typedef SAGE::MPED::family<GI,FI,SI,PI,MI>        family_type; // See note in sptypes.h about this
    typedef SAGE::MPED::subpedigree<GI,FI,SI,PI,MI>   subpedigree_type; // Ditto
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
    /// \param name The name of the member
    /// \param gender The sex of the member
    /// \param g The geninfo_type instance associated with this individual
    member(const string& name, SexCode gender, const geninfo_type& g);

    ///
    /// Destructor.
    ~member() {}

  //@}

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
    /// Returns \c true if the member is male, \c false if not.
    bool is_male() const;

    ///
    /// Returns \c true if the member is female, \c false if not.
    bool is_female() const;
    
    ///
    /// Returns the index number of the member within the multipedigree.
    uint mpindex() const;

    ///
    /// Returns the index number of the member within the pedigree.
    uint index() const;

    ///
    /// Returns the index number of the member within the subpedigree.
    uint subindex() const;

    ///
    /// Returns the number of mates this member has.
    uint mate_count() const;

    ///
    /// Returns the number of offspring of this member.
    uint offspring_count() const;

    ///
    /// Returns the number of siblings of this member.
    uint sibling_count() const;

    ///
    /// Returns the non-const geninfo_type instance associated with this member.
    geninfo_type & info();

    ///
    /// Returns the const geninfo_type instance associated with this member.
    const geninfo_type & info() const;

    ///
    /// Returns the name of this member.
    const string & name() const;

  //@}
 
  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to the multipedigree to which this member belongs.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to the multipedigree to which this member belongs.
    multipedigree_const_pointer multipedigree() const;

    ///
    /// Returns a non-const pointer to the pedigree to which this member belongs.
    pedigree_pointer pedigree();

    ///
    /// Returns a const pointer to the pedigree to which this member belongs.
    pedigree_const_pointer pedigree() const;

    ///
    /// Returns a non-const pointer to the subpedigree to which this member belongs.
    subpedigree_pointer subpedigree();

    ///
    /// Returns a const pointer to the subpedigree to which this member belongs.
    subpedigree_const_pointer subpedigree() const;

    ///
    /// Returns a non-const pointer to the family to which this member belongs.
    family_pointer family();

    ///
    /// Returns a const pointer to the family to which this member belongs.
    family_const_pointer family() const;

    ///
    /// Returns a non-const pointer to the first parent of this member.
    member_pointer parent1();

    ///
    /// Returns a const pointer to the first parent of this member.
    member_const_pointer parent1() const;

    ///
    /// Returns a non-const pointer to the second parent of this member.
    member_pointer parent2();

    ///
    /// Returns a const pointer to the second parent of this member.
    member_const_pointer parent2() const;

    member_pointer       get_mother();
    member_const_pointer get_mother() const;
    member_pointer       get_father();
    member_const_pointer get_father() const;
    
  //@}

  /// @name Mate iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this member's mate list.
    mate_iterator mate_begin();

    ///
    /// Returns a non-const end iterator for this member's mate list.
    mate_iterator mate_end();

    ///
    /// Returns a const begin iterator for this member's mate list.
    mate_const_iterator mate_begin() const;

    ///
    /// Returns a const end iterator for this member's mate list.
    mate_const_iterator mate_end() const;

  //@}

  /// @name Mate iterators
  //@{

    ///
    /// Given a mate iterator, returns a non-const begin iterator for this member's offspring list.
    /// \param i The member in question
    offspring_iterator offspring_begin(const mate_const_iterator& i);

    ///
    /// Given a member iterator, returns a non-const begin iterator for this member's offspring list.
    /// \param i The member in question
    offspring_iterator offspring_begin(const member_const_iterator& i);

    ///
    /// Given a member, returns a non-const begin iterator for this member's offspring list.
    /// \param m The member in question
    offspring_iterator offspring_begin(const member& m);

    ///
    /// Given a member's name, returns a non-const begin iterator for this member's offspring list.
    /// \param m The member in question
    offspring_iterator offspring_begin(const string& m);

    ///
    /// Returns a non-const end operator for comparison with a member's offspring_iterator.
    ///
    /// Please note that this function is non-parameterized because the offspring_iterator has cursor-like
    /// abilities.
    offspring_iterator offspring_end();

    ///
    /// Given a mate iterator, returns a non-const begin iterator for this member's offspring list.
    /// \param i The member in question
    offspring_const_iterator offspring_begin(const mate_const_iterator& i) const;

    ///
    /// Given a member iterator, returns a non-const begin iterator for this member's offspring list.
    /// \param i The member in question
    offspring_const_iterator offspring_begin(const member_const_iterator& i) const;

    ///
    /// Given a member, returns a non-const begin iterator for this member's offspring list.
    /// \param m The member in question
    offspring_const_iterator offspring_begin(const member& m) const;

    ///
    /// Given a member's name, returns a non-const begin iterator for this member's offspring list.
    /// \param m The member in question
    offspring_const_iterator offspring_begin(const string& m) const;

    ///
    /// Returns a non-const end operator for comparison with a member's offspring_iterator.
    ///
    /// Please note that this function is non-parameterized because the offspring_iterator has cursor-like
    /// abilities.
    offspring_const_iterator offspring_end() const;

  //@}

  /// @name Parent iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this member's parent list.
    parent_iterator parent_begin();

    ///
    /// Returns a non-const end iterator for this member's parent list.
    parent_iterator parent_end();

    ///
    /// Returns a const begin iterator for this member's parent list.
    parent_const_iterator parent_begin() const;

    ///
    /// Returns a const end iterator for this member's parent list.
    parent_const_iterator parent_end() const;

  //@}

  /// @name Progeny iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this member's progeny list.
    progeny_iterator progeny_begin();

    ///
    /// Returns a non-const end iterator for this member's progeny list.
    progeny_iterator progeny_end();

    ///
    /// Returns a const begin iterator for this member's progeny list.
    progeny_const_iterator progeny_begin() const;

    ///
    /// Returns a const end iterator for this member's progeny list.
    progeny_const_iterator progeny_end() const;

  //@}

  /// @name Sibling iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this member's sibling list.
    sibling_iterator sibling_begin();

    ///
    /// Returns a non-const end iterator for this member's sibling list.
    sibling_iterator sibling_end();

    ///
    /// Returns a const begin iterator for this member's sibling list.
    sibling_const_iterator sibling_begin() const;

    ///
    /// Returns a const end iterator for this member's sibling list.
    sibling_const_iterator sibling_end() const;

  //@}

  private:
    geninfo_type    my_info;

    member(const member&);
    void    operator=(const member&);
};

}
}

#include "mped/member.ipp"

#endif
