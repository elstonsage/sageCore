#ifndef PEDIGREE_BASE_H
#define PEDIGREE_BASE_H

#include "mped/spbaseiter.h"
#include "mped/PedigreeBuilder.h"

namespace SAGE {
namespace MPED {

/** \ingroup BaseStorageClasses
  * \brief Represent a single pedigree without any genetic information attached to its indvidual members
  *
  * This class is coupled to its derived class 'pedigree' in
  * order to eliminate bloat, minimize memory usage, and
  * maximize speed.
  */
class pedigree_base
{
  public:
    friend  class PedigreeBuilder;
    friend  class cursor_proxy;
    friend  class subpedigree_base;
    friend  class multipedigree_base;

  /** \name Typedefs
    * Typdef-ed local to this class
    */
  //@{

    typedef error_list::const_iterator          error_iterator;

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
    /// \param name The name of the pedigree
    pedigree_base(const string& name);

    ///
    /// Destructor
    virtual ~pedigree_base();

  //@}

  /// @name Adding members & relationships
  //@{

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param child The name of the child
    /// \param parent The name of the parent
    void add_lineage(const string& child, const string& parent);

    ///
    /// Amends the internal pedigree storage to reflect a parent-child relationship.
    /// \param child The name of the child
    /// \param parent1 The name of the first parent
    /// \param parent2 The name of the second parent
    void add_lineage(const string& child, const string& parent1, const string& parent2);

    ///
    /// Amends the internal pedigree storage to reflect a spousal relationship.
    /// \param spouse1 The name of the first spouse
    /// \param spouse2 The name of the second spouse
    void add_marriage(const string& spouse1, const string& spouse2);

    ///
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param name The name of the new member
    /// \param s The sex of the new member
    void add_member(const string& name, SexCode s=SEX_MISSING);

    ///
    /// Amends the internal pedigree storage to reflect a sibling relationship.
    /// \param sib1 The name of the first sibling
    /// \param sib2 The name of the second sibling
    void add_sibship(const string& sib1, const string& sib2);

  //@}

  /// @name Building functions
  //@{

    /// 
    /// Having added entries for all members as well as structural relationships, this function finalizes
    /// the internal storage for this pedigree.
    void build();

    ///
    /// Having invoked build(), you can then use freeze() to disallow any additional changes to the
    /// data. Any subsequent attempt to alter the data will cause an error.
    void freeze();

  //@}

  /// @name Attributes
  //@{

    ///
    /// Returns the index number of this pedigree within the multipedigree_base.
    uint index() const;

    ///
    /// Returns the number of families present in this pedigree.
    uint family_count() const;

    ///
    /// Returns the number of members present in this pedigree.
    uint member_count() const;

    ///
    /// Returns the number of subpedigrees present in this pedigree.
    uint subpedigree_count() const;

    ///
    /// Returns the number of unconnected individuals (singletons) present in this pedigree.
    uint unconnected_count() const;

    ///
    /// Returns the name of this pedigree (as indicated in the original source data).
    const string & name() const;

  //@}

  /// @name Constituent pointers
  //@{

    ///
    /// Returns a non-const pointer to this pedigree_base's multipedigree_base.
    multipedigree_pointer multipedigree();

    ///
    /// Returns a const pointer to this pedigree_base's multipedigree_base.
    multipedigree_const_pointer multipedigree() const;

  //@}

  /// @name Error information
  //@{

    ///
    /// Returns the number of errors generated in the course of importing data into this pedigree_base.
    uint error_count() const;

    ///
    /// Returns a non-const begin iterator for this pedigree_base's error list.
    error_iterator error_begin() const;

    ///
    /// Returns a non-const end iterator for this pedigree_base's error list.
    error_iterator error_end() const;

  //@}

  /// @name Warning information
  //@{

    ///
    /// Returns the number of warnings generated in the course of importing data into this pedigree_base.
    uint warning_count() const;

    ///
    /// Returns a non-const begin iterator for this pedigree_base's warning list.
    error_iterator warning_begin() const;

    ///
    /// Returns a non-const end iterator for this pedigree_base's warning list.
    error_iterator warning_end() const;

  //@}

  /// @name Index-based lookup functions
  //@{

    ///
    /// Given a family index number, returns the corresponding non-const family reference.
    /// \param i The index number of the family
    family_base& family_index(uint i);

    ///
    /// Given a family index number, returns the corresponding const family reference.
    /// \param i The index number of the family
    const family_base & family_index(uint i) const;

    ///
    /// Given a member index number, returns the corresponding non-const member reference.
    /// \param i The index number of the member
    member_base& member_index(uint i);

    ///
    /// Given a member index number, returns the corresponding const member reference.
    /// \param i The index number of the member
    const member_base & member_index(uint i) const;

    ///
    /// Given a subpedigree_base index number, returns the corresponding non-const subpedigree_base reference.
    /// \param i The index number of the subpedigree_base
    subpedigree_base& subpedigree_index(uint i);

    ///
    /// Given a subpedigree_base index number, returns the corresponding const subpedigree_base reference.
    /// \param i The index number of the subpedigree_base
    const subpedigree_base & subpedigree_index(uint i) const;

  //@}

  /// @name Index swapping functions
  //@{

    ///
    /// Swaps the index numbers of the two indicated families.
    /// \param i The index number of the first family
    /// \param j The index number of the second family
    void family_index_swap(uint i, uint j);

    ///
    /// Swaps the index numbers of the two indicated members.
    /// \param i The index number of the first member
    /// \param j The index number of the second member
    void member_index_swap(uint i, uint j);

    ///
    /// Swaps the index numbers of the two indicated subpedigrees
    /// \param i The index number of the first subpedigree
    /// \param j The index number of the second subpedigree
    void subpedigree_index_swap(uint i, uint j);

  //@}

  /// @name Name-, member-reference-, and member-iterator-based lookup functions
  //@{

    ///
    /// Given the names of the two constituent parents, return a non-const pointer to the corresponding family.
    /// \param p1 The name of the first parent
    /// \param p2 The name of the second parent
    family_pointer family_find(const string& p1, const string& p2);

    ///
    /// Given the names of the two constituent parents, return a const pointer to the corresponding family.
    /// \param p1 The name of the first parent
    /// \param p2 The name of the second parent
    family_const_pointer family_find(const string& p1, const string& p2) const;

    ///
    /// Given the member_base objects that represent the two constituent parents, returns a non-const pointer to the 
    /// corresponding family.
    /// \param p1 The member_base of the first parent
    /// \param p2 The member_base of the second parent
    family_pointer family_find(const member_base& p1, const member_base& p2);

    ///
    /// Given the member_base objects that represent the two constituent parents, returns a const pointer to the 
    /// corresponding family.
    /// \param p1 The member_base of the first parent
    /// \param p2 The member_base of the second parent
    family_const_pointer family_find(const member_base& p1, const member_base& p2) const;

    ///
    /// Given the member_base iterators that point to the two constituent parents, returns a 
    /// non-const pointer to the corresponding family.
    /// \param p1 The member_const_iterator of the first parent
    /// \param p2 The member_const_iterator of the second parent
    family_pointer family_find(const member_const_iterator& p1, const member_const_iterator& p2);

    ///
    /// Given the member_base iterators that point to the two constituent parents, returns a 
    /// const pointer to the corresponding family.
    /// \param p1 The member_const_iterator of the first parent
    /// \param p2 The member_const_iterator of the second parent
    family_const_pointer family_find(const member_const_iterator& p1, const member_const_iterator& p2) const;

    ///
    /// Given a member's name, returns the corresponding non-const member pointer.
    /// \param name The name of the member
    member_pointer member_find(const string& name);

    ///
    /// Given a member's name, returns the corresponding const member pointer.
    /// \param name The name of the member
    member_const_pointer member_find(const string& name) const;

    ///
    /// Given a subpedigree's name, returns the corresponding non-const subpedigree pointer.
    /// \param name The name of the subpedigree
    subpedigree_pointer subpedigree_find(const string& name);

    ///
    /// Given a subpedigree's name, returns the corresponding const subpedigree pointer.
    /// \param name The name of the subpedigree
    subpedigree_const_pointer subpedigree_find(const string& name) const;

  //@}

  /// @name Family iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this pedigree's family list.
    family_iterator family_begin();

    ///
    /// Returns a non-const end iterator for this pedigree's family list.
    family_iterator family_end();

    ///
    /// Given a subpedigree, returns a non-const begin iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_iterator family_begin(const subpedigree_base& s);

    ///
    /// Given a subpedigree, returns a non-const end iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_iterator family_end(const subpedigree_base& s);

    ///
    /// Given a subpedigree iterator, returns a non-const begin iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_iterator family_begin(const subpedigree_const_iterator& s);

    ///
    /// Given a subpedigree iterator, returns a non-const end iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_iterator family_end(const subpedigree_const_iterator& s);

    ///
    /// Returns a const begin iterator for this pedigree's family list.
    family_const_iterator family_begin() const;

    ///
    /// Returns a const end iterator for this pedigree's family list.
    family_const_iterator family_end() const;

    ///
    /// Given a subpedigree, returns a const begin iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_const_iterator family_begin(const subpedigree_base& s) const;

    ///
    /// Given a subpedigree, returns a const end iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_const_iterator family_end(const subpedigree_base& s) const;

    ///
    /// Given a subpedigree iterator, returns a const begin iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_const_iterator family_begin(const subpedigree_const_iterator& s) const;

    ///
    /// Given a subpedigree iterator, returns a const end iterator for that subpedigree's family list.
    /// \param s The subpedigree in question
    family_const_iterator family_end(const subpedigree_const_iterator& s) const;

  //@}
  /// @name Mate iterators
  //@{

    ///
    /// Given a member's name, returns a non-const begin iterator for that member's mate list.
    /// \param m The name of the member
    mate_iterator mate_begin(const string & m);

    ///
    /// Given a member's name, returns a non-const end iterator for that member's mate list.
    /// \param m The name of the member
    mate_iterator mate_end(const string& m);

    ///
    /// Given a member, returns a non-const begin iterator for that member's mate list.
    /// \param m The member in question
    mate_iterator mate_begin(const member_base & m);

    ///
    /// Given a member, returns a non-const end iterator for that member's mate list.
    /// \param m The member in question
    mate_iterator mate_end(const member_base& m);

    ///
    /// Given a member iterator, returns a non-const begin iterator for that member's mate list.
    /// \param m The member in question
    mate_iterator mate_begin(const member_const_iterator & m);

    ///
    /// Given a member iterator, returns a non-const end iterator for that member's mate list.
    /// \param m The member in question
    mate_iterator mate_end(const member_const_iterator& m);

    ///
    /// Given a member's name, returns a const begin iterator for that member's mate list.
    /// \param m The name of the member
    mate_const_iterator mate_begin(const string & m) const;

    ///
    /// Given a member's name, returns a const end iterator for that member's mate list.
    /// \param m The name of the member
    mate_const_iterator mate_end(const string& m) const;

    ///
    /// Given a member, returns a const begin iterator for that member's mate list.
    /// \param m The member in question
    mate_const_iterator mate_begin(const member_base & m) const;

    ///
    /// Given a member, returns a const end iterator for that member's mate list.
    /// \param m The member in question
    mate_const_iterator mate_end(const member_base& m) const;

    ///
    /// Given a member iterator, returns a const begin iterator for that member's mate list.
    /// \param m The member in question
    mate_const_iterator mate_begin(const member_const_iterator & m) const;

    ///
    /// Given a member iterator, returns a const end iterator for that member's mate list.
    /// \param m The member in question
    mate_const_iterator mate_end(const member_const_iterator& m) const;

  //@}

  /// @name Member iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this pedigree's member list.
    member_iterator member_begin();

    ///
    /// Returns a non-const end iterator for this pedigree's member list.
    member_iterator member_end();

    ///
    /// Given a subpedigree, returns a non-const begin iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_iterator member_begin(const subpedigree_base& s);

    ///
    /// Given a subpedigree, returns a non-const end iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_iterator member_end(const subpedigree_base& s);

    ///
    /// Given a subpedigree iterator, returns a non-const begin iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_iterator member_begin(const subpedigree_const_iterator& s);

    ///
    /// Given a subpedigree iterator, returns a non-const end iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_iterator member_end(const subpedigree_const_iterator& s);

    ///
    /// Returns a const begin iterator for this pedigree's member list.
    member_const_iterator member_begin() const;

    ///
    /// Returns a const end iterator for this pedigree's member list.
    member_const_iterator member_end() const;

    ///
    /// Given a subpedigree, returns a const begin iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_const_iterator member_begin(const subpedigree_base& s) const;

    ///
    /// Given a subpedigree, returns a const end iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_const_iterator member_end(const subpedigree_base& s) const;

    ///
    /// Given a subpedigree iterator, returns a const begin iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_const_iterator member_begin(const subpedigree_const_iterator& s) const;

    ///
    /// Given a subpedigree iterator, returns a const end iterator for that subpedigree's member list.
    /// \param s The subpedigree in question
    member_const_iterator member_end(const subpedigree_const_iterator& s) const;

  //@}

  /// @name Offspring iterators
  //@{

    ///
    /// Given a family, returns a non-const begin iterator for that family's offspring list.
    /// \param f The family in question
    offspring_iterator offspring_begin(const family_base& f);

    ///
    /// Given a family iterator, returns a non-const begin iterator for that family's offspring list.
    /// \param f The family in question
    offspring_iterator offspring_begin(const family_const_iterator& f);

    ///
    /// Given the member objects of two constituent parents, returns the corresponding non-const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The first parent in question
    /// \param p2 The second parent in question
    offspring_iterator offspring_begin(const member_base& p1, const member_base& p2);

    ///
    /// Given the member iterators of two constituent parents, returns the corresponding non-const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The first parent in question
    /// \param p2 The second parent in question
    offspring_iterator offspring_begin(const member_const_iterator& p1, const member_const_iterator& p2);

    ///
    /// Given the names of two constituent parents, returns the corresponding non-const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The name of the first parent in question
    /// \param p2 The name of the second parent in question
    offspring_iterator offspring_begin(const string& p1, const string& p2);

    ///
    /// Returns the non-const end iterator for comparison with a offspring_iterator.
    ///
    /// Please note that this function is non-parameterized because the offspring_iterator has cursor-like
    /// abilities.
    offspring_iterator offspring_end();

    ///
    /// Given a family, returns a const begin iterator for that family's offspring list.
    /// \param f The family in question
    offspring_const_iterator offspring_begin(const family_base& f) const;

    ///
    /// Given a family iterator, returns a const begin iterator for that family's offspring list.
    /// \param f The family in question
    offspring_const_iterator offspring_begin(const family_const_iterator& f) const;

    ///
    /// Given the member objects of two constituent parents, returns the corresponding const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The first parent in question
    /// \param p2 The second parent in question
    offspring_const_iterator offspring_begin(const member_base& p1, const member_base& p2) const;

    ///
    /// Given the member iterators of two constituent parents, returns the corresponding const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The first parent in question
    /// \param p2 The second parent in question
    offspring_const_iterator offspring_begin(const member_const_iterator& p1, const member_const_iterator& p2) const;

    ///
    /// Given the names of two constituent parents, returns the corresponding const
    /// begin iterator for those parents' offspring list.
    /// \param p1 The name of the first parent in question
    /// \param p2 The name of the second parent in question
    offspring_const_iterator offspring_begin(const string& p1, const string& p2) const;

    ///
    /// Returns the const end iterator for comparison with a offspring_iterator.
    ///
    /// Please note that this function is non-parameterized because the offspring_iterator has cursor-like
    /// abilities.
    offspring_const_iterator offspring_end() const;

  //@}

  /// @name Parent iterators
  //@{

    ///
    /// Given a member's name, returns a non-const begin iterator for that member's parent list.
    /// \param name The name of the member in question
    parent_iterator parent_begin(const string& name);

    ///
    /// Given a member, returns a non-const begin iterator for that member's parent list.
    /// \param m The member in question
    parent_iterator parent_begin(const member_base& m);

    ///
    /// Given a member iterator, returns a non-const begin iterator for that member's parent list.
    /// \param m The member in question
    parent_iterator parent_begin(const member_const_iterator& m);

    ///
    /// Given a family, returns a non-const begin iterator for that family's parent list.
    /// \param f The family in question
    parent_iterator parent_begin(const family_base& f);

    ///
    /// Given a family iterator, returns a non-const begin iterator for that family's parent list.
    /// \param f The family in question
    parent_iterator parent_begin(const family_const_iterator& f);

    ///
    /// Returns the non-const end iterator for comparison with a parent_iterator.
    ///
    /// Please note that this function is non-parameterized because the parent_iterator has cursor-like
    /// abilities.
    parent_iterator parent_end();

    ///
    /// Given a member's name, returns a const begin iterator for that member's parent list.
    /// \param name The name of the member in question
    parent_const_iterator parent_begin(const string& name) const;

    ///
    /// Given a member, returns a const begin iterator for that member's parent list.
    /// \param m The member in question
    parent_const_iterator parent_begin(const member_base& m) const;

    ///
    /// Given a member iterator, returns a const begin iterator for that member's parent list.
    /// \param m The member in question
    parent_const_iterator parent_begin(const member_const_iterator& m) const;

    ///
    /// Given a family, returns a const begin iterator for that family's parent list.
    /// \param f The family in question
    parent_const_iterator parent_begin(const family_base& f) const;

    ///
    /// Given a family iterator, returns a const begin iterator for that family's parent list.
    /// \param f The family in question
    parent_const_iterator parent_begin(const family_const_iterator& f) const;

    ///
    /// Returns the const end iterator for comparison with a parent_iterator.
    ///
    /// Please note that this function is non-parameterized because the parent_iterator has cursor-like
    /// abilities.
    parent_const_iterator parent_end() const;

  //@}

  /// @name Progeny iterators
  //@{

    ///
    /// Given a member's name, returns a non-const begin iterator for that member's progeny list.
    /// \param p1 The name of the member in question
    progeny_iterator progeny_begin(const string& p1);

    ///
    /// Given a member's name, returns a non-const end iterator for that member's progeny list.
    /// \param p1 The name of the member in question
    progeny_iterator progeny_end(const string& p1);

    ///
    /// Given a member, returns a non-const begin iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_iterator progeny_begin(const member_base& p1);

    ///
    /// Given a member, returns a non-const end iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_iterator progeny_end(const member_base& p1);

    ///
    /// Given a member iterator, returns a non-const begin iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_iterator progeny_begin(const member_const_iterator& p1);

    ///
    /// Given a member iterator, returns a non-const end iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_iterator progeny_end(const member_const_iterator& p1);

    ///
    /// Given a member's name, returns a const begin iterator for that member's progeny list.
    /// \param p1 The name of the member in question
    progeny_const_iterator progeny_begin(const string& p1) const;

    ///
    /// Given a member's name, returns a const end iterator for that member's progeny list.
    /// \param p1 The name of the member in question
    progeny_const_iterator progeny_end(const string& p1) const;

    ///
    /// Given a member, returns a const begin iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_const_iterator progeny_begin(const member_base& p1) const;

    ///
    /// Given a member, returns a const end iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_const_iterator progeny_end(const member_base& p1) const;

    ///
    /// Given a member iterator, returns a const begin iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_const_iterator progeny_begin(const member_const_iterator& p1) const;

    ///
    /// Given a member iterator, returns a const end iterator for that member's progeny list.
    /// \param p1 The member in question
    progeny_const_iterator progeny_end(const member_const_iterator& p1) const;

  //@}

  /// @name Sibling iterators
  //@{

    ///
    /// Given a member's name, returns a non-const begin iterator for that member's sibling list.
    /// \param name The name of the member in question
    sibling_iterator sibling_begin(const string& name);

    ///
    /// Given a member, returns a non-const begin iterator for that member's sibling list.
    /// \param m The name of the member in question
    sibling_iterator sibling_begin(const member_base& m);

    ///
    /// Given a member iterator, returns a non-const begin iterator for that member's sibling list.
    /// \param m The name of the member in question
    sibling_iterator sibling_begin(const member_const_iterator& m);

    ///
    /// Returns the non-const end iterator for comparison with a sibling_iterator.
    ///
    /// Please note that this function is non-parameterized because the sibling_iterator has cursor-like
    /// abilities.
    sibling_iterator sibling_end();

    ///
    /// Given a member's name, returns a const begin iterator for that member's sibling list.
    /// \param name The name of the member in question
    sibling_const_iterator sibling_begin(const string& name) const;

    ///
    /// Given a member, returns a const begin iterator for that member's sibling list.
    /// \param m The name of the member in question
    sibling_const_iterator sibling_begin(const member_base& m) const;

    ///
    /// Given a member iterator, returns a const begin iterator for that member's sibling list.
    /// \param m The name of the member in question
    sibling_const_iterator sibling_begin(const member_const_iterator& m) const;

    ///
    /// Returns the const end iterator for comparison with a sibling_iterator.
    ///
    /// Please note that this function is non-parameterized because the sibling_iterator has cursor-like
    /// abilities.
    sibling_const_iterator sibling_end() const;

  //@}

  /// @name Subpedigree iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this pedigree's subpedigree list.
    subpedigree_iterator subpedigree_begin();

    ///
    /// Returns a non-const end iterator for this pedigree's subpedigree list.
    subpedigree_iterator subpedigree_end();

    ///
    /// Returns a const begin iterator for this pedigree's subpedigree list.
    subpedigree_const_iterator subpedigree_begin() const;

    ///
    /// Returns a const end iterator for this pedigree's subpedigree list.
    subpedigree_const_iterator subpedigree_end() const;

  //@}

  /// @name Singleton iterators
  //@{

    ///
    /// Returns a non-const begin iterator for this pedigree's singleton list.
    member_iterator unconnected_begin();

    ///
    /// Returns a non-const end iterator for this pedigree's singleton list.
    member_iterator unconnected_end();

    ///
    /// Returns a const begin iterator for this pedigree's singleton list.
    member_const_iterator unconnected_begin() const;

    ///
    /// Returns a const end iterator for this pedigree's singleton list.
    member_const_iterator unconnected_end() const;

  //@}

  protected:
    pedigree_base(const string& name, multipedigree_pointer mp);

    //- Initialization and destruction functions.
    //
    void    init();
    void    clear_all();


    //- Lookup functions.
    //
    fset_iterator       ifind_family(const string& parent1, const string& parent2) const;
    fset_iterator       ifind_family(member_id par1, member_id par2) const;
    mset_iterator       ifind_member(const string& name) const;
    sset_iterator       ifind_subped(const string& name) const;

    family_id           lookup_family(const string& parent1, const string& parent2) const;
    family_id           lookup_family(member_id par1, member_id par2) const;
    member_id           lookup_member(const string& name) const;
    subped_id           lookup_subped(const string& name) const;

    //- These functions implement type-independent functionality that has
    //  been deferred to the derived class.
    //
    virtual subped_id   allocate_subped(const string& name, pedigree_id pped);
    virtual family_id   allocate_family();
    virtual member_id   allocate_member(const string& name, SexCode x, const void* pdata);

    virtual void        deallocate_subped(subped_id id);
    virtual void        deallocate_family(family_id id);
    virtual void        deallocate_member(member_id id);

    //- Functions for adding elements to the pedigree graph.
    //
    family_id           add_family(member_id P1, member_id P2, member_id C);
    member_id           add_member(const string& name, SexCode s, const void* p);
    subped_id           add_subped();


  private:
    //- Pedigree name and state information.
    //
    uint            my_index;
    uint            my_fcount;
    uint            my_mcount;
    uint            my_scount;
    uint            my_ucount;
    string          my_name;
    bool            my_readiness;
    bool            my_fstate;
    subped_id       my_unconnecteds;
    multiped_id     my_mped;

    //- Subpedigrees, members, nuclear families, and mating chains.
    //
    subped_set      my_subpeds;
    member_set      my_members;
    family_set      my_families;
    mate_multimap   my_matechains;

    //- Vectors for random indexing of families, members, and subpedigrees.
    //
    sidx_vector     my_subped_index;
    midx_vector     my_member_index;
    fidx_vector     my_family_index;

    //
    PedigreeBuilder my_builder; ///< Builder object for this pedigree

    private:
    pedigree_base(const pedigree_base&);
    void    operator=(const pedigree_base&);

    void        set_index(uint i);

    /// 
    /// Clears all entries from this pedigree.
    void clear();

    void flush_build_info();
};

}
}

#endif

