#ifndef PEDIGREE_H
#define PEDIGREE_H

namespace SAGE {
namespace MPED {

/** \ingroup DerivedStorageClasses
  * \brief This is the templatized, derived form of the pedigree_base
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
class pedigree : public pedigree_base
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
    typedef pedigree_base                       base;
    typedef SAGE::MPED::member<GI,FI,SI,PI,MI>        member_type;
    typedef SAGE::MPED::family<GI,FI,SI,PI,MI>        family_type;
    typedef SAGE::MPED::subpedigree<GI,FI,SI,PI,MI>   subpedigree_type;
    typedef SAGE::MPED::pedigree<GI,FI,SI,PI,MI>      pedigree_type;
    typedef SAGE::MPED::multipedigree<GI,FI,SI,PI,MI> multipedigree_type;
    typedef base::error_iterator                error_iterator;

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

    typedef void (* callback_type)(pedigree<GI,FI,SI,PI,MI>&);

  //@}

  public:

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor
    /// \param name The name of the pedigree
    pedigree(const string& name);

    ///
    /// Constructor
    /// \param name The name of the pedigree
    /// \param mp The multipedigree to which this pedigree belongs
    pedigree(const string& name, multipedigree_type& mp);

    ///
    /// Destructor
    virtual ~pedigree();

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
    /// Adds an entry for a new member in the internal pedigree storage.
    /// \param name The name of the new member
    /// \param s The sex of the new member
    /// \param info The geninfo_type instance associated with this member
    void add_member(const string& name, SexCode s, const geninfo_type& info);

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
    /// Returns the index number of this pedigree within the multipedigree.
    uint index() const;


    uint family_count() const;

    uint member_count() const;

    uint subpedigree_count() const;

    uint unconnected_count() const;

    pedinfo_type & info();

    const pedinfo_type & info() const;

    const string & name() const;

  //@}

    uint error_count() const;

  /// @name Constituent pointers
  //@{

    multipedigree_pointer   multipedigree();

    multipedigree_const_pointer multipedigree() const;

  //@}

  /// @name Index-based lookup functions
  //@{

    const family_type&          family_index(uint i) const;
    const member_type&          member_index(uint i) const;
    const subpedigree_type&     subpedigree_index(uint i) const;

    family_type&            family_index(uint i);
    member_type&            member_index(uint i);
    subpedigree_type&       subpedigree_index(uint i);

  //@}

  /// @name Index swapping functions
  //@{

    void    family_index_swap(uint i, uint j);
    void    member_index_swap(uint i, uint j);
    void    subpedigree_index_swap(uint i, uint j);

  //@}

  /// @name Name-, member-reference-, and member-iterator-based lookup functions
  //@{

    family_const_pointer        family_find(const string& p1, const string& p2) const;
    family_const_pointer        family_find(const member_type& p1, const member_type& p2) const;
    family_const_pointer        family_find(const member_const_iterator& p1, const member_const_iterator& p2) const;
    member_const_pointer        member_find(const string& name) const;
    subpedigree_const_pointer   subpedigree_find(const string& name) const;

    family_pointer          family_find(const string& p1, const string& p2);
    family_pointer          family_find(const member_type& p1, const member_type& p2);
    family_pointer          family_find(const member_const_iterator& p1, const member_const_iterator& p2);
    member_pointer          member_find(const string& name);
    subpedigree_pointer     subpedigree_find(const string& name);

  //@}

    //- Const iteration.
    //
    error_iterator              error_begin() const;
    error_iterator              error_end() const;

    family_iterator         family_begin();
    family_iterator         family_begin(const subpedigree_type& s);
    family_iterator         family_begin(const subpedigree_const_iterator& s);
    family_iterator         family_end();
    family_iterator         family_end(const subpedigree_type& s);
    family_iterator         family_end(const subpedigree_const_iterator& s);
    family_const_iterator       family_begin() const;
    family_const_iterator       family_begin(const subpedigree_type& s) const;
    family_const_iterator       family_begin(const subpedigree_const_iterator& s) const;
    family_const_iterator       family_end() const;
    family_const_iterator       family_end(const subpedigree_type& s) const;
    family_const_iterator       family_end(const subpedigree_const_iterator& s) const;

    mate_iterator           mate_begin(const string& m);
    mate_iterator           mate_begin(const member_type& m);
    mate_iterator           mate_begin(const member_const_iterator& m);
    mate_iterator           mate_end(const string& m);
    mate_iterator           mate_end(const member_type& m);
    mate_iterator           mate_end(const member_const_iterator& m);
    mate_const_iterator         mate_begin(const string& m) const;
    mate_const_iterator         mate_begin(const member_type& m) const;
    mate_const_iterator         mate_begin(const member_const_iterator& m) const;
    mate_const_iterator         mate_end(const string& m) const;
    mate_const_iterator         mate_end(const member_type& m) const;
    mate_const_iterator         mate_end(const member_const_iterator& m) const;

    member_iterator         member_begin();
    member_iterator         member_begin(const subpedigree_type& s);
    member_iterator         member_begin(const subpedigree_const_iterator& s);
    member_iterator         member_end();
    member_iterator         member_end(const subpedigree_type& s);
    member_iterator         member_end(const subpedigree_const_iterator& s);
    member_const_iterator       member_begin() const;
    member_const_iterator       member_begin(const subpedigree_type& s) const;
    member_const_iterator       member_begin(const subpedigree_const_iterator& s) const;
    member_const_iterator       member_end() const;
    member_const_iterator       member_end(const subpedigree_type& s) const;
    member_const_iterator       member_end(const subpedigree_const_iterator& s) const;

    offspring_iterator      offspring_begin(const family_type& f);
    offspring_iterator      offspring_begin(const family_const_iterator& f);
    offspring_iterator      offspring_begin(const member_type& p1, const member_type& p2);
    offspring_iterator      offspring_begin(const member_const_iterator& p1, const member_const_iterator& p2);
    offspring_iterator      offspring_begin(const string& p1, const string& p2);
    offspring_iterator      offspring_end();
    offspring_const_iterator    offspring_begin(const family_type& f) const;
    offspring_const_iterator    offspring_begin(const family_const_iterator& f) const;
    offspring_const_iterator    offspring_begin(const member_type& p1, const member_type& p2) const;
    offspring_const_iterator    offspring_begin(const member_const_iterator& p1, const member_const_iterator& p2) const;
    offspring_const_iterator    offspring_begin(const string& p1, const string& p2) const;
    offspring_const_iterator    offspring_end() const;

    parent_iterator         parent_begin(const string& name);
    parent_iterator         parent_begin(const member_type& m);
    parent_iterator         parent_begin(const member_const_iterator& m);
    parent_iterator         parent_begin(const family_type& f);
    parent_iterator         parent_begin(const family_const_iterator& f);
    parent_iterator         parent_end();
    parent_const_iterator       parent_begin(const string& name) const;
    parent_const_iterator       parent_begin(const member_type& m) const;
    parent_const_iterator       parent_begin(const member_const_iterator& m) const;
    parent_const_iterator       parent_begin(const family_type& f) const;
    parent_const_iterator       parent_begin(const family_const_iterator& f) const;
    parent_const_iterator       parent_end() const;

    progeny_iterator        progeny_begin(const string& p1);
    progeny_iterator        progeny_begin(const member_type& p1);
    progeny_iterator        progeny_begin(const member_const_iterator& p1);
    progeny_iterator        progeny_end(const string& p1);
    progeny_iterator        progeny_end(const member_type& p1);
    progeny_iterator        progeny_end(const member_const_iterator& p1);
    progeny_const_iterator      progeny_begin(const string& p1) const;
    progeny_const_iterator      progeny_begin(const member_type& p1) const;
    progeny_const_iterator      progeny_begin(const member_const_iterator& p1) const;
    progeny_const_iterator      progeny_end(const string& p1) const;
    progeny_const_iterator      progeny_end(const member_type& p1) const;
    progeny_const_iterator      progeny_end(const member_const_iterator& p1) const;

    sibling_iterator        sibling_begin(const string& name);
    sibling_iterator        sibling_begin(const member_type& m);
    sibling_iterator        sibling_begin(const member_const_iterator& m);
    sibling_iterator        sibling_end();
    sibling_const_iterator      sibling_begin(const string& name) const;
    sibling_const_iterator      sibling_begin(const member_type& m) const;
    sibling_const_iterator      sibling_begin(const member_const_iterator& m) const;
    sibling_const_iterator      sibling_end() const;

    subpedigree_iterator    subpedigree_begin();
    subpedigree_iterator    subpedigree_end();
    subpedigree_const_iterator  subpedigree_begin() const;
    subpedigree_const_iterator  subpedigree_end() const;

    member_iterator         unconnected_begin();
    member_iterator         unconnected_end();
    member_const_iterator       unconnected_begin() const;
    member_const_iterator       unconnected_end() const;


 /// ??? callback functions
    callback_type   set_build_callback(callback_type f);
    callback_type   set_freeze_callback(callback_type f);


  protected:
    virtual subped_id   allocate_subped(const string&, pedigree_base*);
    virtual member_id   allocate_member(const string&, SexCode, const void*);
    virtual family_id   allocate_family();
    virtual void        deallocate_subped(subped_id);
    virtual void        deallocate_member(member_id);
    virtual void        deallocate_family(family_id);

  private:
    pedinfo_type            my_info;
    callback_type           my_build_callback;
    callback_type           my_freeze_callback;

    pedigree(const pedigree&);
    void    operator =(const pedigree&);

    /// 
    /// Clears all entries from this pedigree.
    void clear();

    ///
    /// ???
    void flush_build_info();
};
}
}

#include "mped/pedigree.ipp"

#endif
