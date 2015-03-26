#ifndef _SPBASEITER_HPP
#define _SPBASEITER_HPP

//============================================================================
//  File:       spbaseiter.h
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Notes:      
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "mped/sptypes.h"

namespace SAGE {
namespace MPED {

/** \ingroup BaseIterators 
  * \brief Iterates across a set of family_base
  *
  * Iterates across a set of family_base.
  */
class family_base_iterator //: public std::random_access_iterator<family_base,fidx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class family_base_const_iterator;

  /** \name Typdefs
    * Typdef-ed local to this class
    */
  //@{

    typedef family_base*                        pointer;
    typedef family_base&                        reference;
    typedef fidx_vector::size_type              size_type;
    typedef fidx_vector::difference_type        difference_type;
    typedef fidx_vector::iterator               cursor_type;
    typedef family_base_iterator                this_type;

  //@}

  public:
  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    family_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    family_base_iterator(const cursor_type& iter);

    cursor_type my_iter;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across a family_base
  *
  * Const-iterates across a family_base.
  */
class family_base_const_iterator //: public std::random_access_iterator<family_base,fidx_vector::difference_type>
{
  public:
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typdefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const family_base*                  pointer;
    typedef const family_base&                  reference;
    typedef fidx_vector::size_type              size_type;
    typedef fidx_vector::difference_type        difference_type;
    typedef fidx_vector::iterator               cursor_type;
    typedef fidx_vector::const_iterator         const_cursor_type;
    typedef family_base_const_iterator          this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    family_base_const_iterator();
    family_base_const_iterator(const family_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    family_base_const_iterator(const cursor_type& iter);
    family_base_const_iterator(const const_cursor_type& iter);

    const_cursor_type   my_iter;
};

/** \ingroup BaseIterators
  * \brief Iterates across a mate_base
  *
  * Iterates across a mate_base.
  */
class mate_base_iterator //: public std::forward_iterator<mate_info_base,mate_multimap::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class mate_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef mate_info_base*                     pointer;
    typedef mate_info_base&                     reference;
    typedef mate_multimap::size_type            size_type;
    typedef mate_multimap::difference_type      difference_type;
    typedef mate_multimap::iterator             cursor_type;
    typedef mate_base_iterator                  this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    mate_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    mate_base_iterator(const cursor_type& iter);

    cursor_type my_iter;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across a mate_base
  *
  * Const-iterates across a mate_base.
  */
class mate_base_const_iterator //: public std::forward_iterator<mate_info_base,mate_multimap::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const mate_info_base*               pointer;
    typedef const mate_info_base&               reference;
    typedef mate_multimap::size_type            size_type;
    typedef mate_multimap::difference_type      difference_type;
    typedef mate_multimap::iterator             cursor_type;
    typedef mate_multimap::const_iterator       const_cursor_type;
    typedef mate_base_const_iterator            this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    mate_base_const_iterator();
    mate_base_const_iterator(const mate_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    mate_base_const_iterator(const cursor_type& iter);
    mate_base_const_iterator(const const_cursor_type& iter);

    const_cursor_type   my_iter;
};

/** \ingroup BaseIterators
  * \brief Iterates across a member_base
  *
  * Iterates across a member_base.
  */
class member_base_iterator //: public std::random_access_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class member_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base*                        pointer;
    typedef member_base&                        reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef midx_vector::iterator               cursor_type;
    typedef member_base_iterator                this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    member_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    member_base_iterator(const cursor_type& iter);

    cursor_type my_iter;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across a member_base
  *
  * Const-iterates across a member_base.
  */
class member_base_const_iterator //: public std::random_access_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base                         value_type;
    typedef const member_base*                  pointer;
    typedef const member_base&                  reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef midx_vector::iterator               cursor_type;
    typedef midx_vector::const_iterator         const_cursor_type;
    typedef member_base_const_iterator          this_type;
    typedef std::random_access_iterator_tag     iterator_category;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    member_base_const_iterator();
    member_base_const_iterator(const member_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    member_base_const_iterator(const cursor_type& iter);
    member_base_const_iterator(const const_cursor_type& iter);

    const_cursor_type   my_iter;
};

/** \ingroup BaseIterators
  * \brief Iterates across an offspring_base
  *
  * Iterates across an offspring_base.
  */
class offspring_base_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class offspring_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base*                        pointer;
    typedef member_base&                        reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef offspring_base_iterator             this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    offspring_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    offspring_base_iterator(member_id data);
    void        advance();

    member_id   my_data;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across an offspring_base
  *
  * Const-iterates across an offspring_base.
  */
class offspring_base_const_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member_base*                  pointer;
    typedef const member_base&                  reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef offspring_base_const_iterator       this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    offspring_base_const_iterator();
    offspring_base_const_iterator(const offspring_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    offspring_base_const_iterator(member_id data);
    void        advance();

    member_id   my_data;
};

/** \ingroup BaseIterators
  * \brief Iterates across a parent_base
  *
  * Iterates across a parent_base.
  */
class parent_base_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class parent_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base*                        pointer;
    typedef member_base&                        reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef parent_base_iterator                this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    parent_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    parent_base_iterator(family_id F);
    void        init();
    void        advance();

    family_id   my_family;
    member_id   my_data;
};


/** \ingroup BaseConstIterators
  * \brief Const-iterates across a parent_base
  *
  * Const-iterates across a parent_base.
  */
class parent_base_const_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member_base*                  pointer;
    typedef const member_base&                  reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef parent_base_const_iterator          this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    parent_base_const_iterator();
    parent_base_const_iterator(const parent_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    parent_base_const_iterator(family_id F);
    void        init();
    void        advance();

    family_id   my_family;
    member_id   my_data;
};


/** \ingroup BaseIterators
  * \brief Iterates across a pedigree_base
  *
  * Iterates across a pedigree_base.
  */
class pedigree_base_iterator //: public std::random_access_iterator<member_base,pidx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class multipedigree_base;
    friend class pedigree_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef pedigree_base*                      pointer;
    typedef pedigree_base&                      reference;
    typedef pidx_vector::size_type              size_type;
    typedef pidx_vector::difference_type        difference_type;
    typedef pidx_vector::iterator               cursor_type;
    typedef pedigree_base_iterator              this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{    

    pedigree_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    pedigree_base_iterator(const cursor_type& iter);

    cursor_type my_iter;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across a pedigree_base
  *
  * Const-iterates across a pedigree_base.
  */
class pedigree_base_const_iterator //: public std::random_access_iterator<member_base,pidx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class multipedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const pedigree_base*                pointer;
    typedef const pedigree_base&                reference;
    typedef pidx_vector::size_type              size_type;
    typedef pidx_vector::difference_type        difference_type;
    typedef pidx_vector::iterator               cursor_type;
    typedef pidx_vector::const_iterator         const_cursor_type;
    typedef pedigree_base_const_iterator        this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    pedigree_base_const_iterator();
    pedigree_base_const_iterator(const pedigree_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    pedigree_base_const_iterator(const cursor_type& iter);
    pedigree_base_const_iterator(const const_cursor_type& iter);

    const_cursor_type   my_iter;
};


/** \ingroup BaseIterators
  * \brief Iterates across a progeny_base
  *
  * Iterates across a progeny_base.
  */
class progeny_base_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class progeny_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base*                        pointer;
    typedef member_base&                        reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef mate_multimap::iterator             cursor_type;
    typedef progeny_base_iterator               this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    progeny_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    progeny_base_iterator(const cursor_type& iter, const cursor_type& limit);
    void        init();
    void        advance();

    member_id       my_data;
    cursor_type     my_iter;
    cursor_type     my_limit;
};

/** \ingroup BaseConstIterators
  * \brief Const-iterates across a progeny_base
  *
  * Const-iterates across a progeny_base.
  */
class progeny_base_const_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member_base*                  pointer;
    typedef const member_base&                  reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef mate_multimap::iterator             cursor_type;
    typedef mate_multimap::const_iterator       const_cursor_type;
    typedef progeny_base_const_iterator         this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    progeny_base_const_iterator();
    progeny_base_const_iterator(const progeny_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    progeny_base_const_iterator(const cursor_type& iter, const cursor_type& limit);
    progeny_base_const_iterator(const const_cursor_type& iter, const const_cursor_type& limit);
    void        init();
    void        advance();

    member_id           my_data;
    const_cursor_type   my_iter;
    const_cursor_type   my_limit;
};


/** \ingroup BaseIterators
  * \brief Iterates across a sibling_base
  *
  * Iterates across a sibling_base.
  */
class sibling_base_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class sibling_base_const_iterator;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member_base*                        pointer;
    typedef member_base&                        reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef sibling_base_iterator               this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    sibling_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    sibling_base_iterator(member_id self);
    sibling_base_iterator(member_id self, member_id data);
    void        advance();

    member_id   my_self;
    member_id   my_data;
};


/** \ingroup BaseConstIterators
  * \brief Const-iterates across a sibling_base
  *
  * Const-iterates across a sibling_base.
  */
class sibling_base_const_iterator //: public std::forward_iterator<member_base,midx_vector::difference_type>
{
  public:
    friend class member_base;
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member_base*                  pointer;
    typedef const member_base&                  reference;
    typedef midx_vector::size_type              size_type;
    typedef midx_vector::difference_type        difference_type;
    typedef sibling_base_const_iterator         this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    sibling_base_const_iterator();
    sibling_base_const_iterator(const sibling_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}

  private:
    sibling_base_const_iterator(member_id self);
    sibling_base_const_iterator(member_id self, member_id data);
    void        advance();

    member_id   my_self;
    member_id   my_data;
};


/** \ingroup BaseIterators
  * \brief Iterates across a subpedigree_base
  *
  * Iterates across a subpedigree_base.
  */
class subpedigree_base_iterator //: public std::random_access_iterator<subpedigree_base,sidx_vector::difference_type>
{
  public:
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;
    friend class subpedigree_base_const_iterator;

  /** \name Typedefs
    * Typdef-ed local to this class
    */
  //@{

    typedef subpedigree_base*                   pointer;
    typedef subpedigree_base&                   reference;
    typedef sidx_vector::size_type              size_type;
    typedef sidx_vector::difference_type        difference_type;
    typedef sidx_vector::iterator               cursor_type;
    typedef subpedigree_base_iterator           this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    subpedigree_base_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    subpedigree_base_iterator(const cursor_type& iter);

    cursor_type my_iter;
};


/** \ingroup BaseConstIterators
  * \brief Const-iterates across a subpedigree_base
  *
  * Const-iterates across a subpedigree_base.
  */
class subpedigree_base_const_iterator //: public std::random_access_iterator<subpedigree_base,sidx_vector::difference_type>
{
  public:
    friend class family_base;
    friend class subpedigree_base;
    friend class pedigree_base;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const subpedigree_base*             pointer;
    typedef const subpedigree_base&             reference;
    typedef sidx_vector::size_type              size_type;
    typedef sidx_vector::difference_type        difference_type;
    typedef sidx_vector::iterator               cursor_type;
    typedef sidx_vector::const_iterator         const_cursor_type;
    typedef subpedigree_base_const_iterator     this_type;

  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{

    subpedigree_base_const_iterator();
    subpedigree_base_const_iterator(const subpedigree_base_iterator& iter);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;
    bool        operator  <(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;
    reference   operator [](difference_type n) const;

    this_type&  operator +=(difference_type n);
    this_type&  operator -=(difference_type n);
    this_type&  operator ++();
    this_type&  operator --();
    this_type   operator ++(int);
    this_type   operator --(int);

  //@}

  private:
    subpedigree_base_const_iterator(const cursor_type& iter);
    subpedigree_base_const_iterator(const const_cursor_type& iter);

    const_cursor_type   my_iter;
};

} // End namespace MPED
} // End namespace SAGE

#include "mped/spbaseiter.ipp"

#endif  //- _SPBASEITER_HPP
