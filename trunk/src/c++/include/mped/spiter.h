#ifndef _SPITER_HPP
#define _SPITER_HPP

//============================================================================
//  File:       spiter.h
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
#include "mped/mpbase.h"

namespace SAGE {
namespace MPED {

/** \ingroup DerivedIterators
  * \brief Iterates across a set of family s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of family s
  */
template <class GI, class FI, class SI, class PI, class MI>
class family_iterator : public family_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef family<GI,FI,SI,PI,MI>*                     pointer;
    typedef family<GI,FI,SI,PI,MI>&                     reference;
    typedef family_iterator<GI,FI,SI,PI,MI>             this_type;
    typedef family_base_iterator                        base_type;
       
  //@}

  public:

  /** \name Standard iterator interface
    * Standard iterator interface
    */
  //@{

    family_iterator();

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
    family_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of family s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of family s
  */
template <class GI, class FI, class SI, class PI, class MI>
class family_const_iterator : public family_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const family<GI,FI,SI,PI,MI>*               pointer;
    typedef const family<GI,FI,SI,PI,MI>&               reference;
    typedef family_const_iterator<GI,FI,SI,PI,MI>       this_type;
    typedef family_base_const_iterator                  base_type;
       
  //@}
  
  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    family_const_iterator();
    family_const_iterator(const family_iterator<GI,FI,SI,PI,MI>& i);

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
    family_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of mate s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of mate s
  */
template <class GI, class FI, class SI, class PI, class MI>
class mate_iterator : public mate_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef mate_info<GI,FI,SI,PI,MI>*                  pointer;
    typedef mate_info<GI,FI,SI,PI,MI>&                  reference;
    typedef mate_iterator<GI,FI,SI,PI,MI>               this_type;
    typedef mate_base_iterator                          base_type;
       
  //@}
  
  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    mate_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    mate_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set a mate s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set a mate s
  */
template <class GI, class FI, class SI, class PI, class MI>
class mate_const_iterator : public mate_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const mate_info<GI,FI,SI,PI,MI>*            pointer;
    typedef const mate_info<GI,FI,SI,PI,MI>&            reference;
    typedef mate_const_iterator<GI,FI,SI,PI,MI>         this_type;
    typedef mate_base_const_iterator                    base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    mate_const_iterator();
    mate_const_iterator(const mate_iterator<GI,FI,SI,PI,MI>& i);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    mate_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of member s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of member s
  */
template <class GI, class FI, class SI, class PI, class MI>
class member_iterator : public member_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;
    friend class member_const_iterator<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member<GI,FI,SI,PI,MI>*                     pointer;
    typedef member<GI,FI,SI,PI,MI>&                     reference;
    typedef member_iterator<GI,FI,SI,PI,MI>             this_type;
    typedef member_base_iterator                        base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    member_iterator();

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
    member_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of member s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of member s
  */
template <class GI, class FI, class SI, class PI, class MI>
class member_const_iterator : public member_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member<GI,FI,SI,PI,MI>*               pointer;
    typedef const member<GI,FI,SI,PI,MI>&               reference;
    typedef member_const_iterator<GI,FI,SI,PI,MI>       this_type;
    typedef member_base_const_iterator                  base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    member_const_iterator();
    member_const_iterator(const member_iterator<GI,FI,SI,PI,MI>& i);

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
    member_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of offspring s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of offspring s
  */
template <class GI, class FI, class SI, class PI, class MI>
class offspring_iterator : public offspring_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member<GI,FI,SI,PI,MI>*                     pointer;
    typedef member<GI,FI,SI,PI,MI>&                     reference;
    typedef offspring_iterator<GI,FI,SI,PI,MI>          this_type;
    typedef offspring_base_iterator                     base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    offspring_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    offspring_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of offspring s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of offspring s
  */
template <class GI, class FI, class SI, class PI, class MI>
class offspring_const_iterator : public offspring_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member<GI,FI,SI,PI,MI>*               pointer;
    typedef const member<GI,FI,SI,PI,MI>&               reference;
    typedef offspring_const_iterator<GI,FI,SI,PI,MI>    this_type;
    typedef offspring_base_const_iterator               base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    offspring_const_iterator();
    offspring_const_iterator(const offspring_iterator<GI,FI,SI,PI,MI>& i);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    offspring_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of parent s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of parent s
  */
template <class GI, class FI, class SI, class PI, class MI>
class parent_iterator : public parent_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member<GI,FI,SI,PI,MI>*                     pointer;
    typedef member<GI,FI,SI,PI,MI>&                     reference;
    typedef parent_iterator<GI,FI,SI,PI,MI>             this_type;
    typedef parent_base_iterator                        base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    parent_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    parent_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of parent s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of parent s
  */
template <class GI, class FI, class SI, class PI, class MI>
class parent_const_iterator : public parent_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member<GI,FI,SI,PI,MI>*               pointer;
    typedef const member<GI,FI,SI,PI,MI>&               reference;
    typedef parent_const_iterator<GI,FI,SI,PI,MI>       this_type;
    typedef parent_base_const_iterator                  base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    parent_const_iterator();
    parent_const_iterator(const parent_iterator<GI,FI,SI,PI,MI>& i);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    parent_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of pedigree s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of pedigree s
  */
template <class GI, class FI, class SI, class PI, class MI>
class pedigree_iterator : public pedigree_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;
    friend class multipedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef pedigree<GI,FI,SI,PI,MI>*                   pointer;
    typedef pedigree<GI,FI,SI,PI,MI>&                   reference;
    typedef pedigree_iterator<GI,FI,SI,PI,MI>           this_type;
    typedef pedigree_base_iterator                      base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    pedigree_iterator();

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
    pedigree_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of pedigree s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of pedigree s
  */
template <class GI, class FI, class SI, class PI, class MI>
class pedigree_const_iterator : public pedigree_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;
    friend class multipedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const pedigree<GI,FI,SI,PI,MI>*             pointer;
    typedef const pedigree<GI,FI,SI,PI,MI>&             reference;
    typedef pedigree_const_iterator<GI,FI,SI,PI,MI>     this_type;
    typedef pedigree_base_const_iterator                base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    pedigree_const_iterator();
    pedigree_const_iterator(const pedigree_iterator<GI,FI,SI,PI,MI>& i);

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
    pedigree_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of progeny s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of progeny s
  */
template <class GI, class FI, class SI, class PI, class MI>
class progeny_iterator : public progeny_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member<GI,FI,SI,PI,MI>*                     pointer;
    typedef member<GI,FI,SI,PI,MI>&                     reference;
    typedef progeny_iterator<GI,FI,SI,PI,MI>            this_type;
    typedef progeny_base_iterator                       base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    progeny_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    progeny_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of progeny s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of progeny s
  */
template <class GI, class FI, class SI, class PI, class MI>
class progeny_const_iterator : public progeny_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member<GI,FI,SI,PI,MI>*               pointer;
    typedef const member<GI,FI,SI,PI,MI>&               reference;
    typedef progeny_const_iterator<GI,FI,SI,PI,MI>      this_type;
    typedef progeny_base_const_iterator                 base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    progeny_const_iterator();
    progeny_const_iterator(const progeny_iterator<GI,FI,SI,PI,MI>& i);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    progeny_const_iterator(const base_type& i);
};


/** \ingroup DerivedIterators
  * \brief Iterates across a set of sibling s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of sibling s
  */
template <class GI, class FI, class SI, class PI, class MI>
class sibling_iterator : public sibling_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef member<GI,FI,SI,PI,MI>*                     pointer;
    typedef member<GI,FI,SI,PI,MI>&                     reference;
    typedef sibling_iterator<GI,FI,SI,PI,MI>            this_type;
    typedef sibling_base_iterator                       base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    sibling_iterator();

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    sibling_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of sibling s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of sibling s
  */
template <class GI, class FI, class SI, class PI, class MI>
class sibling_const_iterator : public sibling_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const member<GI,FI,SI,PI,MI>*               pointer;
    typedef const member<GI,FI,SI,PI,MI>&               reference;
    typedef sibling_const_iterator<GI,FI,SI,PI,MI>      this_type;
    typedef sibling_base_const_iterator                 base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    sibling_const_iterator();
    sibling_const_iterator(const sibling_iterator<GI,FI,SI,PI,MI>& i);

    bool        operator ==(const this_type& c) const;
    bool        operator !=(const this_type& c) const;

    pointer     operator ->() const;
    reference   operator  *() const;

    this_type&  operator ++();
    this_type   operator ++(int);

  //@}
  

  private:
    sibling_const_iterator(const base_type& i);
};



/** \ingroup DerivedIterators
  * \brief Iterates across a set of subpedigree s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Iterates across a set of subpedigree s
  */
template <class GI, class FI, class SI, class PI, class MI>
class subpedigree_iterator : public subpedigree_base_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef subpedigree<GI,FI,SI,PI,MI>*                pointer;
    typedef subpedigree<GI,FI,SI,PI,MI>&                reference;
    typedef subpedigree_iterator<GI,FI,SI,PI,MI>        this_type;
    typedef subpedigree_base_iterator                   base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    subpedigree_iterator();

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
    subpedigree_iterator(const base_type& i);
};


/** \ingroup DerivedConstIterators
  * \brief Const-iterates across a set of subpedigree s
  *
  * \par Template parameters
  *
  * \c GI Info object associated with individuals (SAGE::MPED::member)
  * \c FI Info object associated with families (SAGE::MPED::family)
  * \c SI Info object associated with subpedigrees (SAGE::MPED::subpedigree)
  * \c PI Info object associated with pedigrees (SAGE::MPED::pedigree)
  * \c MI Info object associated with multipedigrees (SAGE::MPED::multipedigree)
  * 
  * Const-iterates across a set of subpedigree s
  */
template <class GI, class FI, class SI, class PI, class MI>
class subpedigree_const_iterator : public subpedigree_base_const_iterator
{
  public:
    friend class member<GI,FI,SI,PI,MI>;
    friend class family<GI,FI,SI,PI,MI>;
    friend class subpedigree<GI,FI,SI,PI,MI>;
    friend class pedigree<GI,FI,SI,PI,MI>;

  /** \name Typedefs
    * Typedef-ed local to this class
    */
  //@{

    typedef const subpedigree<GI,FI,SI,PI,MI>*          pointer;
    typedef const subpedigree<GI,FI,SI,PI,MI>&          reference;
    typedef subpedigree_const_iterator<GI,FI,SI,PI,MI>  this_type;
    typedef subpedigree_base_const_iterator             base_type;
       
  //@}
  

  public:
    
  /** \name Standard iterator interface
    * Standard iterator interface      
    */
  //@{
    
    subpedigree_const_iterator();
    subpedigree_const_iterator(const subpedigree_iterator<GI,FI,SI,PI,MI>& i);

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
    subpedigree_const_iterator(const base_type& i);
};

} // End namespace MPED
} // End namespace SAGE

#include "mped/spiter.ipp"

#endif  //- _SPITER_HPP
