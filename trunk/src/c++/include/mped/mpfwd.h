#ifndef _MPFWD_HPP
#define _MPFWD_HPP

//============================================================================
//  File:       mpfwd.h
//
//  Purpose:    This header file provides forward definitions and lays the
//              foundation for the new SAGE multi-pedigree class library.
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include <assert.h>
#include <iostream>
#include <string>
#include "globals/config.h"
#include "globals/data_types.h"

// This really should go away!

using std::string;

namespace SAGE {
namespace MPED {

/** \defgroup BaseStorageClasses Base storage classes
  *
  * The base storage classes keep track of pedigree structure only. Although you will generally
  * never instantiate them directly, you may very well write algorithms that process them.
  */

/** \defgroup DerivedStorageClasses Derived storage classes
  *
  * The derived storage classes extend the base storage classes by adding additional data storage
  * at each organizational level. The data storage is templatized, so that you can associate 
  * whatever data type you want with each structural level of data organization.
  */

/** \defgroup BaseIterators Base non-const iterators
  *
  * Non-const iterators for the base classes
  */

/** \defgroup DerivedIterators Derived non-const iterators
  *
  * Non-const iterators for the derived classes
  */

/** \defgroup BaseConstIterators Base const iterators
  *
  * Const iterators for the base classes
  */

/** \defgroup DerivedConstIterators Derived const iterators
  *
  * Const iterators for the derived classes
  */

/** \defgroup UtilityFunctions Utility functions
  *
  * Various utility functions (declared within mp_utilities) for processing multipedigrees.
  */


//- These are the base classes (or engines) and associated iterators that will 
//  drive the parametrized representations.
//
class family_base;
class member_base;
class subpedigree_base;
class pedigree_base;
class multipedigree_base;

class ancestor_base_iterator;
class descendant_base_iterator;
class family_base_iterator;
class mate_base_iterator;
class member_base_iterator;
class offspring_base_iterator;
class parent_base_iterator;
class relation_base_iterator;
class sibling_base_iterator;
class subpedigree_base_iterator;
class pedigree_base_iterator;

class ancestor_base_const_iterator;
class descendant_base_const_iterator;
class family_base_const_iterator;
class mate_base_const_iterator;
class member_base_const_iterator;
class offspring_base_const_iterator;
class parent_base_const_iterator;
class relation_base_const_iterator;
class sibling_base_const_iterator;
class subpedigree_base_const_iterator;
class pedigree_base_const_iterator;


//- These are the parameterized classes and associated iterators that provide
//  specific genetic information to our pedigree representation.
//
template <class GI, class FI, class SI, class PI, class MI> class member;
template <class GI, class FI, class SI, class PI, class MI> class subpedigree;
template <class GI, class FI, class SI, class PI, class MI> class pedigree;
template <class GI, class FI, class SI, class PI, class MI> class multipedigree;
template <class GI, class FI, class SI, class PI, class MI> class family;

template <class GI, class FI, class SI, class PI, class MI> class ancestor_iterator;
template <class GI, class FI, class SI, class PI, class MI> class descendant_iterator;
template <class GI, class FI, class SI, class PI, class MI> class family_iterator;
template <class GI, class FI, class SI, class PI, class MI> class mate_iterator;
template <class GI, class FI, class SI, class PI, class MI> class member_iterator;
template <class GI, class FI, class SI, class PI, class MI> class offspring_iterator;
template <class GI, class FI, class SI, class PI, class MI> class parent_iterator;
template <class GI, class FI, class SI, class PI, class MI> class relation_iterator;
template <class GI, class FI, class SI, class PI, class MI> class sibling_iterator;
template <class GI, class FI, class SI, class PI, class MI> class subpedigree_iterator;
template <class GI, class FI, class SI, class PI, class MI> class pedigree_iterator;

template <class GI, class FI, class SI, class PI, class MI> class ancestor_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class descendant_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class family_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class mate_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class member_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class offspring_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class parent_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class relation_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class sibling_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class subpedigree_const_iterator;
template <class GI, class FI, class SI, class PI, class MI> class pedigree_const_iterator;


/** \brief Empty struct for derived storage classes
  *
  * This is a special dummy class for creating a derived storage class with an uninformative
  * info object.
  *
  * \code
  * MPED::pedigree<no_info, no_info, no_info, no_info, no_info> m(...);
  * \endcode
  */
struct no_info {};

std::istream&   operator >>(std::istream& is, no_info&);
std::ostream&   operator <<(std::ostream& os, const no_info&);


/// Defines a sex code
///
/// Sex codes are designed such that the first bit (1) indicates an inferred
/// or inferrable status.  This allows us to mask off that bit when returning 
/// the effective sex, or include it for more detailed sex information.
enum SexCode 
{ 
    SEX_MALE     = 2,      /*!< Male (direct from source data)     */
    SEX_XMALE    = 3,      /*!< Male (inferred from source data)   */
    SEX_FEMALE   = 4,      /*!< Female (direct from source data)   */
    SEX_XFEMALE  = 5,      /*!< Female (inferred from source data) */
    SEX_MISSING  = 6,      /*!< Missing (source data has no info)  */
    SEX_ARB      = 7,      /*!< Sex may be assigned arbitrarily (assuming consistency) when inferring sexes */
};

inline SexCode get_effective_sex(SexCode c);

inline bool is_male        (SexCode c);
inline bool is_female      (SexCode c);
inline bool is_sex_unknown (SexCode c);

static const SexCode INFER_MASK = (SexCode) 6;

//- This struct provides information about errors that have occurred during
//  the construction of a pedigree.
//
struct error_info
{
    enum ErrorStatus
    {
        bad_sibship = 16, 
        bad_marriage, 
        bad_lineage,
        bad_gender,
        same_sex_marriage,
        gender_inferred,
        bad_marriage_loop,
        no_sex_parents
    };

    error_info();
    error_info(ErrorStatus s);
    error_info(ErrorStatus s, const string& a);
    error_info(ErrorStatus s, const string& a, const string& b);
    error_info(ErrorStatus s, const string& a, const string& b, const string& c);

    ErrorStatus  state;
    string       name1;
    string       name2;
    string       name3;
};

} // End namespace MPED
} // End namespace SAGE

#include "mped/mpfwd.ipp"

#endif
