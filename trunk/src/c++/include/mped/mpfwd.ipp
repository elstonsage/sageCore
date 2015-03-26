//============================================================================
//  File:       mpfwd.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//
//============================================================================
// IMPLEMENTATION:  no_info
//============================================================================
//
namespace SAGE {
namespace MPED {

inline std::istream&
operator >>(std::istream& is, no_info&)
{
    return is;
}

inline std::ostream&
operator <<(std::ostream& os, const no_info&)
{
    return os;
}

/// Returns the effective sex, removing details of inferrence.
///
/// \param c the SexCode to test.
inline SexCode get_effective_sex(SexCode c)
{
  return (SexCode) (c & INFER_MASK);
}

/// Test for if a sex code is one labeled as male (SEX_MALE or SEX_XMALE)
///
/// \param c the SexCode to test.
inline bool is_male(SexCode c)
{
  return get_effective_sex(c) == SEX_MALE;
}

/// Test for if a sex code is one labeled as female (SEX_FEMALE or SEX_XFEMALE)
///
/// \param c the SexCode to test.
inline bool is_female(SexCode c)
{
  return get_effective_sex(c) == SEX_FEMALE;
}

/// Test for if a sex code is one labeled as unknown (SEX_MISSING or SEX_ARB)
///
/// \param c the SexCode to test.
inline bool is_sex_unknown(SexCode c)
{
  return get_effective_sex(c) == SEX_MISSING;
}

//============================================================================
// IMPLEMENTATION:  error_info
//============================================================================
//
inline
error_info::error_info() 
{}

inline
error_info::error_info(ErrorStatus s)
  : state(s) 
{}

inline
error_info::error_info(ErrorStatus s, const string& a)
  : state(s), name1(a) 
{}

inline
error_info::error_info(ErrorStatus s, const string& a, const string& b) 
  : state(s), name1(a), name2(b) 
{}

inline
error_info::error_info
(ErrorStatus s, const string& a, const string& b, const string& c) 
  : state(s), name1(a), name2(b), name3(c) 
{}

} // End namespace MPED
} // End namespace SAGE
