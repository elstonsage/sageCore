//============================================================================
//  File:       sptypes.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
// 
//============================================================================
// IMPLEMENTATION:  name_pair
//============================================================================
//
namespace SAGE {
namespace MPED {

//============================================================================
// IMPLEMENTATION:  mate_info_base
//============================================================================
//
inline 
mate_info_base::mate_info_base(member_id M, family_id F)
  : the_mate(M), the_family(F) 
{}

inline const mate_info_base::member_type&
mate_info_base::mate() const
{
    return *the_mate;
}

inline const mate_info_base::family_type&
mate_info_base::family() const
{
    return *the_family;
}

inline mate_info_base::member_type&
mate_info_base::mate()
{
    return *the_mate;
}

inline mate_info_base::family_type&
mate_info_base::family()
{
    return *the_family;
}


//============================================================================
// IMPLEMENTATION:  mate_info<GI,FI,SI,PI,MI>
//============================================================================
//
template <class GI, class FI, class SI, class PI, class MI> inline 
const typename mate_info<GI,FI,SI,PI,MI>::member_type&
mate_info<GI,FI,SI,PI,MI>::mate() const
{
    return reinterpret_cast<const member<GI,FI,SI,PI,MI>&>(base_type::mate());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
const typename mate_info<GI,FI,SI,PI,MI>::family_type&
mate_info<GI,FI,SI,PI,MI>::family() const
{
    return reinterpret_cast<const SAGE::MPED::family<GI,FI,SI,PI,MI>&>(base_type::family());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_info<GI,FI,SI,PI,MI>::member_type&
mate_info<GI,FI,SI,PI,MI>::mate()
{
    return reinterpret_cast<member<GI,FI,SI,PI,MI>&>(base_type::mate());
}

template <class GI, class FI, class SI, class PI, class MI> inline 
typename mate_info<GI,FI,SI,PI,MI>::family_type&
mate_info<GI,FI,SI,PI,MI>::family()
{
    return reinterpret_cast<SAGE::MPED::family<GI,FI,SI,PI,MI>&>(base_type::family());
}

} // End namespace MPED
} // End namespace SAGE
