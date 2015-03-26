  //============================================================================
//  File:       spbase1.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <mpfwd.h>
    #pragma hdrstop
#endif

#include <algorithm>
#include "mped/spbase.h"

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: mate_info
//============================================================================
//
const string&
mate_info_base::name() const
{
    return the_mate->name();
}


} // End namespace MPED
} // End namespace SAGE
