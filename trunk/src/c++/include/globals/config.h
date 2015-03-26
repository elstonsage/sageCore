#ifndef CONFIG_H
#define CONFIG_H
//============================================================================
//  File:       config.h
//
//  Purpose:    This header file provides forward definitions and lays the
//              foundation for the SAGE configuration.  Originally part of
//              several config files in several directories and combined here
//              for common access
//
//  Author:     Geoff Wedig
//
//  History:    Version 0.01  Initial creation  Aug 2004
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//============================================================================


// This file is primarily for system specific or global configuration options.
// These global options are used throughout the SAGE code.

// Sun's computation need to be specifically made ieee floating point
// compliant.

#if defined(sun)
#  include "ieeefp.h"
#endif

// Windows needs isnan to be made global

#if defined(WIN32)
#include <cmath>

using std::isnan;
#endif

#endif

