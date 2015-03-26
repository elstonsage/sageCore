#ifndef DATA_TYPES_H
#define DATA_TYPES_H
//============================================================================
//  File:       data_types.h
//
//  Purpose:    This header file provides typedefs of globally used data types
//              within the SAGE code.  Originally part of several config files 
//              in several directories and combined here for common access
//
//  Author:     Geoff Wedig
//
//  History:    Version 0.01  Initial creation  Aug 2004
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//============================================================================
//

#include <string>

namespace SAGE
{
  
// These are common-sense typedefs to cut down on unnecessary typing.
//
// They are only necessary on some platforms, so we define them only as needed.

#if defined(MINGW)

/// Basic typedef for the unsigned integer
typedef unsigned int    uint;

/// Basic typedef for an unsigned short
typedef unsigned short  ushort;
#endif

/// Basic typedef for an unsigned char
typedef unsigned char   uchar;

#if defined(__MACOS__)

typedef unsigned int    uint;
typedef unsigned short  ushort;
#endif

namespace GLOBALS {

#ifdef __WIN32__
  const std::string path_delimiter = "\\";
#else
  const std::string path_delimiter = "/";
#endif

} // End namepsace GLOBALS
} // End namesapce SAGE

#endif  

