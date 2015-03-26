#ifndef NUMERICS_PRINT_UTIL_H
#define NUMERICS_PRINT_UTIL_H

//**************************************************************************
// File:           print_util.h
//
// Author:
//
// History:  Version 0.0 Initial implementation
//           Modified by:    Alexandra Borosova, Jul. 2004
//           Modification 1: "Fixing problem with asterisk."
//           Last change:    Alexandra Borosova, 07/08/04
//
// Copyright (c) 	2001 R.C. Elston
//   All Rights Reserved
//****************************************************************************

#include <string>

namespace SAGE
{

std::string fp(double d, size_t w, size_t p=4, char invalid='-');
std::string pval(double p, size_t w, int prec=-1, int max_stars=2);

//* Modification 1: "Fixing problem with asterisk."
//	Last change: 	Alexandra Borosova, 07/08/04
std::string fp_scientific(double d, size_t w, size_t p=4, char invalid='-');
std::string pval_scientific(double p, size_t w, int prec=-1, int max_stars=2);

}

#endif
