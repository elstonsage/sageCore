#ifndef PAIRS_STRINGBANK
#define PAIRS_STRINGBANK

//****************************************************************************
//* File:      relnamebank.h                                                 *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     Maps conventional names onto relationship objects.            *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>

namespace SAGE
{

extern const char* ordinal_names[];
extern const int   max_ordinal_names;
extern const char* ordinal_th_names[];
extern const int    max_ordinal_th_names;
extern const char* ordinal_tens_names[];
extern const int   max_ordinal_tens_names;
extern const char* ordinal_power_names[];
extern const int   max_ordinal_power_names;
extern const char* ordinal_st_names[];
extern const int   max_ordinal_st_names;
extern const char* relative_words[];
extern const int   max_relative_words;

} // end of namespace SAGE

#endif
