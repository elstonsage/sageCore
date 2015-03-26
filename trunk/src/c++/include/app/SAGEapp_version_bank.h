#ifndef SAGE_VERSION_H
#define SAGE_VERSION_H

//****************************************************************************
//* File:      SAGEapp_version_bank.h                                        *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     Maps app index onto their app name & micro version numbers.   *
//*                                                                          *
//* Copyright (c) 2002 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <vector>
#include <string>

namespace SAGE {
namespace APP  {

/// Identifies a program with a unique ID number.  Applications are listed
/// alphabetically.
enum app_index_type
{
  APP_AGEON,
  APP_ASSOC,
  APP_DECIPHER,
  APP_FCOR,
  APP_FREQ,
  APP_GENIBD,
  APP_LODLINK,
  APP_LODPAL,
  APP_MARKERINFO,
  APP_MLOD,
  APP_PEDINFO,
  APP_RELPAL,
  APP_RELTEST,
  APP_SEGREG,
  APP_SIBPAL,
  APP_TDTEX,
  TEST_PROGRAM,
  APP_COUNT
};

/// String identifying the release version of SAGE.
extern const char* sage_release;

/** \brief Static class providing conversion functions for APP::app_index_type and micro version number.
 *
 */
class app_names
{
public:

    ///
    /// Populates the application name and micro version lists with data.
    static void construct();

    ///
    /// Returns a string version of the application name.
    /// \param t The application code to convert into a string.
    static const std::string& get_application_name(app_index_type t);

    ///
    /// Returns the micro version number of the application.
    /// \param t The applicatino code whose micro version will be returned.
    static const std::string& get_micro_version_number(app_index_type t);

protected:

    app_names();

    static std::vector<std::string> application_names;
    static std::vector<std::string> micro_version_numbers;
};

} // End namespace APP
} // End namespace SAGE

#endif
