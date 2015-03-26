#ifndef RELATIONTYPENAME_H
#define RELATIONTYPENAME_H

//****************************************************************************
//* File:      reltypename.h                                                 *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     Maps conventional names to relationship objects.              *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>
#include "pairs/reltype.h"

namespace SAGE
{

std::string                        relationship_name(const RelationType& r);
std::pair<std::string,
          std::string>             gender_strings(const CompleteRelationType& c,
                                                  bool show_inferred = true);

std::string                        long_relationship_name(const RelationType& r);
std::string                        long_relationship_gender_name(const CompleteRelationType& c);

std::string                        gender_lists(const CompleteRelationType& c);
std::string                        relationship_number(const RelationType& r);

} // end of namespace SAGE

#endif
