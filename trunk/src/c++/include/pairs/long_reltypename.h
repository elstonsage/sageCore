#ifndef LONGRELATIONTYPENAME_H
#define LONGRELATIONTYPENAME_H

//****************************************************************************
//* File:      long_reltypename.h                                            *
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

std::string      long_relationship_name(const RelationType& r);
std::string      long_relationship_gender_name(const CompleteRelationType& c);

} // end of namespace SAGE

#endif
