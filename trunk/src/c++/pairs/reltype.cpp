//****************************************************************************
//* File:      reltype.cpp                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song & Kevin Jacobs                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     This source file defines a relationtype class to uniquely     *
//*            represent relationships between pairs of individuals          *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "pairs/reltype.h"

// ---------------------------------------------------------------------------
// Out-of-line implementation of RelationType
// ---------------------------------------------------------------------------

namespace SAGE
{

std::string RelationType::str() const
{
  char val[25];

  char p1 = ' ';
  char p2 = ' ';
  char b  = ' ';

  if     ( phase1() == PARENT1_PHASE )
    p1 = '-';
  else if( phase1() == PARENT2_PHASE )
    p1 = '+';

  if     ( phase2() == PARENT1_PHASE )
    p2 = '-';
  else if( phase2() == PARENT2_PHASE )
    p2 = '+';

  switch( bridge() )
  {
    case     SELF_BRIDGE: b = 'i'; break;
    case   MZTWIN_BRIDGE: b = 'm'; break;
    case   DZTWIN_BRIDGE: b = 'd'; break;
    case FULL_SIB_BRIDGE: b = 'f'; break;
    case HALF_SIB_BRIDGE: b = 'h'; break;
    default             : break;
  }

  if( has_offspring() )
    if(b == ' ')
      b = '!';
    else
      b = toupper(b);

  sprintf(val, "%c%d %c %d%c", p1, distance1(), b, distance2(), p2);

  return val;
}

// end of RelationType Implementation

} // end of namespace SAGE
