//****************************************************************************
//* File:      reltypename.cpp                                               *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     Maps conventional names onto relationship objects.            *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>
#include "pairs/reltype.h"
#include "pairs/reltypename.h"

namespace SAGE
{

std::string relationship_name(const RelationType& r)
{
  char val[25];

  if( r.distance1() + r.distance2() == 0 )
  {
    switch(r.bridge())
    {
      case SELF_BRIDGE:
        if( r.has_offspring() )
          return "invalid";
        else
          return "SELF";

      case NO_BRIDGE:
        if( r.has_offspring() )
          return "mate";
        return "none";

      default:
        return "invalid";
    }
  }

  if( r.distance1() == 1 && r.distance2() == 1 )
    switch(r.bridge())
    {
      case MZTWIN_BRIDGE:   return "MZTWIN";
      case DZTWIN_BRIDGE:   return "DZTWIN";
      case FULL_SIB_BRIDGE: return "SIB";
      case HALF_SIB_BRIDGE: return "H.SIB";
      default:              return "invalid";
    }

  int dl = std::min(r.distance1(), r.distance2() );
  int dh = std::max(r.distance1(), r.distance2() );

  if( dl == 0 && r.bridge() != NO_BRIDGE )
    return "invalid";

  if( dl == 0 && dh == 1 )
  {
    if( r.distance1() < r.distance2() )
      return "P-O";

    return "O-P";
  }

  if( dl == 0 && dh == 2 )
    return "GRANDP";

  if( dl == 0 && dh == 3 )
    return "G.GRANDP";

  if( dl == 0 )
  {
    sprintf(val, "G*%d GRANDP", (int)dh-2);
    return val;
  }

  if( dl == dh && (r.bridge() == SELF_BRIDGE || r.bridge() == NO_BRIDGE) )
    return "invalid";

  if( dl == 2 && dh == 2 && r.bridge() == HALF_SIB_BRIDGE)
    return "H.COUSIN";

  if( dl == 2 && dh == 2 )
    return "COUSIN";

  if( dl == dh && r.bridge() == HALF_SIB_BRIDGE )
  {
    sprintf(val, "H.COUSIN+%d", (int)dl-1);
    return val;
  }

  if( dl == dh )
  {
    sprintf(val, "%d COUSIN", (int)dl-1);
    return val;
  }

  if( dl == 1 && dh == 2 )
    switch(r.bridge())
    {
      case MZTWIN_BRIDGE:
      case DZTWIN_BRIDGE:
      case FULL_SIB_BRIDGE: return "AVUNC";
      case HALF_SIB_BRIDGE: return "H.AVUNC";
                   default: return "invalid";
    }

  if( dl == 1 && dh == 3 )
    switch(r.bridge())
    {
      case   MZTWIN_BRIDGE:
      case   DZTWIN_BRIDGE:
      case FULL_SIB_BRIDGE: return "G.AVUNC";
      case HALF_SIB_BRIDGE: return "G.H.AVUNC";
                   default: return "invalid";
    }

  if( dl == 1 )
    switch(r.bridge())
    {
      case   MZTWIN_BRIDGE:
      case   DZTWIN_BRIDGE:
      case FULL_SIB_BRIDGE:
                            sprintf(val, "G*%d AVUNC", (int)dh-2);
                            return val;
      case HALF_SIB_BRIDGE:
                            sprintf(val, "G*%d H.AVUNC", (int)dh-2);
                            return val;
                   default: return "invalid";
    }

  if( dl > 1 )
    switch(r.bridge())
    {
      case   MZTWIN_BRIDGE:
      case   DZTWIN_BRIDGE:
      case FULL_SIB_BRIDGE:
                             sprintf(val, "%d COUSIN+%d", (int)dl-1, (int)dh-(int)dl);
                             return val;

      case HALF_SIB_BRIDGE:
                             sprintf(val, "%d H.COUSIN+%d", (int)dl-1, (int)dh-(int)dl);
                             return val;

                   default:  return "invalid";
    }

  if( r.bridge() == NO_BRIDGE || r.bridge() == SELF_BRIDGE )
    return "invalid";

  return r.str();
}

static string gender_string(const pedigree_member_list& l)
{
  string s;

  for( pedigree_member_list::const_iterator j  = l.begin();
                                            j != l.end();
                                          ++j)
    switch( (*j)->get_detailed_sex() )
    {
      case MPED::SEX_MALE    : s += 'M'; break;
      case MPED::SEX_XMALE   : s += 'M'; break;
      case MPED::SEX_FEMALE  : s += 'F'; break;
      case MPED::SEX_XFEMALE : s += 'F'; break;
           default           : s += '?'; break;
    }

  return s;
}

static string gender_string_assumed(const pedigree_member_list& l)
{
  string s;

  for( pedigree_member_list::const_iterator j  = l.begin();
                                            j != l.end();
                                          ++j)
    switch( (*j)->get_effective_sex() )
    {
      case MPED::SEX_MALE    : s += 'M'; break;
      case MPED::SEX_FEMALE  : s += 'F'; break;
           default           : s += '?'; break;
    }

  return s;
}

std::pair<std::string, std::string> gender_strings(const CompleteRelationType& c,
                                                   bool show_inferred)
{
  std::pair<std::string, std::string> gs;

  if(show_inferred)
  {
    gs.first  = gender_string(c.lineage1);
    gs.second = gender_string(c.lineage2);
  }
  else
  {
    gs.first  = gender_string_assumed(c.lineage1);
    gs.second = gender_string_assumed(c.lineage2);
  }

  return gs;
}

} // end of namespace SAGE
