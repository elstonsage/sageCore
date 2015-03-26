//****************************************************************************
//* File:      long_reltypename.cpp                                          *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                    1.1  yjs Added gender_list()                          *
//*                                 & relationship_number()        Mar 01 01 *
//*                                                                          *
//* Notes:     Maps conventional names onto relationship objects.            *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include <string>
#include "pairs/reltype.h"
#include "pairs/stringbank.h"
#include "pairs/reltypename.h"

namespace SAGE
{

static void get_gender        (string& val, char gender);
static void get_phase         (string& val, char gender);
static void get_parent        (string& val, char gender);
static void get_offspring     (string& val, char gender);
static void get_sibling       (string& val, char gender);
static void get_avuncular     (string& val, char gender);
static void get_nephewship    (string& val, char gender);
static void get_great         (string& val);
static void get_half          (string& val);
static void get_half_gender   (string& s1,    string& s2, char gender);
static void get_ancestor_phase(string& phase, string& gr, const string& gs, size_t dh);
static void get_ancestor_phase(string& phase,             const string& gs, size_t dh);
static void get_ancestor      (string& phase, string& gr, const string& gs, size_t dh);

std::string long_relationship_name(const RelationType& r)
{
  string val;

  // Self or Spouse relationship
  if( r.distance1() + r.distance2() == 0 )
    switch(r.bridge())
    {
      case SELF_BRIDGE:
        if( r.has_offspring() )
          return relative_words[29]; // invalid

        // Self
        return relative_words[6]; // self

      case NO_BRIDGE:
        // Spouse
        if( r.has_offspring() )
        {
          val  = relative_words[8]; // mother
          val += ':';
          val += relative_words[7]; // father
          return val;
        }

        // None
        return relative_words[28]; // none

      default:
        return relative_words[29]; // invalid
    }

  // Sibling relationship
  if( r.distance1() == 1 && r.distance2() == 1 )
  {
    if( r.bridge() == FULL_SIB_BRIDGE )
      return relative_words[30]; // sibling

    if( r.bridge() == MZTWIN_BRIDGE )
      return relative_words[15]; // mztwin

    if( r.bridge() == DZTWIN_BRIDGE )
      return relative_words[16]; // dztwin

    if( r.bridge() == HALF_SIB_BRIDGE )
    { 
      get_half(val);

      val += relative_words[30]; // sibling
      return val;
    }

    return relative_words[29]; // invalid
  }

  size_t dl = std::min(r.distance1(), r.distance2() );
  size_t dh = std::max(r.distance1(), r.distance2() );

  if( dl == 0 && r.bridge() != NO_BRIDGE )
    return relative_words[29]; //invalid

  // Parent-Offspring relationship
  if( dl == 0 && dh == 1 )
  {
    val  = relative_words[9]; // parent
    val += ':';
    val += relative_words[12]; // offspring
    return val;
  }

  // Grandparent-Grandchild relationship
  if( dl == 0 && dh == 2 )
  {
    val  = relative_words[18]; // grand
    val += relative_words[9]; // parent
    return val;
  }

  // Great* Grandparent-Great* Grandchild relationship
  if( dl == 0 )
  {
    for( size_t d = 0; d < dh-2; ++d )
      get_great(val);

    val += relative_words[18]; // grand
    val += relative_words[9]; // parent
    return val;
  }

  if( dl == dh && (r.bridge() == SELF_BRIDGE || r.bridge() == NO_BRIDGE) )
    return relative_words[29]; // invalid

  // Cousin relationship
  if( dl == 2 && dh == 2 )
  {
    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half(val);

    val += relative_words[26];    //COUSIN
    return val;
  }

  // *th Cousin relationship
  if( dl == dh )
  {
    val += ordinal_th_names[dl-1];
    val += '-';

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half(val);

    val += relative_words[26];     //COUSIN
    return val;
  }

  // Avuncular relationship
  if( dl == 1 && dh == 2 )
  {
    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half(val);

    val += relative_words[22];     //Avuncular
    return val;
  }

  // Great* Avucular relationship
  if( dl == 1 )
  {
    for( size_t d = 0; d < dh-2; ++d )
      get_great(val);

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half(val);

    val += relative_words[22];     //Avuncular
    return val;
  }

  if( dl > 1 )
  {
    val += ordinal_th_names[dl-1];
    val += '-';

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half(val);

    val += relative_words[26];     //Cousin
    val += '-';
    if( dh - dl  < 4 ) 
      val += ordinal_st_names[dh-dl];
    else
    {
      val += ordinal_names[dh-dl];
      val += '-';
      val += ordinal_st_names[4];
    }
    val += '-';
    val += relative_words[27];     //Removed
 
    return val;
  }

  if( r.bridge() == NO_BRIDGE || r.bridge() == SELF_BRIDGE )
    return relative_words[29];

  return r.str();
}

std::string long_relationship_gender_name(const CompleteRelationType& c)
{
  string val;

  const RelationType& r = c.relationship;
  std::pair<std::string, std::string> gs = gender_strings(c, false);

  // Self or Spouse relationship
  if( r.distance1() + r.distance2() == 0 )
    switch(r.bridge())
    {
      case SELF_BRIDGE:
        if( r.has_offspring() )
          return relative_words[29]; // invalid

        // Male self
        if( gs.first[0] == 'M' || gs.second[0] == 'M' )
        {
          val  = relative_words[0]; // male
          val += '-';
          val += relative_words[6]; // self
          return val;
        }
        // Female self
        if( gs.first[0] == 'F' || gs.second[0] == 'F' )
        {
          val  = relative_words[1]; // female
          val += '-';
          val += relative_words[6]; // self
          return val;
        }
        // Unknown_gender self
        val  = relative_words[2]; // unknown_gender
        val += '-';
        val += relative_words[6]; // self
        return val;

      case NO_BRIDGE:
        // Spouse
        if( r.has_offspring() )
        {
          if( gs.first[0] == 'M' || gs.first[0] == 'F' )
          {
            val  = relative_words[8]; // mother
            val += ':';
            val += relative_words[7]; // father
          }
          else
          {
            val  = relative_words[9]; // parent
            val += ':';
            val += relative_words[9]; // parent
          }

          return val;
        }

        // None
        return relative_words[28]; // none

      default:
        return relative_words[29]; // invalid
    }

  // Sibling relationship
  if( r.distance1() == 1 && r.distance2() == 1 )
  {
    string s1;
    string s2;

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half_gender(s1, s2, gs.first[1]);

    get_sibling(s1, gs.first[0]);
    get_sibling(s2, gs.second[0]);

    val = s2 + ":" + s1;

    if( r.bridge() == MZTWIN_BRIDGE )
    {
      val += '-';
      val += relative_words[15];
      return val;
    }    
    if( r.bridge() == DZTWIN_BRIDGE )
    {
      val += '-';
      val += relative_words[16];
      return val;
    }

    return val;
  }

  size_t dl = std::min(r.distance1(), r.distance2() );
  size_t dh = std::max(r.distance1(), r.distance2() );

  if( dl == 0 && r.bridge() != NO_BRIDGE )
    return relative_words[29];

  // Parent-Offspring relationship
  if( dl == 0 && dh == 1 )
  {
    string s1;
    string s2;

    get_offspring(s1, gs.first[0]);
    get_parent   (s2, gs.first[1]);

    val = s2 + ":" + s1;
    return val;
  }

  // Grandparent-Grandchild relationship
  if( dl == 0 && dh == 2 )
  {
    string s1 = relative_words[18]; //grand
    string s2 = relative_words[18];

    get_offspring(s1, gs.first[0]);
    get_parent   (s2, gs.first[2]);

    string phase = relative_words[32]; //through
    phase += '-';
    get_parent(phase, gs.first[1]);

    val = s2 + '-' + phase + ":" + s1;
    return val;
  }

  // Great* Grandparent-Great* Grandchild relationship
  if( dl == 0 )
  {
    string gr;
    get_great(gr);

    string phase;
    get_ancestor_phase(phase, gr, gs.first, dh);

    string s1 = gr + relative_words[18]; //grand
    string s2 = gr + relative_words[18]; //grand

    get_offspring(s1, gs.first[0]);
    get_parent   (s2, gs.first[dh]);

    val = s2 + '-' + phase + ":" + s1;
    return val;
  }

  if( dl == dh && (r.bridge() == SELF_BRIDGE || r.bridge() == NO_BRIDGE) )
    return relative_words[29];

  // Cousin relationship
  if( dl == 2 && dh == 2 )
  {
    string s1;
    string s2;

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half_gender(s1, s2, gs.first[2]);

    get_gender(s1, gs.first[0]);
    get_gender(s2, gs.second[0]);

    s1 += '-';
    s2 += '-';

    s1 += relative_words[26]; //COUSIN
    s2 += relative_words[26]; //COUSIN

    string phase1 = relative_words[32];// through
    string phase2 = relative_words[32];// through

    phase1 += '-';
    phase2 += '-';

    get_parent(phase1, gs.first[1]);
    get_parent(phase2, gs.second[1]);

    val = s2 + '-' + phase2 + ":" + s1 + '-' + phase1;
    return val;
  }

  // *th Cousin relationship
  if( dl == dh )
  {
    string s1;
    string s2;

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half_gender(s1, s2, gs.first[dh]);

    get_gender(s1, gs.first[0]);
    get_gender(s2, gs.second[0]);

    s1 += '-';
    s2 += '-';

    s1 += ordinal_th_names[dl-1];
    s2 += ordinal_th_names[dl-1];

    s1 += '-';
    s2 += '-';

    s1 += relative_words[26]; //COUSIN
    s2 += relative_words[26]; //COUSIN

    string phase1;
    string phase2;

    get_ancestor_phase(phase1, gs.first,  dh);
    get_ancestor_phase(phase2, gs.second, dh);

    val = s2 + '-' + phase2 + ":" + s1 + '-' + phase1;
    return val;
  }

  // Avuncular relationship
  if( dl == 1 && dh == 2 )
  {
    string s1;
    string s2;

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half_gender(s1, s2, gs.first[2]);

    get_nephewship(s1, gs.first[0]);
    get_avuncular(s2, gs.second[0]);

    string phase = relative_words[32]; //through
    phase += '-';
    get_parent(phase, gs.first[1]);

    val = s2 + '-' + phase + ":" + s1;
    return val;
  }

  // Great* Avucular relationship
  if( dl == 1 )
  {
    string s1;
    string s2;

    if( r.bridge() == HALF_SIB_BRIDGE )
      get_half_gender(s1, s2, gs.first[dh]);

    string gr;
    get_great(gr);

    string phase;
    get_ancestor_phase(phase, gr, gs.first, dh);

    s1 += gr;
    s2 += gr;

    get_nephewship(s1, gs.first[0]);
    get_avuncular (s2, gs.second[0]);

    val = s2 + '-' + phase + ":" + s1;
    return val;
  }

  // *th cousin *times removed
  if( dl > 1 )
  {
    string removed_sex;
    pedigree_member_list first  = c.lineage1;
    CompleteRelationType rel;

    for( size_t i = 0; i < dh-dl; ++i )
    {
      first.pop_front();
      removed_sex += gs.first[i];
    }

    rel.relationship.set_bridge(r.bridge());
    rel.relationship.set_distance1(dl);
    rel.relationship.set_distance2(dl);
    rel.relationship.set_has_offspring(r.has_offspring());
    rel.lineage1 = first;
    rel.lineage2 = c.lineage2;

    string cousin = long_relationship_gender_name(rel);

    string removed = ",";
    if( removed_sex.size() == 1 )
    {
      removed += relative_words[34]; //removed
      removed += '-';

      get_offspring(removed, removed_sex[0]);
    }
    else if( removed_sex.size() > 1 )
    {
      removed += relative_words[34]; //removed
      removed += '-';

      size_t size = removed_sex.size();

      get_offspring(removed, removed_sex[size-1]);
      removed += relative_words[33]; //",s"
      removed += '-';

      for( size_t d = size-2; d > 0; --d )
      {
        get_offspring(removed, removed_sex[d]);

        removed += relative_words[33]; //",s"
        removed += '-';
      }
      get_offspring(removed, removed_sex[0]);
    }
    else
      removed = "invalid";

    val  = cousin + removed;
    return val;
  }

  if( r.bridge() == NO_BRIDGE || r.bridge() == SELF_BRIDGE )
  {
    val = relative_words[29];
    return val;
  }

  return r.str();
}

static void get_gender(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[0]; //male
  else if( gender == 'F' )
    val += relative_words[1]; //female
  else
    val += relative_words[2]; //unknown_gender
}

static void get_phase(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[3]; //paternal
  else if( gender == 'F' )
    val += relative_words[4]; //maternal
  else
    val += relative_words[5]; //unknown_phase
}


static void get_parent(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[7]; //father
  else if( gender == 'F' )
    val += relative_words[8]; //mother
  else
    val += relative_words[9]; //parent
}

static void get_offspring(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[10]; //son
  else if( gender == 'F' )
    val += relative_words[11]; //daughter
  else
    val += relative_words[12]; //offspring
}

static void get_sibling(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[13]; //brother
  else if( gender == 'F' )
    val += relative_words[14]; //sister
  else
    val += relative_words[30]; //sibling
}

static void get_avuncular(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[20]; //uncle
  else if( gender == 'F' )
    val += relative_words[21]; //aunt
  else
    val += relative_words[22]; //avuncular
}

static void get_nephewship(string& val, char gender)
{
  if( gender == 'M' )
    val += relative_words[23]; //nephew
  else if( gender == 'F' )
    val += relative_words[24]; //niece
  else
    val += relative_words[25]; //nephewship
}

static void get_great(string& gr)
{
  gr += relative_words[19]; //great
  gr += '-';
}

static void get_half(string& val)
{
  val += relative_words[17]; //Half
  val += "-";
}

static void get_half_gender(string& s1, string& s2, char gender)
{
  get_phase(s2, gender);

  s2 += '-';
  
  get_half(s1);
  get_half(s2);
}

static void get_ancestor_phase(string& phase, const string& gs, size_t dh)
{
  string gr;
  get_ancestor_phase(phase, gr, gs, dh);
}

static void get_ancestor_phase(string& phase, string& gr, const string& gs, size_t dh)
{
  phase += relative_words[32]; //through
  phase += '-';

  get_parent(phase, gs[1]);    //parent phase

  phase += relative_words[33]; //'s
  phase += '-';

  get_ancestor(phase, gr, gs, dh);

  get_parent(phase, gs[dh-1]); //farthest ancestor phase
}

static void get_ancestor(string& phase, string& gr, const string& gs, size_t dh)
{
  for( size_t d = 2; d < dh-1; ++d )
  {
    get_great(gr);
    get_parent(phase, gs[d]);

    phase += relative_words[33]; //'s
    phase += '-';
  }
}

std::string
gender_lists(const CompleteRelationType& c)
{
  std::pair<string, string> gs = gender_strings(c);

  reverse(gs.first.begin(),  gs.first.end() );

  string gender_list;
  if( gs.first.size() && gs.second.size() )
  {
    if( c.relationship.bridge() == SAGE::HALF_SIB_BRIDGE )
    {
      char phase = gs.first[0];
      gs.first.erase(gs.first.begin());
      string::iterator last = gs.second.end();
      --last;
      gs.second.erase(last);
      gender_list = gs.second + "," + phase + "," + gs.first;
    }
    else
      gender_list = gs.second + "," + gs.first;
  }
  else
    gender_list = gs.second + gs.first;

  return gender_list;
}

std::string
relationship_number(const RelationType& r)
{
  string rel_num;

  if( !r.distance1() && !r.distance2() )
  {
    if( r.bridge() == SAGE::SELF_BRIDGE )
      return "0";
    else if( r.has_offspring() )
      return "1";
    else
      return "?";
  }
  else
  {
    rel_num  = long2str(r.distance1());
    rel_num += long2str(r.distance2());
  }

  if( r.bridge() == SAGE::HALF_SIB_BRIDGE )
    rel_num += "h";

  return rel_num;
}

} // end of namespace SAGE
