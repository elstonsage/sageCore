//
//  Name Manager 1.0 -- Manager for string id mappings
//  
//  Purpose: Allows lists of strings and id mappings to be stored and
//           updated.
//  
//  Limitations:  Current implementation is very simplistic.
//
//  Authors: Kevin Jacobs (jacobs@darwin.cwru.edu)
//           Geoff Wedig  (wedig@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation (Attributes)	Mar  5 1996
//                   gcw  Made Generic				Jul 11 1996  
//
//  Copyright (c) 1996  R.C. Elston
//

#include "LSF/NameMgr.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

id_handle NameManager::query(const string &s)
{
  NameMap::iterator i;

  if(!s.length()) return 0;

  for(i=name_map.begin(); i != name_map.end() && 
                         !((*i).second.first == s); i++);
  if( i != name_map.end() ) return (*i).first;
  return (id_handle) -1;
}

const string &NameManager::query(id_handle id)
{
  NameMap::iterator i = name_map.find(id);
  if( i != name_map.end() )
    return (*i).second.first;
  return nilString;
}

void NameManager::release(id_handle id)
{
  NameMap::iterator i = name_map.find(id);
  if ( i != name_map.end() )
    if ( ( --(*i).second.second ) <= 0 )
      name_map.erase(i);
}

id_handle NameManager::set(const string &s, bool b)
{
  int c=set(s,count, b);
  if( c == (int)count ) 
    count++;
  return c;
}

id_handle NameManager::set(const string &s, id_handle c, bool b)
{
  if (!s.length()) return 0;

  NameMap::iterator i = name_map.begin();

  for(; i != name_map.end() && !((*i).second.first == s); i++);

  if( i != name_map.end() )
  {
    (*i).second.second++;
    return (*i).first;
  }
  if (!b) return (id_handle) -1;

  name_map[c] = npair(s,1);
  return c;
}

int NameManager::ref(id_handle id)
{
  NameMap::iterator i = name_map.find(id);
  if ( i != name_map.end() )
    return (*i).second.second;
  return 0;
}
