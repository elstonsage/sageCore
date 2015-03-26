#ifndef NAME_MGR_H
#define NAME_MGR_H
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


#include <string>
#include <utility>
#include <map>
#include <iterator>
#include <vector>
#include <stdio.h>

using namespace std;

typedef unsigned int  id_handle;

static std::string nilString;

class NameManager
{
private:
  typedef pair<string, int>                       npair;
  typedef std::map<id_handle, npair, less<id_handle> > NameMap;
public:
  NameManager(id_handle c = 0x10000 ) { count = c; }
  virtual ~NameManager() { }

  virtual id_handle query(const string &s);
  virtual const string &query(id_handle id);
  virtual id_handle set(const string &s, bool b);
  virtual id_handle set(const string &s) { return set(s, false); }
  virtual id_handle add(const string &s) { return set(s, true); }
  void release(id_handle);
  int ref(id_handle);

protected:
  id_handle set(const string &s, id_handle i, bool b);
  NameMap name_map;
  id_handle count;
};

#endif
