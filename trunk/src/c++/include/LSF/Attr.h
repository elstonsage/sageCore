#ifndef __ATTR_LIST_H
#define __ATTR_LIST_H
//
//  ATTRIBUTE LISTS 1.0 -- Low level attribute list object
//
//  Purpose: Implements a list of arbitrary attributes of fixed
//           symbol and value types.
//
//  Limitations:  Current implementation is very simplistic.
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation	Mar 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//


#include <string>
#include <map>
#include <iterator>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <limits>
#include "LSF/NameMgr.h"
#include "LSF/parse_ops.h"

using namespace std;

typedef id_handle attr_id;

class AttrManager : public NameManager
{
public:
  AttrManager(id_handle c = 0x100 ): NameManager(c)
  { (void) NameManager::set("VALUE",0, true); }

  virtual id_handle query(const string &s)
          { return NameManager::query(toUpper(s)); }
  virtual const string &query(id_handle id)
          { return NameManager::query(id); }
  virtual id_handle set(const string &s, bool b)
          { return NameManager::set( toUpper(s), b); }
  virtual id_handle set(const string &s)
          { return NameManager::set( toUpper(s), false); }
  virtual id_handle add(const string &s)
          { return NameManager::set( toUpper(s), true); }
};

extern AttrManager AttrNameMgr;

class AttrVal
{
public:
  enum attr_t { Attr_None, Attr_Integer, Attr_Real, Attr_String, Attr_Ptr };

private:
  attr_t attr_type;
  union attr_val_t
  {
    int   i;
    double f;
    string *s;
    const void *p;
  } attr_val;

public:
  AttrVal()                { clear(); attr_type  = Attr_None; clear(); }
  AttrVal(const     int i) { clear(); attr_type  = Attr_Integer; attr_val.i = i; }
  AttrVal(const  double f) { clear(); attr_type  = Attr_Real;    attr_val.f = f; }
  AttrVal(const string &s) { clear(); attr_type  = Attr_String;
                             attr_val.s = new string(s); }
  AttrVal(const   char *s) { clear(); attr_type  = Attr_String;
                             if(s) attr_val.s = new string(s);
                             else  attr_val.s = new string();  }
  AttrVal(const   char  s) { clear(); attr_type  = Attr_String;
                             attr_val.s = new string();
                            *attr_val.s = s; }
  AttrVal(const   void *p) { clear(); attr_type  = Attr_Ptr;     attr_val.p = p; }
  AttrVal(const AttrVal &a){ clear(); attr_type  = Attr_None;    copy(a); }

  ~AttrVal() { FreeData(); attr_type = Attr_None; }

  AttrVal &operator=(const AttrVal &a)
  {
    if(this == &a) return *this; 
    copy(a);
    return *this;
  }

  void copy(const AttrVal &a)
  {
    if(this == &a) return; 
    FreeData();
    attr_type = a.attr_type;
    if(attr_type == Attr_String)
      attr_val.s = new string(*a.attr_val.s);
    else
      attr_val  = a.attr_val;
  }

  void clear()
  {
    memset( &attr_val, 0, sizeof(attr_val_t) );
  }

  AttrVal &operator=(const int i)
      { FreeData(); attr_type = Attr_Integer; attr_val.i = i; return *this;}
  AttrVal &operator=(const double f)
      { FreeData(); attr_type = Attr_Real; attr_val.f = f; return *this;}
  AttrVal &operator=(const string &s)
      { FreeData(); attr_type = Attr_String;  attr_val.s = new string(s);
          return *this; }
  AttrVal &operator=(char s)
      { if(attr_type == Attr_String) *attr_val.s = s;
        else
        {
          FreeData();
          attr_type = Attr_String;
          attr_val.s = new string();
          *attr_val.s = s;
        }
        return *this;
      }
  AttrVal &operator=(const char *s)
      { if(attr_type == Attr_String) *attr_val.s = s;
        else
        { FreeData(); attr_type = Attr_String;  attr_val.s = new string(s); }
        return *this;
      }

  attr_t Type() const { return attr_type; }
  bool has_value() const { return attr_type != Attr_None; }
  bool operator==(const AttrVal &rhs) const;
  bool operator!=(const AttrVal &rhs) const { return !(*this == rhs); }
  bool operator< (const AttrVal &rhs) const;
  bool operator<=(const AttrVal &rhs) const;
  bool operator> (const AttrVal &rhs) const;
  bool operator>=(const AttrVal &rhs) const;
  AttrVal operator-() const;
  AttrVal operator+(const AttrVal &rhs) const;
  AttrVal operator-(const AttrVal &rhs) const;
  AttrVal operator*(const AttrVal &rhs) const;
  AttrVal operator/(const AttrVal &rhs) const;
  AttrVal operator%(const AttrVal &rhs) const;

  int         Int() const;
  double     Real() const;
  string   String(int width  = 0,
                  int pres   = -1,
                  long flags = 0,
                  char fch   = '\0') const;
  const void *Ptr()               const;

  void FreeData()
  {
    if(attr_type == Attr_String && attr_val.s)
      delete attr_val.s;
    attr_type = Attr_None;
    clear();  
  }
};

class AttrList
{
public:
  typedef pair<attr_id, AttrVal> Attr;
  typedef vector<Attr> ListType;

  typedef ListType::iterator               iterator;
  typedef ListType::const_iterator         const_iterator;
  typedef ListType::reverse_iterator       reverse_iterator;
  typedef ListType::const_reverse_iterator const_reverse_iterator;
  typedef ListType::size_type              size_type;
  typedef ListType::reference              reference;
  typedef ListType::const_reference        const_reference;

  AttrList& copy(const AttrList& rhs)
  { attrs = rhs.attrs; return *this; }
  AttrList& operator=(const AttrList& rhs)
  { attrs = rhs.attrs; return *this; }

  void                    clear()        { attrs.clear();           }
  iterator                begin()        { return attrs.begin();    }
  const_iterator          begin()  const { return attrs.begin();    }
  reverse_iterator       rbegin()        { return attrs.rbegin();   }
  const_reverse_iterator rbegin()  const { return attrs.rbegin();   }
  iterator                  end()        { return attrs.end();      }
  const_iterator            end()  const { return attrs.end();      }
  reverse_iterator         rend()        { return attrs.rend();     }
  const_reverse_iterator   rend()  const { return attrs.rend();     }
  bool                    empty()  const { return attrs.empty();    }
  size_type                size()  const { return attrs.size();     }
  reference               front()        { return attrs.front();    }
  const_reference         front()  const { return attrs.front();    }
  reference                back()        { return attrs.back();     }
  const_reference          back()  const { return attrs.back();     }

private:
  ListType attrs;

public:
//const ListType &attr() const   { return attrs; }

  static const string &Name( attr_id l) { return AttrNameMgr.query(l);    }
  static attr_id  AttrID( const string &s ) { return AttrNameMgr.query(s);}

  iterator attr(const string &s) { return attr( AttrNameMgr.query(s) ); }
  iterator attr(const attr_id l) { return(find(l)); }

  int IntAttr(const string &s)   { return IntAttr( AttrNameMgr.query(s) ); }
  int IntAttr(const attr_id l)
  {
    iterator i = find(l);
    return (i != attrs.end())? i->second.Int() : 0;
  }

  double FloatAttr(const string &s)
                              { return FloatAttr( AttrNameMgr.query(s) ); }
  double FloatAttr(const attr_id l)
  {
    iterator i = find(l);
    return (i != attrs.end())? i->second.Real() : std::numeric_limits<double>::quiet_NaN();
  }
  string StringAttr(const string &s) const
                              { return StringAttr( AttrNameMgr.query(s) ); }
  string StringAttr(const attr_id l) const
  {
    const_iterator i = find(l);
    return (i != attrs.end())? i->second.String() : string();
  }
  const void *PtrAttr(const string &s)
                              { return PtrAttr( AttrNameMgr.query(s) ); }
  const void *PtrAttr(const attr_id l)
  {
    iterator i = find(l);
    return (i != attrs.end())? i->second.Ptr() : NULL;
  }

  void set(const Attr &a)
  {
    iterator j = find(a.first);
    if(j != attrs.end())
      (*j).second = a.second;
    else
      attrs.push_back(a);
  }
  void set(const string  &s, const AttrVal &r )
                               { set( AttrNameMgr.add(s), r); }
  void set(const attr_id &l, const AttrVal &r )
                               { set(Attr(l,r)); }


  void set(const string  &id, const int i )
                              { set( AttrNameMgr.add(id), i); }
  void set(const attr_id l, const int i)
  {
    iterator j = find(l);
    if( j != attrs.end() )
      (*j).second  = i;
    else
      attrs.push_back( Attr(l, AttrVal(i)) );
  }

  void set(const string  &id, const double f )
                              { set( AttrNameMgr.add(id), f); }
  void set(const attr_id l, const double f)
  {
    iterator j = find(l);
    if( j != attrs.end() )
      (*j).second  = f;
    else
      attrs.push_back( Attr(l, AttrVal(f)) );
  }
  void set(const string  &id, const string &s )
                              { set( AttrNameMgr.add(id), s); }
  void set(const attr_id l, const string &s)
  {
    iterator j = find(l);
    if( j != attrs.end() )
      (*j).second  = s;
    else
      attrs.push_back( Attr(l, AttrVal(s)) );
  }
  void set(const string  &id, const char *s )
                              { set( AttrNameMgr.add(id), s); }
  void set(const attr_id l, const char *s)
  {
    iterator j = find(l);
    if( j != attrs.end() )
      (*j).second  = s;
    else
      attrs.push_back( Attr(l, AttrVal(s)) );
  }

  void set(const string  &id, const void *p )
                              { set( AttrNameMgr.add(id), p); }
  void set(const attr_id l, const void *p)
  {
    iterator j = find(l);
    if( j != attrs.end() )
      (*j).second  = p;
    else
      attrs.push_back( Attr(l, AttrVal(p)) );
  }

  iterator find(const string &s)
                              { return find( AttrNameMgr.add(s) ); }

  iterator find(const attr_id l)
  {
    for(iterator i=attrs.begin(); i != attrs.end(); i++)
      if(l == i->first )
        return i;
    return attrs.end();
  }

  const_iterator find(const string &s) const
                              { return find( AttrNameMgr.add(s) ); }

  const_iterator find(const attr_id l) const
  {
    for(const_iterator i=attrs.begin(); i != attrs.end(); i++)
      if(l == i->first )
        return i;
    return attrs.end();
  }

  bool has_attr(const string &s) const { return has_attr(AttrNameMgr.add(s)); }
  bool has_attr(const attr_id l) const { return (find(l) != attrs.end());     }
  
  void erase(iterator a) { attrs.erase(a); }
  void erase(iterator b,iterator e) { attrs.erase(b,e); }
  void erase(const string &id)  { erase( AttrNameMgr.add(id) ); }
  void erase(const attr_id l)
  {
    iterator i=find(l);
    if(i!=attrs.end())
      attrs.erase(i);
  }
};

#endif
