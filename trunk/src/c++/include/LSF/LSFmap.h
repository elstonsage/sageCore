#ifndef __LSF_MAP_H_
#define __LSF_MAP_H_

#include <map>
#include "LSF/LSF.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

class LSFmap : public LSFBase
{
public:
  typedef LSFDef::pmap                 pmap;

  typedef pmap::key_type               key_type;
  typedef pmap::key_compare            key_compare;
  typedef pmap::reference              reference;
  typedef pmap::const_reference        const_reference;
  typedef pmap::iterator               iterator;
  typedef pmap::const_iterator         const_iterator;
  typedef pmap::reverse_iterator       reverse_iterator;
  typedef pmap::const_reverse_iterator const_reverse_iterator;
  typedef pmap::size_type              size_type;
  typedef pmap::difference_type        difference_type;

  LSFmap(const char *n = "") : LSFBase(n) { _type = LSF_MAP; }

  LSFmap(const LSFmap &m)
  {
    name( m.name() );
    _type = LSF_MAP;
     for(const_iterator i=m.begin(); i != m.end(); ++i)
       add( i->first, i->second );
  }
  
  virtual ~LSFmap()  { }

  virtual void Accept(LSFVisitor *v) { v->Process(this); }
  virtual void AcceptAll(LSFVisitor *v) 
  { for(iterator i=begin(); i!=end(); i++) (*i).second->Accept(v); }
  
// accessors:

  key_compare          key_comp() const { return p_map.key_comp(); }
  iterator                begin()       { return p_map.begin();    }
  const_iterator          begin() const { return p_map.begin();    }
  iterator                  end()       { return p_map.end();      }
  const_iterator            end() const { return p_map.end();      }
  reverse_iterator       rbegin()       { return p_map.rbegin();   }
  const_reverse_iterator rbegin() const { return p_map.rbegin();   }
  reverse_iterator         rend()       { return p_map.rend();     }
  const_reverse_iterator   rend() const { return p_map.rend();     }
  bool                    empty() const { return p_map.empty();    }
  size_type                size() const { return p_map.size();     }

  LSFBase *find(const string &s) const
  {
    if(!s.length()) return NULL;
    const_iterator i=p_map.find(s);
    if(i!=p_map.end()) return (*i).second;
    return NULL;
  }

  void   add(const string &s, LSFBase *v) { if (v) set(s, v, true); }
  bool   set(const string &s, LSFBase *v, bool b = false);
  void remove(const string &s);

  LSFBase *operator[](const string &s)    { return find(s); }
  size_type count(const string& x) const { return p_map.count(x); }
  
private:
  pmap p_map;
};

#endif
