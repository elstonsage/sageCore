
#ifndef LSF_H
#define LSF_H

#include <string>
#include <list>
#include <iterator>
#include <iostream>
#include <map>
#include "globals/config.h"
#include "LSF/Attr.h"
#include "LSF/LSFtypes.h"

using namespace std;

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

class LSFBase;
class LSFList;
class LSFmap;
class LSFIterator;
class LSF_list_iterator;
class GenericFactory;

typedef unsigned long lsf_t;

template <class T = LSFBase>
class LSF_ptr
{
public:

  LSF_ptr(const LSF_ptr& p)  : my_ptr(p.my_ptr) { if(my_ptr) my_ptr->grab(); }
  LSF_ptr()                  : my_ptr(NULL)     { }
  LSF_ptr(T* p)              : my_ptr(p)        { if(my_ptr) my_ptr->grab(); }
  
  ~LSF_ptr() { if(my_ptr) my_ptr->release(); }

  LSF_ptr& operator=(const LSF_ptr& p)
  {
    if(my_ptr == p.my_ptr) return *this;

    T* t = const_cast<T*>(p.my_ptr);
    
    if(t)        t->grab();
    if(my_ptr)   my_ptr->release();
    my_ptr = t;
    
    return *this;
  }
  
  LSF_ptr& operator=(T* p)
  {
    if(my_ptr == p) return *this;
    
    if(p)      p->grab();
    if(my_ptr) my_ptr->release();
    my_ptr = p;
    
    return *this;
  }
  
  T& operator*() { return *my_ptr; }
  const T& operator*() const { return *my_ptr; }

  T* operator->() { return my_ptr; }
  const T* operator->() const { return my_ptr; }

  operator T*() const { return my_ptr; }

  bool operator<(const LSF_ptr<T>& p) const
  { return my_ptr < p.my_ptr; }
  
private:

  T* my_ptr;
};

class LSFDef
{
public:
  typedef AttrList                AList;  // List of attributes
  typedef AList::Attr             Attr;   // A single attribute
  typedef std::map<string, LSF_ptr<LSFBase>,
          less<string> >      pmap;   // String, LSF pointer map;
  typedef list< LSF_ptr<LSFBase> >       CList;  // List of LSF pointers

  virtual ~LSFDef() { }
};

class LSFVisitor
{
public:
  virtual ~LSFVisitor() { }
  virtual void Process(LSFBase      *p);
protected:
  virtual void def(LSFBase *) { }
  LSFVisitor() {}
};

class LSF_list_iterator : public LSFDef
{
public:
  typedef list<LSF_ptr<LSFBase> > CList;
  typedef CList::iterator                              iterator;
  typedef CList::const_iterator                  const_iterator;
  typedef CList::reverse_iterator              reverse_iterator;
  typedef CList::const_reverse_iterator  const_reverse_iterator;

  LSF_list_iterator(LSFBase *c);
  LSF_list_iterator(LSFBase *c, iterator &x); 
  virtual bool eq(const LSF_list_iterator *x) const 
                                     { return (i == x->i); }
  virtual bool eq(const iterator &x) const { return (x==i); };
  virtual LSFBase *deref()           { return *i;     }
  virtual LSFBase *deref() const     { return *i;     }
  virtual LSFBase *operator->()      { return *i;     }
  virtual LSFBase *operator->() const{ return *i;     }
  virtual iterator current()         { return  i;     }
  virtual iterator next()            { ++i; return i; } 
  virtual iterator prev()            { --i; return i; }
  virtual void set(const iterator &x) {i=x;}
  virtual void set(iterator &x) {i=x;}

protected:
  iterator i;
  LSFList *l;
};

class LSFIterator //: public bidirectional_iterator<LSFBase *, ptrdiff_t>
{
private:
  LSF_list_iterator *it;
public:
  typedef LSFBase *& reference;
  typedef list<LSF_ptr<LSFBase> >::iterator iterator;
  LSFIterator();
  LSFIterator(LSFBase *x);
  LSFIterator(const LSFIterator &x);
  LSFIterator(LSF_list_iterator *x);
  LSFIterator &operator= (const iterator &x)  { it->set(x); return *this; }
  bool operator==(const LSFIterator &x) const { return (it->eq(x.it));    }
  bool operator==(const LSFDef::CList::iterator &x) const
                                              { return (it->eq(x));       }
  LSFBase     *operator*()                    { return it->deref();       }
  LSFBase     *operator*() const              { return it->deref();       }
  LSFBase     *operator->()                   { return it->deref();       }
  LSFBase     *operator->() const             { return it->deref();       }
  LSFIterator &operator++()                   { it->next(); return *this; }
  LSFIterator  operator++(int)                { LSFIterator tmp = *this;
                                                ++*this;
                                                return tmp; }
  LSFIterator& operator--()                   { it->prev(); return *this; }
  LSFIterator  operator--(int)                { LSFIterator tmp = *this;
                                                --*this;
                                                return tmp; }
};

// Virtual base for composite objects w/ memory management and simple
// error flags
class LSFBase
{
public:
  typedef LSFDef::AList AList;
  typedef LSFDef::Attr  Attr;

  // Constructor/Destructor
  LSFBase(const char *n = "", lsf_t t = 0) { init(n, t); }
  LSFBase(const string &n, lsf_t t = 0)    { init(n, t); }
  virtual ~LSFBase();

  // Memory Management
  void    grab() { ++_ref;  }
  void release() { if(--_ref <= 1 && _ref != -99 ) 
                   { _ref  = -99; delete this; } }
  int      ref() { return _ref-1; }

  bool selfish()       { return _selfish; }
  bool selfish(bool s) { return (_selfish = s); }
  static long int existing()  { return existing_components; }
  static LSFBase *getLSFObjs();
  static bool isValidLSF(LSFBase *);
      
  // Name & Type handlers
  const string &name(const char *n) 
                                  { if(n) _name = n;
                                    return _name; }
  const string &name(const string &n) 
                                  { _name = n;
                                    return _name; }
  const string &name()      const { return _name; }
  lsf_t type()              const { return (_type << LSF_SHIFT) |
                                           (_local_type & LOCAL_MASK); }
  lsf_t type(const lsf_t t)       { _type = (t>>LSF_SHIFT); 
                                    _local_type = t & LOCAL_MASK; 
                                    return(t); }
  lsf_t lsf_type()          const { return _type; }
  lsf_t lsf_type(const lsf_t t)   { return (_type = t); }
  lsf_t local_type()        const { return _local_type; }
  lsf_t local_type(const lsf_t t) { return (_local_type = t); }

  virtual AList           *attrs(bool create = FALSE);
  virtual AList           *attrs(bool = FALSE) const { return _attrs; }
  virtual LSFList          *List(bool create = FALSE);
  virtual const LSFList    *List(bool = FALSE) const { return _list;  }
  virtual LSF_list_iterator *local_iterator()  { return new LSF_list_iterator(this); }

  virtual void free_attrs();
  virtual void free_list();
  
  // Visitor entry points
  virtual void Accept(LSFVisitor *v)   { v->Process(this); }
  virtual void AcceptAll(LSFVisitor *);

  enum m_state { goodbit, failbit, badbit };

  void     clear()                { _state  = 0;      }
  void     clear(const int state) { _state &= ~state; }
  void  setstate(const int flag)  { _state |= flag;   }

  int       good() const { return _state == 0; }
  int       fail() const { return _state & (badbit|failbit); }
  int        bad() const { return _state & badbit; }
  int    rdstate() const { return _state; }
  operator void*() const { return fail() ? (void*)0 : (void*)(-1); }
  int  operator!() const { return fail(); }

protected:
  unsigned int    _type;               // Object type
  unsigned int    _local_type;         // local type tag
  AList          *_attrs;
  LSFList        *_list;
private:
  void init(const string &s, lsf_t t);
  static long int existing_components;

  int    _ref;                // Number of times ref'd
  bool   _selfish;            // Can we be shared.
  string _name;
  int    _state;              // Flags
};

class LSFList
{
public:
  typedef list<LSF_ptr<LSFBase> >        CList;

  typedef CList::iterator                              iterator;
  typedef CList::const_iterator                  const_iterator;
  typedef CList::reverse_iterator              reverse_iterator;
  typedef CList::const_reverse_iterator  const_reverse_iterator;

  typedef CList::reference       reference;
  typedef CList::const_reference const_reference;
  typedef CList::difference_type difference_type;
  typedef CList::value_type      value_type;
  typedef CList::size_type       size_type;

  LSFList() {}
  LSFList(const LSFList& lst) { copy(lst); }

  virtual ~LSFList() { }

  LSFList &copy(const LSFList& rhs)
  { clist = rhs.clist; return *this; }
  LSFList& operator=(const LSFList& rhs)
  { clist = rhs.clist; return *this; }

  iterator                begin()       { return clist.begin(); }
  const_iterator          begin() const { return clist.begin(); }
  reverse_iterator       rbegin()       { return clist.rbegin();}
  const_reverse_iterator rbegin() const { return clist.rbegin();}
  iterator                  end()       { return clist.end();   }
  const_iterator            end() const { return clist.end();   }
  reverse_iterator         rend()       { return clist.rend();  }
  const_reverse_iterator   rend() const { return clist.rend();  }

  bool                    empty() const { return clist.empty();  }
  size_type                size() const { return clist.size();   }
//size_type            max_size() const { return clist.max_size();}
  reference               front()       { return clist.front();  }
  const_reference         front() const { return clist.front();  }
  reference                back()       { return clist.back();   }
  const_reference          back() const { return clist.back();   }

  void erase(iterator position)
  {
    clist.erase(position);
  }
  
  void erase(iterator first, iterator last)
  {
    for(; first != last; ++first)
      erase(first);
  }
      
  void add(LSFBase *item)
  { 
    if(!item) return;
    clist.push_back(item);
  }

  void erase(LSFBase *item)
  {
    if(!item) return;
    iterator i = std::find (clist.begin(), clist.end(), item);
    if ( i != clist.end() )
      clist.erase( i );
  }

  void remove(LSFBase *item) { erase(item); }

  LSFBase *find(LSFBase *item)
  {
    if(!item) return NULL;
    iterator i = std::find (clist.begin(), clist.end(), item);
    if ( i != clist.end() )
      return *i;
    return NULL;
  }

  iterator insert(iterator position, LSFBase *item) 
  { if(!item) return end(); return clist.insert(position, item); }

  iterator insert(iterator position) 
  { LSFBase *item = new LSFBase(); return clist.insert(position, item);  }

  void push_front(LSFBase *item)
  { 
    if(!item) return;
    clist.push_front(item);
  }
  
  void push_back(LSFBase *item)
  { 
    if(!item) return;
    clist.push_back(item);
  }
  
  void pop_front() { erase(begin()); }
  void pop_back()  { iterator tmp = end(); erase(--tmp); }
  void clear() { erase( begin(), end() ); }

protected:
  CList clist;
};

class GenericFactory : public LSFBase
{
public:
  typedef LSFDef::AList   AList;
 
  GenericFactory(const char *n = "") : LSFBase(n) { _type = LSF_FACTORY; }

  virtual ~GenericFactory();
  virtual const string &type_map( lsf_t );
  virtual void add_factory( GenericFactory *f, Type_mapping *m );
  virtual void add_factory( GenericFactory *f, LSF_mapping *m );
  virtual void add_factory( GenericFactory *f, lsf_t t, Local_mapping *m );
  virtual void dump_map(ostream &);

  LSFBase *build( lsf_t t, const char *name = NULL, const char *type = NULL,
                   AList *l = NULL, LSFBase *p = NULL); 
  LSFBase *build( lsf_t lsf, lsf_t local, const char *name = NULL,
                  const char *type = NULL, AList *l = NULL,
                  LSFBase *p = NULL)
               { return build (lsf << LSF_SHIFT | local, name, type, l, p); } 
  LSFBase *build(const char* name, const char *type = NULL, AList *l = NULL,
                 LSFBase* p = NULL)
               { return build( (lsf_t) -1, name, type, l, p); }
        
protected:
  typedef AList::iterator                            a_iter;
  typedef std::map<lsf_t, GenericFactory *, less<lsf_t> > factorymap;
  typedef std::map<string, lsf_t, less<string> >          typemap;

  virtual LSFBase *build_local( lsf_t t, const char* name = NULL,
                                const char *type = NULL, AList *s = NULL, 
                                LSFBase* p = NULL );

  factorymap _factory_map;
  typemap _type_map;
};

inline AttrVal attr_value(const LSFBase *b, const string &name, attr_id attr)
{
  if(!b) return AttrVal();

  if(name.size() && name != toUpper(b->name()))
    return AttrVal();

  AttrList::const_iterator a;
  if( b->attrs() && (a=b->attrs()->find(attr)) != b->attrs()->end())
    return a->second;
  return AttrVal();
}

inline AttrVal attr_value(const LSFBase *b, const string &name, const string &attr)
{
  return attr_value(b,name, AttrNameMgr.query(attr));
}

inline AttrVal attr_value(const LSFBase *b, const string &attr)
{
  return attr_value(b, "", AttrNameMgr.query(attr));
}

inline AttrVal attr_value(const LSFBase *b, attr_id attr)
{
  return attr_value(b, "", attr);
}

inline bool has_attr(const LSFBase *b, const string &name, attr_id attr)
{
  if(!b) return false;

  if(name.size() && name != toUpper(b->name()))
    return false;

  if( b->attrs() && b->attrs()->find(attr) != b->attrs()->end())
    return true;
  return false;
}

inline bool has_attr(const LSFBase *b, const string &name, const string &attr)
{
  return has_attr(b,name, AttrNameMgr.query(attr));
}

inline bool has_attr(const LSFBase *b, const string &attr)
{
  return has_attr(b, "", AttrNameMgr.query(attr));
}

#endif
