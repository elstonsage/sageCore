#ifndef _DISPLAY_ITER_H
#define _DISPLAY_ITER_H

#include <stack>
#include "LSF/LSF.h"
#include "LSF/LSFsymbol.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

  
class ELSF_list_iterator : public LSF_list_iterator
{
public:
  ELSF_list_iterator(LSFBase *c, SymbolTable *s = NULL) 
                        : LSF_list_iterator(c), syms(s) { current(); }
  virtual ~ELSF_list_iterator() {}
  virtual iterator next();
protected:
  struct iter_state 
  { 
  	int start, stop, current; 
  	LSFBase *looper;
  	LSFList *l; 
  	LSF_list_iterator::iterator i; 
	// Needed since VC++ is rather dumb
  	bool operator==(const iter_state &) const { return false; }
  	bool operator!=(const iter_state &) const { return false; }
  	bool operator< (const iter_state &) const { return false; }
  	bool operator> (const iter_state &) const { return false; }
  };
  virtual iterator current();
  virtual iterator initIterator( LSFBase *n );
  void saveState();

  SymbolTable *syms;

#ifdef OLD_STL
  typedef stack< list< iter_state > > StateStack;
#else
  typedef stack< iter_state, list<iter_state> > StateStack;
#endif
  StateStack sstack;
};

class ELSFIterator //: public bidirectional_iterator<LSFBase *, ptrdiff_t>
{
private:
  ELSF_list_iterator *it;
public:
  typedef LSFBase *& reference;
  typedef list<LSFBase *>::iterator iterator;
  ELSFIterator();
  ELSFIterator(LSFBase *x, SymbolTable *s);
  ELSFIterator(ELSF_list_iterator *x);
  bool operator==(const ELSFIterator &x) const   { return  (it->eq(x.it)); }
  bool operator==(const LSFList::iterator &x) const { return  (it->eq(x));    }
  bool operator!=(const ELSFIterator &x) const   { return !(it->eq(x.it)); }
  bool operator!=(const LSFList::iterator &x) const { return !(it->eq(x));    }

  LSFBase     *operator->() const    { return it->deref();         }
  LSFBase     *operator*() const     { return it->deref();         }
  ELSFIterator& operator++()         { it->next(); return *this;   }
  ELSFIterator  operator++(int)      { ELSFIterator tmp = *this;
                                          ++*this;
                                          return tmp; }
  ELSFIterator& operator--()         { it->prev(); return *this;   }
  ELSFIterator  operator--(int)      { ELSFIterator tmp = *this;
                                          --*this;
                                          return tmp; }
};

#endif
