#include <cassert>
#include "LSF/LSFiter.h"
#include "LSF/LSFexpr.h"

#define DEBUG(x)

LSFList::iterator ELSF_list_iterator::current()
{
  if(!syms || !l) return i;
  if(!sstack.size() && i != l->end() && (*i) 
     && ( (*i)->lsf_type() == LSF_ITER || is_lsfcond(*i)))
  {
    saveState();
    return initIterator( *i );
  }
  else if(sstack.size() && i != sstack.top().l->end() && *i 
     && ( (*i)->lsf_type() == LSF_ITER || is_lsfcond(*i)))
    return initIterator( *i );
  
  return i;
}

LSFList::iterator ELSF_list_iterator::next()           
{ 
  ++i;
  if(!syms) return i;

  if(!sstack.size()) 
    return current();
  iter_state &cstate = sstack.top();
  
  i=++cstate.i; 

  if ( i != cstate.l->end() && *i 
     && ( (*i)->lsf_type() == LSF_ITER || is_lsfcond( *i )))
    return current();

  if(cstate.i != cstate.l->end()) return i; // not the end,

  if(cstate.i == cstate.l->end() && sstack.size() <= 1 &&
    cstate.current >= cstate.stop) 
    return (i=l->end());

  if(cstate.current < cstate.stop)     // Another turn around the track
  {
    cstate.current++;
    DEBUG( cout << "Iteration #" << cstate.current << "<BR>" << endl; )
    if(cstate.looper) 
      cstate.looper->attrs(TRUE)->set(0, cstate.current);
    i=cstate.i=cstate.l->begin();
    return current();
  }

  // Or else it is time to return to the previous iteration
  sstack.pop();
  i = sstack.top().i;
  next();
  return current();
}   

LSFList::iterator ELSF_list_iterator::initIterator( LSFBase *n )
{
  if (!n || !n->attrs() ) return( next() );

  iter_state nstate;

  if(n->lsf_type() == LSF_ITER)
  {
    if (!n->List() || !n->List()->size() || !n->attrs() || !n->name().size() ) 
      return( next() );

    nstate.start =
      atoi(syms->resolve_string(n->attrs()->StringAttr("START")).c_str());
    nstate.stop  = 
      atoi(syms->resolve_string(n->attrs()->StringAttr("STOP" )).c_str());
    if(nstate.start > nstate.stop) return( next() );
      nstate.current = nstate.start;
    nstate.looper  = syms->find( n->name() );
    if(!nstate.looper) 
      nstate.looper  = syms->add(n->name() );
    nstate.looper->attrs(TRUE)->set(0, nstate.start);
    nstate.i = n->List()->begin();
    nstate.l = n->List();
    DEBUG( cout << "Starting iteration for " << n->name() << "  ("; )
    DEBUG( cout << nstate.start << "," << nstate.stop << ")<BR>" << endl; )
  }
  else if( is_lsfcond(n) )
  {
    LSFList *lst = evaluate_lsfcond(n, syms);
    if (!lst || !lst->size() ) return( next() );
    DEBUG( cout << "Adding guarded expr " << n->name() << "<BR>"; )
    nstate.start = 0;
    nstate.stop = 0;
    nstate.current = 0;
    nstate.looper = NULL;
    nstate.i = lst->begin();
    nstate.l = lst;
  }
  else
    return next();

  sstack.push( nstate );
  i=nstate.i;
  return current();
}

void ELSF_list_iterator::saveState()
{
  iter_state nstate;
  nstate.start = 0;
  nstate.stop  = 0;
  nstate.current = 1;
  nstate.looper  = NULL;
  nstate.l       = l;
  nstate.i       = i;
  sstack.push( nstate );
}

ELSFIterator::ELSFIterator(LSFBase *x, SymbolTable *s)
{ 
  it = new ELSF_list_iterator(x,s); 
}

ELSFIterator::ELSFIterator(ELSF_list_iterator *x)  
{ 
  it = x;
}
	
#undef DEBUG
