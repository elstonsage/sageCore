#include <string>
#include "LSF/LSF.h"
#include "LSF/LSFmap.h"
#include "LSF/LSFtypes.h"
#include "LSF/LSFfactory.h"


LSFFactory* Factory;

LSFBase *LSFFactory::build_local( lsf_t t, const char* name,
                                  const char *, AList *l, LSFBase*)
{
  LSFBase *c = NULL;

  // Do storage type mapping
  if((int) t == -1 && l)
  {
    AList::iterator i = l->find(0);
    if( i != l->end() )
      switch( (*i).second.Type() )
      {
        case AttrVal::Attr_Integer: t = LSF_INT    << LSF_SHIFT; break;
        case AttrVal::Attr_Real:    t = LSF_REAL   << LSF_SHIFT; break;
        case AttrVal::Attr_String:  t = LSF_STRING << LSF_SHIFT; break;
        case AttrVal::Attr_Ptr:   
        case AttrVal::Attr_None:
        default:                  break;
      }
  }

  // Create the object
  switch (t>>LSF_SHIFT)
  {
    case LSF_MAP       :     c = new LSFmap ();       break;
    case LSF_REF       :   //c = new LSFref ();       break;
    case LSF_ITER      :
    case LSF_COMPOSITE :     
    case LSF_COMPONENT :     
    case LSF_INT       :
    case LSF_REAL      :
    case LSF_STRING    :
    case LSF_EXPR      :
    case LSF_SWITCH    :
    case LSF_GUARD     :     c = new LSFBase();
                             if(c) c->type(t);        break;
    case LSF_BASE      :     c = new LSFBase();       break;
    default            :     if(l || name) 
                               c = new LSFBase();     break;
  }

  if(!c) return NULL;

  if(name && *name)
    c->name( name );

  if(l)
    for(AList::iterator i = l->begin(); i != l->end(); i++)
      c->attrs(TRUE)->set( (*i) );

  return c;
}
