#include <string>
#include <stdlib.h>
#include <cassert>
#include "LSF/LSF.h"
#include "LSF/LSFsymbol.h"
#include "LSF/SymbolInput.h"
#include "LSF/var_ops.h"
#include "LSF/LSFfactory.h"

#define DEBUG(x)

using namespace std;

void expandArray(LSFBase *root, LSFBase *tmpl, int size)
{
  if( !root || !tmpl || size < 1) return;

  int start = root->List(TRUE)->size();

  LSFBase *t;
  int j;

  if ( start > size ) 
  {
    // We don't touch elements that we don't use
#if 0 
    for(j=1; j<size;  j++,i++); // Seek end of usable array
    for( ;   j<start; j++,i++)  // Set the correct default values
      setDef( (*i), true );     // Array is already sized
#endif
    return;
  }
  
  for(j = start; j < size; j++)
  {
    t = Factory->build( tmpl->type() );
    root->List()->add(t);
  }
}

void addIndexSymbols( const LSFBase *root, SymbolTable *syms)
{
  // If NULLs or no list elements, don't do anything.
  if( !root || !syms || !root->List() || !root->List()->size()  ) return;

  LSFList::const_iterator i = root->List()->begin();

  for( ; i != root->List()->end() ; i++ )  // For each element in the List
  {
    if( !(*i)->name().size() ) continue;       // Unnamed object

    syms->add( (*i)->name(), *i );
  }
}

LSFBase *extractHead(LSFBase *b)
{
  LSFList::iterator i;
  
  if(!b || !b->List() || !b->List()->size() ) return b;
  i = b->List()->begin();
  (*i)->grab();
  assert( b->List()->find( (*i) ) != NULL );
  LSFBase *head = *i;
  b->List()->remove( head );
  assert( b->List()->find( head ) == NULL );
  return head;
}

void addSymbol(LSFBase *v, SymbolTable *vars)
{
  if(!v || !v->name().size()) return;     // Can't do anything with this
  
  string n = v->name();

  DEBUG( cout << "addSymbol(" << n << ")<BR>" << endl; )
  SymbolTable::AttrPair parts = vars->resolveID( n );

  LSFBase *oldvar = parts.first;

  // If there isn't a variable create one
  if(!oldvar)
  {
    DEBUG( cout << "1<BR>" << endl; )
    oldvar = Factory->build( v->type() );

    // Last resort
    if(!oldvar) 
      oldvar = new LSFBase();
    oldvar = vars->add( n, oldvar );
    DEBUG( cout << "2<BR>" << endl; )
  }

  // This means that something is horribly wrong with the variable
  // and it cannot be created
  if(!oldvar) return;

  // Copy attributes to new variable if needed
  if(v->attrs())
  {
    DEBUG( cout << "3<BR>" << endl; )
    bool val_saved = false;              // Do we save the value
    AttrVal oldval;                      // Storage for saved value
    DEBUG( cout << "3.1<BR>" << endl; )

    // Preserve attribute VALUE if it exists
    if( oldvar->attrs() ) 
    {
      DEBUG( cout << "4<BR>" << endl; )
      AttrList::iterator o;
      o = oldvar->attrs()->find(0); // Find the old value
      if( o != oldvar->attrs()->end() )
      {
        oldval = (*o).second;            // Save and mark
        val_saved = true;
      }
      DEBUG( cout << "5<BR>" << endl; )
    }
 
    DEBUG( cout << "6<BR>" << endl; )
    // Remove old attributes
    oldvar->attrs(true)->clear();  
    DEBUG( cout << "7<BR>" << endl; )

    AttrList::iterator i;
    for( i = v->attrs()->begin(); i != v->attrs()->end(); i++)
    {
      if((*i).second.Type() == AttrVal::Attr_String)
        oldvar->attrs()->set( (*i).first, 
                              vars->resolve_string((*i).second.String()) );
      else
        oldvar->attrs()->set(*i);
      DEBUG( cout << "8<BR>" << endl; )
    }

    if(val_saved)
      SymbolInput::set( oldvar, 0, oldval );  // Restore value
    else
      setDef( oldvar );                       // else throw it to setDef
      DEBUG( cout << "9<BR>" << endl; )
  }
}

AttrVal resolveAttrVal(const AttrVal &v, const SymbolTable* sym_tab)
{
  if(!sym_tab) return v;

  if(v.Type()==AttrVal::Attr_String && v.String().size() && v.String()[0]=='$')
    return sym_tab->query(v.String().substr(1, v.String().length()-1));

  return v;
}

AttrVal find_and_resolveAttrVal(const string &name, const SymbolTable *sym_tab)
{
  if(!sym_tab) return AttrVal();

  SymbolTable::AttrPair r = sym_tab->resolveID(name);

  if(!r.first || !r.first->attrs() || (signed)r.second < 0)
    return AttrVal();

  AttrList::iterator a = r.first->attrs()->find(r.second);

  if( a == r.first->attrs()->end() )
    return AttrVal();

  return resolveAttrVal( (*a).second, sym_tab );
}

bool setDef(LSFBase *b, bool clear )
{
  if(!b || !b->attrs()) return false;
  AttrList::iterator a = b->attrs()->find("DEF");
  if(a != b->attrs()->end() )
  {
    b->attrs()->set( 0, (*a).second );
    return true;
  }
  else if(clear)
    b->attrs()->erase( (unsigned)0 );

  return false;
}

#undef DEBUG
