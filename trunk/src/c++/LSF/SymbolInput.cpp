#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "LSF/LSFsymbol.h"
#include "LSF/SymbolInput.h"
#include "LSF/LSFiter.h"
#include "LSF/var_ops.h"

SymbolInput::SymbolInput(SymbolTable* s) : sym(s)
{
  if(sym) sym->grab();
}

SymbolInput::~SymbolInput()
{
  if(sym) sym->release();
}

SymbolInput::result SymbolInput::set(const string& symbol, const string& value)
{
  if(!sym) return NO_SYM_TAB;
  SymbolTable::AttrPair a = sym->resolveID(symbol);

  return set(a.first, a.second, value);
}

SymbolInput::result SymbolInput::set(const string& symbol, const AttrVal& attr)
{
  if(!sym) return NO_SYM_TAB;
  SymbolTable::AttrPair a = sym->resolveID(symbol);

  return set(a.first, a.second, attr);
}

SymbolInput::result SymbolInput::set(LSFBase* b, const attr_id& aid,
                      const string& value)
{
  if(!b || (signed) aid == -1) return NOT_IN_ST;

  string val = value;

  if(aid != 0)
  {
    (*(b->attrs(true)->attr(aid))).second = val;
    return VALID;
  }

  if( !val.size() && b->attrs() )
  {
    AttrList::iterator i = b->attrs()->attr("FALSE");
    if( i  != b->attrs()->end())
    {
      val = (*i).second.String();
    }
  }

  if (MapString(b, aid, val)) return VALID;

  AttrVal v;

  switch(b->lsf_type())
  {
    case LSF_INT:
      if(!val.size() ) 
        if(setDef(b, true)) return VALID; else return test_required(b);

      if((v = stringtol(val)).Type() == AttrVal::Attr_None)
        return CONV_ERROR;
      break;

    case LSF_REAL :
      if(!val.size() ) 
        if(setDef(b, true)) return VALID; else return test_required(b);

      if((v = stringtof(val)).Type() == AttrVal::Attr_None)
        return CONV_ERROR;
      break;

    default       :
      v = val; 
      break;
  }
  return set(b, aid, v);
}

SymbolInput::result SymbolInput::set(LSFBase* b, const attr_id& aid,
                      const AttrVal& value)
{
  if(!b || (signed) aid == -1) return NOT_IN_ST;

  AttrVal v = value;

  if(aid != 0)
  {
    (*(b->attrs(true)->attr(aid))).second = v;
    return VALID;
  }

  // If no attrs
  if(!b->attrs())
  {
    if(v.String().size())
      { b->attrs(true)->set(0,v); return VALID; }
  }    
  else
    if(v.Type() == AttrVal::Attr_String && !v.String().size() )
    {
      AttrList::iterator i = b->attrs()->attr("FALSE");
      if ( i != b->attrs()->end() )
        v = (*i).second;
    }

  AttrList::iterator i = b->attrs(TRUE)->end();
  AttrList::iterator end_iter = i;

  switch(b->lsf_type())
  {
    case LSF_INT  :
      if(v.Type() != AttrVal::Attr_Integer)
        return set(b, aid, v.String()); 
      break;
    case LSF_REAL :
      if(v.Type() != AttrVal::Attr_Real)
        return set(b, aid, v.String());
      break;
    default       :
      if(v.String().size())                         // If there's a string.
      {
        i = b->attrs()->attr("LEN");
        if(i != end_iter && (unsigned) (*i).second.Int() < v.String().size()) 
          b->attrs()->set(0, v.String().substr(0,(*i).second.Int()));
        else
        {
          b->attrs()->set(0, v);
        }
        return VALID;
      }
      
      if(setDef(b, false))
        return VALID;
      else
      {
        b->attrs()->set(0, v);     // If nothing there, blank it.
        return test_required(b);
      }
  }

  bool valid = true;

  i = b->attrs()->attr("MIN");
  if(i != end_iter && (*i).second.Real() >  v.Real()) valid = false;

  i = b->attrs()->attr("EMIN");
  if(i != end_iter && (*i).second.Real() >= v.Real()) valid = false;

  i = b->attrs()->attr("MAX");
  if(i != end_iter && (*i).second.Real() <  v.Real()) valid = false;

  i = b->attrs()->attr("EMAX");
  if(i != end_iter && (*i).second.Real() <= v.Real()) valid = false;

  if(valid)
  {
    b->attrs()->set(0,v);
    return VALID;
  }
  else
  {
    if(b->attrs()->attr(0) == end_iter)
      setDef(b, false);
    return OUT_OF_BOUND;
  }
}

bool SymbolInput::MapString(LSFBase* b, const attr_id& aid, string& val)
{
  if(!b->attrs() || !b->List() || !b->List()->size() ||
      b->attrs()->attr("MAP") == b->attrs()->end())
    return false;

  LSFList::iterator item = b->List()->begin();
  if((*item)->name().size())                   // If string mapped
  {
    for( ; item != b->List()->end() && !((*item)->name() == val); item++);

    if(item != b->List()->end() && (*item)->attrs() &&  
       (*item)->attrs()->attr(0) != (*item)->attrs()->end())
    {
      set(b, aid, (*((*item)->attrs()->attr(0))).second);
      return true;
    }
    setDef(b, false);
    return true;
  }

  // Value mapped
  AttrVal v = stringtol(val);
  if(v.Type() == AttrVal::Attr_Integer &&
     (unsigned) v.Int() <= b->List()->size())
  {
    for(int counter = 1;
        counter < v.Int() && item != b->List()->end();
        counter++, item++);

    if(item != b->List()->end() && (*item) && (*item)->attrs() &&  
       (*item)->attrs()->attr(0) != (*item)->attrs()->end())
    {
      set(b, aid, (*((*item)->attrs()->attr(0))).second);
      return true;
    }
  }
  setDef(b, false);
  return true;
}

SymbolInput::result SymbolInput::test_required(LSFBase* b)
{
  if(!b->attrs() || b->attrs()->find("REQUIRED") == b->attrs()->end())
    return NO_DEFAULT;
  else
    return REQUIRED;
}
