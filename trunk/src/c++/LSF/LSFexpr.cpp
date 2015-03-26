#include "LSF/LSFsymbol.h"
#include "LSF/LSFexpr.h"
#include "LSF/var_ops.h"
#include "LSF/SymbolInput.h"

#define DEBUG(x) 

LSFList* evaluate_lsfexpr(LSFBase* b, SymbolTable* sym_tab)
{
  // If b is NULL, can't do much
  if(!b) return NULL;

  DEBUG( cout << "Entered evaluate_lsfexpr(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  switch( b->lsf_type() ) 
  {
    case  LSF_GUARD: return evaluate_guard(b, sym_tab);
    case LSF_SWITCH: return evaluate_switch(b, sym_tab);
    case   LSF_EXPR: return evaluate_expr(b, sym_tab);
    default:         return b->List();
  }
}

bool is_lsfexpr(LSFBase* b)
{
  // If b is NULL, can't do much
  if(!b) return false;

//DEBUG( cout << "Entered is_lsfexpr(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  switch( b->lsf_type() ) 
  {
    case  LSF_GUARD: 
    case LSF_SWITCH: 
    case   LSF_EXPR: return true;
  }
  return false;
}

LSFList* evaluate_lsfcond(LSFBase* b, SymbolTable* sym_tab)
{
  // If b is NULL, can't do much
  if(!b) return NULL;

  DEBUG( cout << "Entered evaluate_lsfcond(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  switch( b->lsf_type() ) 
  {
    case LSF_GUARD:  return evaluate_guard(b, sym_tab);
    case LSF_SWITCH: return evaluate_switch(b, sym_tab);
    default:         return b->List();
  }
}

bool is_lsfcond(LSFBase* b)
{
  // If b is NULL, can't do much
  if(!b) return false;

//DEBUG( cout << "Entered is_lsfcond(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  switch( b->lsf_type() ) 
  {
    case LSF_GUARD: 
    case LSF_SWITCH: return true;
  }
  return false;
}

LSFList* evaluate_guard(LSFBase* b, SymbolTable* sym_tab)
{
  // If b is NULL, can't do much
  if(!b) return NULL;

  DEBUG( cout << "Entered evaluate_guard(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  if(b->lsf_type() != LSF_GUARD) return b->List();

  // if b is not a list or has no attributes don't bother
  if(!b->List() || !b->attrs()) return NULL;

  // First, find variables and make sure they're all around.
  AttrList::const_iterator i;

  // Get the test criteria and do basic testing on it.

  if((i = b->attrs()->find("TEST")) == b->attrs()->end()) return NULL;
  string s=resolveAttrVal((*i).second, sym_tab).String();

  // Get test variables.  If they're missing, the guard fails.
  AttrVal v1, v2;

  if((i = b->attrs()->find("VAR1")) == b->attrs()->end()) return NULL;
  v1 =resolveAttrVal((*i).second, sym_tab);

  if((i = b->attrs()->find("VAR2")) == b->attrs()->end()) return NULL;
  v2 =resolveAttrVal((*i).second, sym_tab);

  // Do Comparison

  bool t = false;

  if     (s == "==" || s == "EQ")  t = v1 == v2;
  else if(s == "!=" || s == "NE")  t = v1 != v2;
  else if(s == "<"  || s == "LT")  t = v1 <  v2;
  else if(s == ">"  || s == "GT")  t = v1 >  v2;
  else if(s == "<=" || s == "LE")  t = v1 <= v2;
  else if(s == ">=" || s == "GE")  t = v1 >= v2;

  DEBUG( cout << "LSFguard(" << b->name() << ", " << v1.String(); )
  DEBUG( cout << " " << s << " " << v2.String() << ") == "; )
  DEBUG( cout << (t? "true" : "false") << "<BR>" << endl; ) 
  if(t) return b->List();

  return NULL;
}

LSFList* evaluate_expr(LSFBase* b, SymbolTable* sym_tab)
{
  // If b is NULL, can't do much
  if(!b) return NULL;

  DEBUG( cout << "Entered evaluate_expr(" << b->name() << ")<BR>" << endl; )

  // If b isn't an expr, return b's List().
  if(b->lsf_type() != LSF_EXPR) return b->List();

  // if b is has no attributes don't bother
  if(!b->attrs()) return NULL;

  // First, find variables and make sure they're all around.
  AttrList::const_iterator i;

  // Get the test criteria and do basic testing on it.

  if((i = b->attrs()->find("OP")) == b->attrs()->end()) return NULL;
  string s=toUpper(resolveAttrVal((*i).second, sym_tab).String());

  // Get test variables.  If they're missing, the guard fails.
  AttrVal v1, v2, r;
  string result;

  if((i = b->attrs()->find("RESULT")) != b->attrs()->end())
  {
    result=resolveAttrVal((*i).second, sym_tab).String();
    if(!result.size()) return NULL;
  }

  if((i = b->attrs()->find("VAR1")) == b->attrs()->end()) return NULL;
  v1=resolveAttrVal((*i).second, sym_tab);

  DEBUG( if(v1.Type() == AttrVal::Attr_None) cout << "v1 missing! (" << (*i).second.String() << ")<BR>" << endl; )
  if(v1.Type() == AttrVal::Attr_None) return NULL;
  DEBUG( if(v1.Type() == AttrVal::Attr_Integer) cout << "v1 INT! (" << (*i).second.String() << ")<BR>" << endl; )
  DEBUG( if(v1.Type() == AttrVal::Attr_String) cout << "v1 STRING! (" << (*i).second.String() << ")<BR>" << endl; )

  if(!result.size() && (*i).second.Type() == AttrVal::Attr_String)
  {
    string n = (*i).second.String();
    if( n.size() && n[0] == '$')
      result = n.substr(1, n.size()-1);
    else
      return NULL;
  }

  if((i = b->attrs()->find("VAR2")) != b->attrs()->end())
    v2=resolveAttrVal((*i).second, sym_tab);

  DEBUG( if(v2.Type() == AttrVal::Attr_Integer) cout << "v2 INT!<BR>" << endl; )
  DEBUG( if(v2.Type() == AttrVal::Attr_String) cout << "v2 STRING!<BR>" << endl; )

  // Evaluate unary cases
  if( v2.Type() == AttrVal::Attr_None )
  {
    if     (s == "-" || s == "NEG") r = -v1;
    else if(s == "LEN")             r = (int)v1.String().size();
    else if(s == "UPCASE")          r = toUpper(v1.String());
    else if(s == "DOWNCASE")        r = toLower(v1.String());
    else if(s == "=" || s == "EQ" 
                     || s == "SET") r = v1;
    else if(s == "PRINT")           r = sym_tab->resolve_string(v1.String());
  }
  else
  {
    if     (s == "+" || s == "PLUS" ) r = v1 + v2;
    else if(s == "-" || s == "MINUS") r = v1 - v2;
    else if(s == "*" || s == "TIMES") r = v1 * v2;
    else if(s == "/" || s == "DIV"  ) r = v1 / v2;
    else if(s == "%" || s == "MOD"  ) r = v1 % v2;
    else if(s == "MAX")               r = max(v1,v2);
    else if(s == "MIN")               r = min(v1,v2);
  }

  DEBUG( if(r.Type() == AttrVal::Attr_Integer) cout << "r INT!<BR>" << endl; )
  DEBUG( if(r.Type() == AttrVal::Attr_String) cout << "r STRING!<BR>" << endl; )

  DEBUG(  cout << "LSFExpr: " << v1.String() << " " << s << " " 
               << v2.String() << " " << r.String() << " (" 
               << result << ")" << endl; )


  DEBUG( cout << "LSFexpr: setting " << result << " = " << r.String() << endl; )
#ifndef __GNUG__
  SymbolInput si(sym_tab);
  SymbolInput::result sr = si.set(result, r);
#else
  //assert(sr == SymbolInput::VALID);
  sym_tab->set(result,r);
#endif
  DEBUG( v1=sym_tab->query(result); )
  DEBUG( if(v1.Type() == AttrVal::Attr_Integer) cout << "r2 INT!<BR>" << endl; )
  DEBUG( if(v1.Type() == AttrVal::Attr_String) cout << "r2 STRING!<BR>" << endl; )

  return NULL;
}

LSFList* evaluate_switch(LSFBase* b, SymbolTable* sym_tab)
{
  // If b is NULL, can't do much
  if(!b) return NULL;

  DEBUG( cout << "Entered evaluate_switch(" << b->name() << ")<BR>" << endl; )

  // If b isn't a guard, return b's List().
  if(b->lsf_type() != LSF_SWITCH) return b->List();

  // if b is not a list or has no attributes don't bother
  if(!b->List() || !b->attrs()) return NULL;

  // First, find variables and make sure they're all around.
  AttrList::const_iterator i;

  // Get test variable.  If its missing, the switch fails.
  string v1, v2;

  if((i = b->attrs()->find("VAR1")) == b->attrs()->end()) return NULL;
  v1 = resolveAttrVal((*i).second, sym_tab).String();

  LSFList::iterator j;
  for( j = b->List()->begin(); j != b->List()->end(); ++j )
  {
    if( !*j ) continue;
    if( !(*j)->name().size() )
      return (*j)->List();
    v2 = sym_tab->resolve_string( (*j)->name() );
    DEBUG( cout << v1 << " == " << v2 << ' '; )
    if( !v2.size() || v1 == v2 ) 
    {
      DEBUG( cout << "passed<BR>" << endl; )
      return (*j)->List();
    }
    DEBUG( cout << "failed<BR>" << endl; )
  }

  return NULL;
}


#undef DEBUG
