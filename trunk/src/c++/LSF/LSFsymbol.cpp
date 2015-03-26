#include <utility>
#include "LSF/Attr.h"
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "LSF/LSFsymbol.h"

LSFBase* SymbolTable::replace(const string& s, LSFBase* b, bool create)
{
  unsigned int i = 0;
  return evalString(s, i, create, b).first;
}

void SymbolTable::remove(const string& s)
{
  unsigned int i = 0;

  string v = getIdentifier(s,i);

  if(v.size()) 
    st.remove(v);
}

LSFBase* SymbolTable::find(const string& s) const
{
  unsigned int i = 0;

  return getPrimary(s, i);
}

LSFBase* SymbolTable::set(const string& s, const AttrVal& v, bool create)
{
  unsigned int i = 0;

  AttrPair p;

  top_level = false;                    // Don't worry about garbage at end

  p.first = getPrimary(s, i);           // Find end of identifier string.

  if(s.size() <= i || s[i] != '.')      // No attribute specified
    p.second = 0;
  else
    p.second = getAttrID(s, ++i, create);

  top_level = true;

  if(!robust && i < s.size())            // Stuff at end problem again.
    return NULL;
  return setPair(s, p, v, create);
}

LSFBase* SymbolTable::set(const string& s, const string& attr,
                          const AttrVal& v, bool create)
{
  unsigned int i;

  AttrPair p;

  top_level = false;                    // Don't worry about garbage at end

  p.first = getPrimary(s, i = 0);

  p.second = getAttrID(attr, i = 0, create);

  top_level = true;

  if(!robust && i < s.size())            // Stuff at end problem again.
    return NULL;

  return setPair(s, p, v, create);
}

AttrVal SymbolTable::query(const string& s) const
{
  return evalPair(resolveID(s));
}

AttrVal SymbolTable::query(const string& var, const string& attr) const
{
  return evalPair(resolveID(var, attr));
}

SymbolTable::AttrPair SymbolTable::resolveID(const string& s) const
{
  unsigned int i = 0;

  return evalString(s, i);
}

SymbolTable::AttrPair SymbolTable::resolveID
   (const string& s, const string& attr) const
{
  unsigned int i = 0;

  AttrPair a;

  top_level = false;                    // Don't worry about garbage at end

  a.first = getPrimary(s, i);

  a.second = getAttrID(attr, i = 0);

  top_level = true;

  if(!robust && i < s.size())            // Stuff at end problem again.
  {
    a.first  = NULL;
    a.second = (attr_id) -1;
  }

  return a;
}
 
string SymbolTable::resolve_string(const string& s) const
{
  string n;

  top_level = false;

  unsigned int j;
  AttrVal temp;

  for(unsigned int i = 0; i < s.size(); i++)
    switch(s[i])
    {
      case '\\' : n += s[++i];                                    break;
      case '$'  : j = i;
                  temp = evalPair(evalString(s,++i));
                  if(temp.Type() == AttrVal::Attr_None)
                    n += s.substr(j, i - j);
                  else
                    n += temp.String();
                  --i; break;
      default   : n += s[i];
    }

  top_level = true;

  return n;
}
  

// Private member functions


SymbolTable::AttrPair SymbolTable::evalString
   (const string& s, unsigned int& i, bool create, LSFBase* tmpl)
{
  AttrPair a;
  a.first  = NULL;
  a.second = (attr_id) -1;

  if(s.size() <= i) return a;

  if(s[i] == '$' && s[i+1] != '(')               // We have a $VAR expansion
  {
    a = evalString(s,++i);                       // Evaluate the VARIABLE

    unsigned int j = 0;
    return evalString(evalPair(a).String(), j); // Return eval of VARIABLE
  }

  bool tl = top_level;
  top_level = false;

  a.first = getPrimary(s, i, create, tmpl);

  if(s.size() <= i || s[i] != '.')      // No attribute specified
    a.second = 0;
  else
    a.second = getAttrID(s, ++i, create);

  if(!robust &&                               // If we're not robust and
     (top_level = tl) &&                      // we're at the top level and
     (i < s.size()) )                         // and there is more string
  {
    a.first  = NULL;                          // It isn't a valid identifier.
    a.second = (attr_id) -1;
  }

  return a;
}

LSFBase* SymbolTable::getPrimary
   (const string& s, unsigned int& i, bool create, LSFBase* tmpl)
{
  if(s.size() <= i) return NULL;

  string ident;
  LSF_ptr<LSFBase> p, a;                      // p - primary index, a - array index

  p = getPrimary(s, i, &ident);        // Find current value of identifier

  if(p)
  {
    a = getIndexObj(p, s, i, create, tmpl);

    if(!tmpl || p != a)                 // If we're not replacing or there
      return a;                         // we're not replacing the primary.
  }
  else
  {
    if(!create)
      return NULL;

    if(tmpl) p = Factory->build(tmpl->type()); // This may be temporary!

    if(!p)   p = new LSFBase();         // If previous doesn't help.

    p->name(ident.c_str());             // Set the name

    a = getIndexObj(p, s, i, create, tmpl);

    if(p == a && tmpl)                  // If temporary not used
    {
      a = tmpl;
    }
    else
    {
      tmpl = p;                         // use temporary for insertion
    }
  }

  if(!tmpl->name().size())
    tmpl->name(ident.c_str());

  st.set(ident, tmpl, true);

  return a;
}

LSFBase* SymbolTable::getIndexObj
    (LSFBase* current, const string& s, unsigned int& i, bool create,
     LSFBase* tmpl)
{
  if(s.size() <= i || s[i] != '[')  // i too big or no indicies
    return current;

  unsigned long index, j;

  LSF_ptr<LSFBase> prev = NULL;
  LSFList::iterator it;

  while(i < s.size() && s[i] == '[')
  {
    index = getNumber(s, ++i);

    if((signed) index < 0)                       // Problem getting index
      return NULL;

    // If current ever becomes NULL, we can't do anything, but
    // we keep parsing anyway, so that we can get the attribute ID
    if(current)                                  // If we still have current
    {
      if(create)                                 // create list and make sure
      {                                          // it holds enough elements.
        LSF_ptr<LSFBase> b;
        while(current->List(true)->size() <= index)
        {
          if(tmpl)
          {
            b = Factory->build( tmpl->type() );
            if(!b)
              b = new LSFBase();
          }
          else
            b = new LSFBase();

          current->List()->add(b);
        }
      }

      if(!current->List() ||
         current->List()->size() <= index)       // There's no object avail.
      {
        current = NULL;
        continue;
      }

      it = current->List()->begin();             // Find the next element
      for(j = 0; j < index; ++j, ++it);

      prev = current;                            // Store previous LSFBase*
      current = *it;                             // Make this index current
    }

    if(s[i] != ']')                              // Error... Do what we can
      if(robust)
        for(; i < s.size() && s[i] != ']'; i++);
      else
        return NULL;

    ++i;                                         // Go to next char.

  }

// NOTE:  This should be part of the LSFList interface!  Fix later
  if(tmpl && prev)                               // If there's a template
  {                                              // and this isn't top level
    *it = tmpl;
  }

  return *it;
}

attr_id  SymbolTable::getAttrID(const string& s, unsigned int& i,
                                bool create)
{
  string attr = getIdentifier(s, i);

  if(!attr.size()) return (attr_id) -1;

  if(!create)
    return AttrNameMgr.query(attr);

  return AttrNameMgr.add(attr);
}

LSFBase* SymbolTable::setPair(const string& s, AttrPair p,
                              const AttrVal& v, bool create)
{
  if(p.first)
  {
    if(p.first->attrs(create) && (signed) p.second != -1)
      p.first->attrs()->set(p.second, v);
  }
  else
  {
    if(!create || (signed) p.second == -1)             // Can't do much here.
      return NULL;

    AttrList al;
    al.set(p.second, v);
    p.first = Factory->build(NULL,NULL, &al);
    if(!p.first) return NULL;
    add(s, p.first);
  }
  return p.first;
}

SymbolTable::AttrPair SymbolTable::evalString
    (const string& s, unsigned int& i) const
{
  AttrPair a;
  a.first  = NULL;
  a.second = (attr_id) -1;

  if(s.size() <= i) return a;

  if(s[i] == '$' && s[i+1] != '(')              // We have a $VAR expansion
  {
    a = evalString(s,++i);                      // Evaluate the VARIABLE

    unsigned int j = 0;
    return evalString(evalPair(a).String(), j); // Return eval of VARIABLE
  }

  bool tl = top_level;
  top_level = false;

  a.first = getPrimary(s,i);

  if(s.size() <= i || s[i] != '.')      // No attribute specified
    a.second = 0;
  else
    a.second = getAttrID(s, ++i);

  if(!robust &&                               // If we're not robust and
     (top_level = tl) &&                      // we're at the top level and
     (i < s.size()) )                         // and there is more string
  {
    a.first  = NULL;                          // It isn't a valid identifier.
    a.second = (attr_id) -1;
  }

  return a;
}

LSFBase* SymbolTable::getPrimary
    (const string& s, unsigned int& i, string* id) const
{
  string ident = getIdentifier(s, i);

  if(!ident.size())                      // We didn't get a valid identifier
  {
    if(id) *id = "";
    return NULL;
  }

  LSFBase* p = st.find(ident);           // Find current value of identifier

  if(id)
    *id = ident;                         // return the identifier
  else
    p = getIndexObj(p, s, i);            // Find the right LSFBase
    
  return p;
}

LSFBase* SymbolTable::getIndexObj(LSFBase* current, 
                                  const string& s, unsigned int& i) const
{
  if(s.size() <= i || s[i] != '[')  // i too big or no indicies
    return current;

  unsigned long index, j;

  LSFList::iterator it;

  while(i < s.size() && s[i] == '[')
  {
    index = getNumber(s, ++i);

    if((signed) index < -1)                      // Problem getting index
      return NULL;

    // If current ever becomes NULL, we can't do anything, but
    // we keep parsing anyway, so that we can get the attribute ID
    if(!current || !current->List() ||
       current->List()->size() <= index)         // There's no object avail.
    {
      current = NULL;
    }
    else
    {
      it = current->List()->begin();
      for(j = 0; j < index; j++, it++);

      current = *it;                             // Make this index current
    }
    if(s[i] != ']')                              // Error... Do what we can
      if(robust)
        for(; i < s.size() && s[i++] != ']'; i++);
      else
        return NULL;

    i++;                                         // Go to next char.
  }

  return current;
}

attr_id  SymbolTable::getAttrID(const string& s, unsigned int& i) const
{
  string attr = getIdentifier(s, i);

  if(!attr.size()) return (attr_id) -1;

  return AttrNameMgr.query(attr);
}

string SymbolTable::getIdentifier(const string& s, unsigned int& i) const
{
  if(s.size() <= i) return string();

  if(s[i] == '$')                                  // Another expression
  {
    if(s[++i] != '(') return string();

    string ident =                                 // Get the parenthetical
        InsideParens(s, i).String();               // value

    if(!isalnum(ident[0])) return string();        // Test for validity

    unsigned int j = 0;
    for(; j < ident.size() && 
          (strchr("-_",ident[j]) || isalnum(ident[j])); j++);

    return ident.substr(0, j);
  }    
  
  unsigned int start = i;

  if(!isalnum(s[i++])) return string();

  for(; i < s.size() && (strchr("-_",s[i]) || isalnum(s[i])); i++);

  return s.substr(start, i-start);
}

unsigned int SymbolTable::getNumber(const string& s, unsigned int& i) const
{
  if(s.size() <= i) return 0;

  if(s[i] == '$')                                 // Another expression
  {
    AttrVal v = evalPair(evalString(s, ++i));
    
    if(!robust)
      switch(v.Type())
      {
        case AttrVal::Attr_Real    :          // Check for valid value
          if(v.Real() != v.Int())
            return (unsigned int) -1;

        case AttrVal::Attr_Integer :
          break;
      
        case AttrVal::Attr_String  :
          if(AttrVal(v.Int()).String() == v.String())
          {
            v = v.Int();
            break;
          }

        default                    :
          return (unsigned int) -1;
      }
    return v.Int();
  }    
  
  if(isdigit(s[i]))
  {
    unsigned int index = atoi(s.c_str() + i);
    if((signed) index < 0) index = 0;            // No negative values!

    for( ; i < s.size() && isdigit(s[i]); i++);  // Read over number
    return index;
  }

  if(robust)
    return 0;
  else
    return (unsigned int) -1;
}

AttrVal SymbolTable::InsideParens(const string& s, unsigned int& i) const
{
  AttrPair a = evalString(s, ++i);

  if(!a.first || (signed) a.second == -1)       // Bad stuff.
    return AttrVal();

  if(s[i++] != ')')                             // Problem.  Do best we can.
    if(robust)
    {
      for( ; i < s.size() && s[i] != ')'; i++);
      i++;
    }
    else
      return AttrVal();

  return evalPair(a);
}

AttrVal SymbolTable::evalPair(const AttrPair& a) const
{
  AttrList*          al;

  if(a.first && (al = a.first->attrs()))
  {
    AttrList::iterator it = al->find(a.second);
    if(it != al->end() )
      return (*it).second;
  }

  return AttrVal();
}

