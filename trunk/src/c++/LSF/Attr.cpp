//
//  ATTRIBUTE LISTS 1.0 -- Low level attribute list object
//  
//  Purpose: Implements a list of arbitrary attributes of fixed
//           symbol and value types.
//  
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation	Mar 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//

#include <utility>
#include "globals/config.h"
#include "LSF/NameMgr.h"
#include "LSF/Attr.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

AttrManager AttrNameMgr;

bool AttrVal::operator>= (const AttrVal &rhs) const
{ 
  if(rhs.attr_type == Attr_None    && attr_type == Attr_None) return true;
  if(rhs.attr_type == Attr_None    || attr_type == Attr_None) return false;
  return !(*this < rhs);
}

bool AttrVal::operator<= (const AttrVal &rhs) const
{ 
  if(rhs.attr_type == Attr_None    && attr_type == Attr_None) return true;
  if(rhs.attr_type == Attr_None    || attr_type == Attr_None) return false;
  return !(rhs < *this);
}

bool AttrVal::operator> (const AttrVal &rhs) const
{ 
  if(rhs.attr_type == Attr_None    || attr_type == Attr_None) return false;
  return (rhs < *this);
}

bool AttrVal::operator< (const AttrVal &rhs) const
{ 
  if(rhs.attr_type == Attr_None    || attr_type == Attr_None) return false;
  if(rhs.attr_type == Attr_String  || attr_type == Attr_String)
    return ( rhs.String() > String() );
  if(rhs.attr_type == Attr_Real    || attr_type == Attr_Real)
    return ( rhs.Real() > Real() );
  if(rhs.attr_type == Attr_Ptr     || attr_type == Attr_Ptr)
    return ( rhs.Ptr() > Ptr() );
  if(rhs.attr_type == Attr_Integer || attr_type == Attr_Integer)
    return ( rhs.Int() > Int() );
  return false;
}
    
bool AttrVal::operator==(const AttrVal &rhs) const
{ 
  if(rhs.attr_type == Attr_None    && attr_type == Attr_None) return true;
  if(rhs.attr_type == Attr_None    || attr_type == Attr_None) return false;
  if(rhs.attr_type == Attr_String  || attr_type == Attr_String)
    return ( rhs.String() == String() );
  if(rhs.attr_type == Attr_Real    || attr_type == Attr_Real)
    return ( rhs.Real() == Real() );
  if(rhs.attr_type == Attr_Ptr     || attr_type == Attr_Ptr)
    return ( rhs.Ptr() == Ptr() );
  if(rhs.attr_type == Attr_Integer || attr_type == Attr_Integer)
    return ( rhs.Int() == Int() );
  return false;
}

AttrVal AttrVal::operator-() const
{
  if(attr_type==Attr_Real)    return AttrVal( -Real() );
  if(attr_type==Attr_Integer) return AttrVal( -Int()  );
  return AttrVal();
}

AttrVal AttrVal::operator+(const AttrVal &rhs) const
{
  if(rhs.attr_type==Attr_None   || attr_type==Attr_None)   return AttrVal();
  if(rhs.attr_type==Attr_Ptr    || attr_type==Attr_Ptr)    return AttrVal();
  if(rhs.attr_type==Attr_String || attr_type==Attr_String)
    return String() + rhs.String();

  if(rhs.attr_type == Attr_Real   || attr_type == Attr_Real)
    return ( rhs.Real() + Real() );
  if(rhs.attr_type == Attr_Integer || attr_type == Attr_Integer)
    return ( rhs.Int() + Int() );
  return AttrVal();
}
  
AttrVal AttrVal::operator-(const AttrVal &rhs) const
{
  if(rhs.attr_type==Attr_None   || attr_type==Attr_None)   return AttrVal();
  if(rhs.attr_type==Attr_String || attr_type==Attr_String) return AttrVal();
  if(rhs.attr_type==Attr_Ptr    || attr_type==Attr_Ptr)    return AttrVal();

  if(rhs.attr_type == Attr_Real   || attr_type == Attr_Real)
    return ( Real() - rhs.Real() );
  if(rhs.attr_type == Attr_Integer || attr_type == Attr_Integer)
    return ( Int() - rhs.Int() );
  return AttrVal();
}

AttrVal AttrVal::operator*(const AttrVal &rhs) const
{
  if(rhs.attr_type==Attr_None  || attr_type==Attr_None)   return AttrVal();
  if(rhs.attr_type==Attr_String|| attr_type==Attr_String) return AttrVal();
  if(rhs.attr_type==Attr_Ptr   || attr_type==Attr_Ptr)    return AttrVal();

  if(rhs.attr_type==Attr_Real  || attr_type==Attr_Real)
    return ( rhs.Real() * Real() );
  if(rhs.attr_type==Attr_Integer || attr_type==Attr_Integer)
    return ( rhs.Int() * Int() );
  return AttrVal();
}

AttrVal AttrVal::operator/(const AttrVal &rhs) const
{
  if(rhs.attr_type==Attr_None  || attr_type==Attr_None)   return AttrVal();
  if(rhs.attr_type==Attr_String|| attr_type==Attr_String) return AttrVal();
  if(rhs.attr_type==Attr_Ptr   || attr_type==Attr_Ptr)    return AttrVal();

  if(rhs.attr_type==Attr_Real  || attr_type==Attr_Real)
  {
    if(Real() != 0.0)
      return ( Real() / rhs.Real() );
    else
      return AttrVal();
  }

  if(rhs.attr_type==Attr_Integer || attr_type==Attr_Integer && Int() != 0)
    return (Int() / rhs.Int() );
  
  return AttrVal();
}

AttrVal AttrVal::operator%(const AttrVal &rhs) const
{
  if(rhs.attr_type==Attr_None  || attr_type==Attr_None)   return AttrVal();
  if(rhs.attr_type==Attr_String|| attr_type==Attr_String) return AttrVal();
  if(rhs.attr_type==Attr_Ptr   || attr_type==Attr_Ptr)    return AttrVal();
  if(rhs.attr_type==Attr_Real  || attr_type==Attr_Real)   return AttrVal();
  
  if(rhs.attr_type == Attr_Integer || attr_type == Attr_Integer && Int() != 0)
    return ( Int() % rhs.Int() );
  
  return AttrVal();
}
    
int AttrVal::Int() const 
{
  switch( attr_type )
  {
    case Attr_Integer: return attr_val.i;
    case Attr_Real:    return (int)attr_val.f;
    case Attr_String:  return str2long(*attr_val.s);
    case Attr_Ptr:     return (long)attr_val.p;
    case Attr_None:    return 0;
  }
  return 0;
}
    
double AttrVal::Real() const 
{
  switch( attr_type )
  {
    case Attr_Integer: return attr_val.i;
    case Attr_Ptr:
    case Attr_Real:    return attr_val.f;
    case Attr_String:  return str2doub(*attr_val.s);
    case Attr_None:    return 0;
  }
  return 0;
}

string AttrVal::String(int width, int pres, long flags, char fch ) const
{
  string s;
  switch( attr_type )
  {
    case Attr_Integer: s = long2str(attr_val.i, width, flags, fch);
                       break;
    case Attr_Real:    s = doub2str(attr_val.f, width, pres, flags, fch);
                       break;
    case Attr_String:  return *attr_val.s;
    case Attr_Ptr:     s = ptr2str(attr_val.p, width, flags, fch);
    case Attr_None:    break;
  }
  return s;
}

const void *AttrVal::Ptr() const 
{
  switch( attr_type )
  {
    case Attr_Integer: return (void *)attr_val.i;
    case Attr_Real:    return (void *)((int)attr_val.f);
    case Attr_String:  return (void *)attr_val.s->c_str();
    case Attr_Ptr:     return attr_val.p;
    case Attr_None:    break;
  }
  return NULL;
}
