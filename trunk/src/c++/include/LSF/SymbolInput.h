#ifndef _SYMBOL_INPUT_H
#define _SYMBOL_INPUT_H

#include "LSF/LSF.h"
#include "LSF/LSFsymbol.h"

#ifdef NAMESPACE_STD
NAMESPACE_STD
#endif

class SymbolInput
{
public:

  enum result { NO_SYM_TAB, NOT_IN_ST, OUT_OF_BOUND, NO_DEFAULT, REQUIRED,
                CONV_ERROR, INVALID, VALID };

  SymbolInput(SymbolTable* = NULL);
  
  ~SymbolInput();
  
  result set(const string& symbol, const string&  attr);
  result set(const string& symbol, const AttrVal& attr);

  // For dumb compilers
  result set(const char *symbol, const AttrVal& attr) { return set( string(symbol), attr ); }
  result set(const char *symbol, const char *attr)    { return set( string(symbol), string(attr) ); }
  
  // Do not use internal data members so may be set static.
  static result set(LSFBase* b, const attr_id& aid, const string& attr);
  static result set(LSFBase* b, const attr_id& aid, const AttrVal& attr);

protected:

  // Do not use internal data so may be set static.
  static bool MapString(LSFBase*, const attr_id&, string&);

  // tests for REQUIRED attribute set 
  static result test_required(LSFBase*);

  SymbolTable* sym;

};

#endif
