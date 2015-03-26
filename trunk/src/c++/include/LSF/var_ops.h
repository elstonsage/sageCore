#ifndef _VAR_OPS_H
#define _VAR_OPS_H

#include <string>
#include "LSF/LSF.h"
#include "LSF/LSFsymbol.h"

using namespace std;

void addIndexSymbols( const LSFBase *root, SymbolTable *syms);
bool setDef(LSFBase *b, bool clear = FALSE);
void expandArray(LSFBase *root, LSFBase *tmpl, int size);
LSFBase *extractHead(LSFBase *b);

// Add a single variable to a symbol table
void addSymbol(LSFBase *v, SymbolTable *vars);

AttrVal resolveAttrVal(const AttrVal &v, const SymbolTable* sym_tab);
AttrVal find_and_resolveAttrVal(const string &name, const SymbolTable *sym_tab);

#endif
