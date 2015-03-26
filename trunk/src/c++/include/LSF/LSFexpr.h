#ifndef __LSF_GUARD_H
#define __LSF_GUARD_H
//
//
// LSFexpr - LSF Expressions
//    includes a Guard node that protects it's members.
//    and a simple expression interpreter
//
// Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
// History: 0.1 gcw Initial Implementation       Nov 20 1997
//
//
// Copyright (c) R. C. Elston.


#include "LSF/LSF.h"
#include "LSF/LSFtypes.h"
#include "LSF/LSFsymbol.h"

// The LSFguard evaluator provides limitations to the access ofguard nodes.
// A guard node is an LSFbase object with its lsf_type() set to LSF_GUARD.
// The evaluator function limits access to the LSFList (if any) contained
// within the guard.  This is accomplished by evaluating the guard
// (optionally, with a SymbolTable). Should the guard's test (see below)
// pass, the LSFList* is returned, otherwise a NULL is returned.  If the
// evaluator should be called on a non-LSF_GUARD type object, the LSFList*
// of the object is returned.
//
// The test function is controlled by the attributes of hte guard node.
// There are three attributes to be aware of:
//
// VAR1 - First variable
// VAR2 - Second variable
// TEST - This is a string set to any of the following or their symbolic
//        equivalents: EQUALS, NOT_EQUALS, LESS_THAN, GREATER_THAN,
//                     LESS_EQUALS (less than or equal to), GREATER_EQUALS
//
// Either one or both of the variables may be Symbolic.  Symbolic values are
// indicated by the $'s as in parsed strings in a SymbolTable.  If a
// SymbolTable* is passed to the evaluator funtion, these Symbolic names are
// parsed.  If there is no SymbolTable, the comparison is done 'as is' in
// which case symbolic variable names are taken as strings to be compared. 
// If the variable is not in the SymbolTable, the base AttrVal() is used in
// the comparison (Note: NOT_EQUALS has a problem in either of these cases. 
// It is recommended to use some other test where possible to avoid this.)

LSFList* evaluate_guard(LSFBase* b, SymbolTable* sym_tab);
LSFList* evaluate_switch(LSFBase* b, SymbolTable* sym_tab);
LSFList* evaluate_expr(LSFBase* b, SymbolTable* sym_tab);

bool is_lsfcond(LSFBase* b);
bool is_lsfexpr(LSFBase* b);

LSFList* evaluate_lsfexpr(LSFBase*, SymbolTable* = NULL);
LSFList* evaluate_lsfcond(LSFBase*, SymbolTable* = NULL);

#endif

