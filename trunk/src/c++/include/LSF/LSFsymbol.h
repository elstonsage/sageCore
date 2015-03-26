#ifndef __SYMBOL_H
#define __SYMBOL_H
//
//  Symbol Table 0.1 -- Symbol Tables using LSF
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   0.0  gcw Initial implementation.                 Oct 17 1996
//             1.0  gcw Redesigned for multi dimensional arrays Jul  2 1997
//       
//  Copyright (c) 1996  R.C. Elston
//


#include <string>
#include <utility>
#include <ctype.h>
#include "LSF/LSF.h"
#include "LSF/LSFmap.h"

// WARNING:  THIS NUMBER SHOULD BE CHANGED!!
#define SYMBOL_TABLE 10010

/**
 *  The SymbolTable is an LSFBase derived class for storing data symbolically.
 *  In essence, LSFBase objects, and attributes set on those objects are
 *  stored within the symbol table in a way in which they can be referred
 *  to efficiently.
 *
 *  <P> There are essentially three sections to the interface to the
 *  SymbolTable.  The first is derived from LSFmap.  The SymbolTable uses
 *  LSFmap for its main storage.  This part of the interface provides
 *  iterators based on the primary index name (see indexing).  It is mainly
 *  used for output of the complete symbol table.
 *
 *  <P> The second section of the interface is used for storing and
 *  retrieving LSFBase pointers in the SymbolTable (see indexing). 
 *  Attribute identifiers should not be used here.
 *
 *  <P> The third section of the interface stores and retrieves variables
 *  (strings, integers and doubles) as AttrVals in the SymbolTable. 
 *  Attribute Identifiers are used for this.  Missing Attribute names are
 *  understood to be VALUE or attr_id 0.
 *
 *  <P><B>Indexing</B>
 *
 *  <P> An index into the symbol table consists of the following parts:
 *  <UL>
 *  <LI>A primary index.  This is an identifier (see below)  This is the
 *      index into the LSFmap inside the SymbolTable.
 *  <LI>(optional) One or more array indices of the form [#][#]...  These
 *      index the primary object identifier's LSFList.  These indices should
 *      resolve to a integer >= 0.
 *  <LI>(optional) A '.' followed by a attribute identifier (see below)
 *  </UL>
 *  <P> An <B>Identifier</B> is a  string, consisting of at least one
 *      alpha-numeric character followed by a combination of alphanumeric
 *      characters or '-' or '_' (see grammar below).
 *  <P> Additionally, a value (primary index, # inside brackets or attribute
 *      identifier) may be itself a reference into the symbol table.  This
 *      is useful for array indicies, etc.  This is done by enclosing the
 *      section in parentheses, or prefacing it by a '$'.  As an example, an
 *      integer i is stored in the SymbolTable.  We may use this to index an
 *      array a as a[$i] or a[$(i)].
 *  <P> The following is a very left non-recursive CFG in EBNF for LSF
 *      variables.  It can be converted into LALR or LL in this form, and
 *      with the addition of semantic information most ambiguities disappear
 *      in actual usage.  The only potential problem is a shift/reduce
 *      conflict due to the left associativity of ']' in a LR parser scheme. 
 *      Also, it is assumed that the syntactic pre-processing has correctly
 *      tokenized the input.  A more complete CFG would include lexical
 *      parsing information as well.
 *  <P> Things contained in {}s are optional.  ()* means zero or more of
 *      things contained in the ().  ()+ means at least one.  | indicates
 *      or.  A \ preceding a ( or ) indicates that parenthesis is a
 *      character, not part of the expression.
 *
 *  <PRE>
 *    VARIABLE ::= VAREXP ([ $VARIABLE | IDENTIFIER ])* {.VARIABLE}
 *               | $VARIABLE 
 *    
 *    // Protected variable expression (protected by order of operation)
 *    VAREXP   ::= $\(VARIABLE\)
 *               | IDENTIFIER
 *    
 *    IDENTIFIER = ALPHANUM(ALPHANUM | '-' | '_')*
 *    
 *    ALPHANUM = alphanumeric character.
 *  </PRE>    
 *
 *  <P> Currently as implemented, primary identifiers are case sensitive.
 *      Attribute identifiers are not case sensitive.
 *  <P><B>Examples of valid identifiers:</B>
 *  <UL>
 *  <LI>F
 *  <LI>$F
 *  <LI>$(F)
 *  <LI>F.a
 *  <LI>F[5]
 *  <LI>F[5].a
 *  <LI>F[$i].$(a) (provided i, a are defined in the symbol table)
 *  <LI>F.$a.i (provided a.i is defined in the symbol table)
 *  <LI>$(F.a).i (provided F.a is defined in the symbol table)
 *  <LI>F[(i)].(a) (provided i, a are defined in the symbol table)
 *  <LI>F[$G.x][5][$i].$N[1] (This is the attribute given as the value of
 *      N[1] on the LSFobject located by finding the G.xth element of F, the
 *      5th element of that, and the ith element of that.  Not that we'd
 *      want something this complex, of course.
 *  </UL>
 *  
 */
class SymbolTable : public LSFBase
{
public:
  SymbolTable(const char *n = "") : LSFBase(n), robust(false), top_level(true)
  { _local_type = SYMBOL_TABLE; }

  SymbolTable(const SymbolTable &s) : st(s.st), top_level(true)
  {
    name( s.name() );
    _local_type = SYMBOL_TABLE;
    set_robust( s.get_robust() );
  }
  
  virtual ~SymbolTable()  { }

  typedef LSFmap::key_type               key_type;
  typedef LSFmap::key_compare            key_compare;
  typedef LSFmap::reference              reference;
  typedef LSFmap::const_reference        const_reference;
  typedef LSFmap::iterator               iterator;
  typedef LSFmap::const_iterator         const_iterator;
  typedef LSFmap::reverse_iterator       reverse_iterator;
  typedef LSFmap::const_reverse_iterator const_reverse_iterator;
  typedef LSFmap::size_type              size_type;
  typedef LSFmap::difference_type        difference_type;
  typedef pair<LSFBase*, attr_id>        AttrPair;

  virtual LSFList       *List(bool create = FALSE) 
  { return st.List(create); }
  virtual const LSFList *List(bool create = FALSE) const  
  { return st.List(create); }
  virtual void Accept(LSFVisitor *v) { v->Process(this); }
  virtual void AcceptAll(LSFVisitor *v) 
  { for(iterator i=begin(); i!=end(); i++) (*i).second->Accept(v); }

// LSFmap accessors:

  key_compare          key_comp() const  { return st.key_comp(); }
  iterator                begin()        { return st.begin();    }
  const_iterator          begin() const  { return st.begin();    }
  iterator                  end()        { return st.end();      }
  const_iterator            end() const  { return st.end();      }
  reverse_iterator       rbegin()        { return st.rbegin();   }
  const_reverse_iterator rbegin() const  { return st.rbegin();   }
  reverse_iterator         rend()        { return st.rend();     }
  const_reverse_iterator   rend() const  { return st.rend();     }
  bool                    empty() const  { return st.empty();    }
  size_type                size() const  { return st.size();     }

// Robustness

   /* sets robust flag as specified. The robustness flag determines how
    * the SymbolTable handles badly formed strings. */
   //@{
    bool get_robust() const        { return robust;     }
    bool set_robust(bool b= false) { return robust = b; }
   //@}

// LSFBase based indexing

  /// Add a LSFBase*, create it if b == NULL.
  LSFBase*        add(const string& s, LSFBase* b = NULL) 
  { return replace(s, b, true); }
  /// Replace the symbol with a new one, creating on NULL and create == true.
  LSFBase*    replace(const string&, LSFBase*, bool create = false);
  /// Removes a symbol (and all its indices and attributes) from the table
  void         remove(const string& s);
  LSFBase*       find(const string& s) const;
  LSFBase* operator[](const string &s         ) const { return find(s); }
  
// AttrVal based indexing.
  
  LSFBase* add(const string& s, const AttrVal& v)
      { return set(s, v, true); }
  LSFBase* add(const string& s, const string& a, const AttrVal& v) 
      { return set(s, a, v, true); }

  /// Set a symbol.  If bool = true, symbol will be created if necessary.
  //@{
  LSFBase* set(const string&, const AttrVal&, bool create = false);
  LSFBase* set(const string& var, const string& attr,
               const AttrVal&, bool create = false);
  //@}
  
  /** Find an AttrVal if it exists.
   *  @return Value if it exists, AttrVal() if it does not. */
  //@{
  AttrVal query(const string&)                              const;
  AttrVal query(const string& var, const string& attr)      const;
  //@}

  /// Return LSFBase*, attr_id pair.
  //@{
  AttrPair resolveID(const string&)                         const;
  AttrPair resolveID(const string& var, const string& attr) const;
  //@}
  
  /// Parse a string.  '$' indicates a symbol to be parsed.
  string  resolve_string(const string& s) const;
  
private:

  LSFmap st;

  /** The robust flag determines behavior when symbols are 
   *  badly formed.  It is normally false requiring symbols to follow
   *  the format strictly. */
  bool   robust;                     

  /// The top_level flag keeps track of where we are in the parse tree.
  mutable bool   top_level;

  /** Non-const functions for adding and modifying the symbol table.
   *  The create flag determines if they create new symbols in the symbol
   *  table.  */
  //@{
  AttrPair evalString   (const string& s, unsigned int& i,
                         bool create, LSFBase* tmpl);
  LSFBase* getPrimary   (const string& s, unsigned int& i,
                         bool create, LSFBase* tmpl);
  LSFBase* getIndexObj  (LSFBase* current, const string& s, unsigned int& i,
                         bool create, LSFBase* tmpl);
  attr_id  getAttrID    (const string& s, unsigned int& i,
                         bool create);

  /** Sets pair's values.  If Pair invalid (NULL) and create true, it tries
   *  to make a vlaid identifier. */
  LSFBase* setPair      (const string& s, AttrPair p, const AttrVal& v,
                         bool create);
  //@}
  
  /// Const functions for lookup
  //@{
  AttrPair evalString   (const string& s, unsigned int& i) const;
  LSFBase* getPrimary   (const string& s, unsigned int& i, 
                         string* id = NULL) const;
  LSFBase* getIndexObj  (LSFBase* current, const string& s,
                         unsigned int& i) const;
  attr_id  getAttrID    (const string& s, unsigned int& i) const;
  
  /// Gets an identifier string - does not assume interpretation
  string  getIdentifier (const string& s, unsigned int& i) const;

  unsigned int getNumber(const string& s, unsigned int& i) const;

  AttrVal InsideParens  (const string& s, unsigned int& i) const;
  AttrVal evalPair      (const AttrPair& a)                const;
  
  //@}
  
};

#endif

