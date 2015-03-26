#ifndef PEN_FUNCTION_H
#define PEN_FUNCTION_H

#include <list>
#include <map>
#include <string>
#include <functional>
#include <iostream>
#include "mped/mp_utilities.h"
#include "fped/fped.h"

namespace SAGE
{
namespace SEGREG
{

/// pg is a struct defining the class for posterior genotype functions for later output.

/// The pg struct never needs instatiated.  It is purely for defines and
/// types regarding the production of posterior genotypes and penetrances
/// for output.  These values are individual specific, and individuals are
/// sorted based upon pedigree name and individual name within a map.
///
/// Values are stored in a list of simple structs, each containing a
/// genotype (string), a posterior probability for that genotype and a
/// penetrance for that genotype.  
struct pg
{
  typedef FPED::MemberConstPointer member_pointer;

  struct member_compare : public std::binary_function<member_pointer, member_pointer, bool>
  {
    bool operator() (member_pointer lhs,member_pointer rhs) const
    {
      const string& lhsped = lhs->pedigree()->name();
      const string& rhsped = rhs->pedigree()->name();
      
      return lhsped < rhsped ||
             (lhsped == rhsped && lhs->name() < rhs->name());
    }
  };


  struct pen_type
  {
    /// Constructors
    //@{

    pen_type()
      : geno      ("" ),
        post_geno (0.0),
        pen       (0.0)
    { }

    pen_type(const string& g, double pg, double pn)
      : geno      (g ),
        post_geno (pg),
        pen       (pn)
    { }

    pen_type(const pen_type& p)
      : geno      (p.geno     ),
        post_geno (p.post_geno),
        pen       (p.pen      )
    { }

    //@}
        
    /// Copy Operator
    
    pen_type& operator=(const pen_type& p)
    {
      geno      = p.geno;
      post_geno = p.post_geno;
      pen       = p.pen;

      return *this;
    }

    string geno;
    double post_geno; //< Penetrance of individual based upon genotype(s)
    double pen;       //< Penetrance of individual based upon genotype(s) and pedigree data
  };

  typedef list<pen_type>                                       member_pen_list;
  typedef std::map<member_pointer, member_pen_list, member_compare> post_geno_map;

  static void dump_post_genos(ostream&, const post_geno_map&);
};

/// pf is a struct defining the class for 3 person penetrance functions for later output.

/// The pf struct never needs instatiated.  It is purely for defines and
/// types regarding the production of penetrances for output.  These values
/// are individual specific, and individuals are sorted based upon pedigree
/// name and individual name within a map.
///
/// Values are stored in a list of simple structs, each containing a
/// genotype (string) for each of the individual, mother and father, and a
/// penetrance for that triple.
struct pf
{
  typedef FPED::MemberConstPointer member_pointer;

  struct member_compare : public std::binary_function<member_pointer, member_pointer, bool>
  {
    bool operator() (member_pointer lhs,member_pointer rhs) const
    {
      const string& lhsped = lhs->pedigree()->name();
      const string& rhsped = rhs->pedigree()->name();
      
      return lhsped < rhsped ||
             (lhsped == rhsped && lhs->name() < rhs->name());
    }
  };


  struct pen_type
  {
    /// Constructors
    //@{

    pen_type()
      : geno      ("" ),
        mo_geno   ("" ),
        fa_geno   ("" ),
        pen       (0.0)
    { }

    pen_type(const string& g, const string& mth, const string& fth, double pn)
      : geno      (g  ),
        mo_geno   (mth),
        fa_geno   (fth),
        pen       (pn )
    { }

    pen_type(const string& g, double pn)
      : geno      (g  ),
        mo_geno   ("" ),
        fa_geno   ("" ),
        pen       (pn )
    { }

    pen_type(const pen_type& p)
      : geno      (p.geno   ),
        mo_geno   (p.mo_geno),
        fa_geno   (p.fa_geno),
        pen       (p.pen    )
    { }

    //@}
        
    /// Copy Operator
    
    pen_type& operator=(const pen_type& p)
    {
      geno      = p.geno;
      mo_geno   = p.mo_geno;
      fa_geno   = p.fa_geno;
      pen       = p.pen;

      return *this;
    }

    string geno;
    string mo_geno;
    string fa_geno;
    double pen;       //< Penetrance of individual based upon genotype(s) and pedigree data
  };

  typedef list<pen_type>                                       member_pen_list;
  typedef std::map<member_pointer, member_pen_list, member_compare> pen_func_map;

  static void dump_pen_func(ostream&, const pen_func_map&);
};

}}

#endif
