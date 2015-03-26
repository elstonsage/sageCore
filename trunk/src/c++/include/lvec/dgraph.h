#ifndef __DESCENT_GRAPH_H
#define __DESCENT_GRAPH_H

//
//  Descent Graphs: Representations of allele sharing from inheritance vectors
//
//  Author: Geoff Wedig
//
//  History:  0.1 Initial Implementation      May  11 1998
//
//  Copyright (c) 1998 R. C. Elston

#include "lvec/inheritance_vector.h"

namespace SAGE
{

class descent_graph
{
public:

  typedef inheritance_vector::equivalence_class equivalence_class;
  typedef long                                  allele_type;
  typedef pair<allele_type, allele_type>        allele_pair;
  typedef meiosis_map::index                    index;
  typedef meiosis_map::member_const_pointer     member_const_pointer;

  descent_graph(const inheritance_vector&, bool a = true);
  descent_graph(equivalence_class, const meiosis_map&, bool a = true);
  
  ~descent_graph();
  
  void move_to(const inheritance_vector&);
  void move_to(equivalence_class, const meiosis_map&);

// return types

  allele_pair alleles(member_const_pointer) const;
  allele_pair alleles(index)  const;

  // Sharing assumes no loops
  int sharing(member_const_pointer, member_const_pointer) const;
  int sharing(index,  index)  const;

  int mother_sharing(member_const_pointer, member_const_pointer) const;
  int mother_sharing(index,  index)  const;

  int father_sharing(member_const_pointer, member_const_pointer) const;
  int father_sharing(index,  index)  const;

  void print_allele_vector() const;

protected:

  typedef vector<allele_pair> avector;

  const meiosis_map* mm;
  avector av;  

  bool autosomal;
};

#include "lvec/dgraph.ipp"

}

#endif
