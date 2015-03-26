#ifndef __NODE_H
#define __NODE_H

//#include <utility>
#include <list>

#ifdef __TESTING_NODES

#define allele int

#else

#include "fped/fped.h"
#include "mlocus/imodel.h"

#define allele SAGE::MLOCUS::allele

#endif

namespace SAGE
{

class edge;

class node
{
public:

  friend class edge;

  typedef int id_type;

  node();
  node(id_type);

  id_type id() const;
  
  bool state() const;
  bool fixed() const;
  
  bool fix(int);
  void unfix();
  
  bool is_sex_null() const;
  void set_sex_null_state(bool);

#ifndef __TESTING_NODES
  double probability(const SAGE::MLOCUS::inheritance_model&, bool);
#endif

  void dump(std::ostream&, int);

  bool add_edge(edge&, int edge_end);

  void remove_edge();

  bool operator==(const node&) const;

protected:

  bool fix(edge&, allele);
  bool propagate(edge&);
  void undo(edge&);
  
  void save(int);
  void restore();

  typedef std::list<edge*>    edge_list;
  typedef edge_list::iterator iterator;

  edge_list edges;

  int ch;

  allele a[2];
  allele f_allele;

  id_type my_id;

  bool my_is_sex_null;
};

class edge
{
public:

  struct genotype
  {
    allele a1;
    allele a2;
    
    allele allele1() const { return a1; }
    allele allele2() const { return a2; }
    
    void flip() { std::swap(a1,a2); } 
 
  };

  friend class node;

#ifndef __TESTING_NODES
  typedef FPED::Multipedigree::member_const_pointer id_type;
#else
  typedef int id_type;
#endif
  
  enum state_type { homozygous, sorted, unsorted, phased };

  edge();
  edge(id_type, genotype, bool phased = false);
  edge(id_type, std::pair<allele, allele>, bool phased = false);

  
  id_type id() const;

  state_type state() const;

  allele get_allele(int);

  void dump(std::ostream&);

protected:


  bool propagate(node&);
  bool fixed_propagate(node&);
  void undo(node&, bool = false);

  id_type             my_id;
  state_type          my_state;
  genotype   my_genotype;

  node* nodes[2];

  int ch;

  static int connection;
};

// ================
// Inline Functions
// ================

inline node::node()
{ f_allele = a[0] = a[1] = allele(); ch = -1; my_id = 0; my_is_sex_null = false; }

inline node::node(id_type i)
{ f_allele = a[0] = a[1] = allele(); ch = -1; my_id = i; my_is_sex_null = false; }

inline node::id_type node::id() const { return my_id; }

inline bool node::state() const
{ return !edges.empty(); }

inline bool node::fixed() const
{ return f_allele != allele(); }

inline bool node::is_sex_null() const
{
  return my_is_sex_null;
}

inline void node::set_sex_null_state(bool s)
{
  my_is_sex_null = s;
}

inline void node::save(int n)
{
  f_allele = a[n];
}

inline void node::restore()
{
  f_allele = allele();
}

inline bool node::operator==(const node& n) const
{ return fixed() && n.fixed() && f_allele == n.f_allele; }

inline edge::edge()
{
  nodes[0] = nodes[1] = NULL;
  ch = -1;
}

inline edge::edge(id_type i, genotype g, bool ph)
  : my_id      (i),
    my_genotype(g)
{
  nodes[0] = nodes[1] = NULL;
  ch = -1;

  if(my_genotype.allele1() == my_genotype.allele2())
  {
    my_state = homozygous;
  }

  else if(ph)
  {
    my_state = phased;
  }
  else
  {
    my_state = unsorted;
  }
}
inline edge::edge(id_type i, std::pair<allele, allele> g, bool ph)
  : my_id(i)
{
  my_genotype.a1 = g.first;
  my_genotype.a2 = g.second;

  nodes[0] = nodes[1] = NULL;
  ch = -1;

  if(my_genotype.allele1() == my_genotype.allele2())
  {
    my_state = homozygous;
  }
  else if(ph)
  {
    my_state = phased;
  }
  else
  {
    my_state = unsorted;
  }
}

inline edge::id_type edge::id() const
{ return my_id; }

inline edge::state_type edge::state() const { return my_state; }

inline allele edge::get_allele(int i)
{
  if(i == 0) return my_genotype.allele1();
  else       return my_genotype.allele2();
}

}

#endif
