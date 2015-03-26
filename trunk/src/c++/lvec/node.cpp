#include "lvec/node.h"

namespace SAGE
{

int edge::connection = 0;

#ifdef __TESTING_NODES
inline int dump_allele(allele i)
{
  return i;
}
#else
inline string dump_allele(allele i)
{
  if( !i.is_valid())
    return "null";

  return i.name();
}
#endif


void node::dump(std::ostream& o, int i)
{
  o << i << '\t' << state() << '\t' << fixed() << '\t' << std::flush;
#if 0
  o << i << '\t';
  if( state() )
    o << "edge(s)";
  else
    o << "no edge";
  o << '\t' << fixed() << '\t' << std::flush;
#endif

  if(ch != -1) o << ch << '\t';
  else         o << '\t';

  o << dump_allele(a[0]) << '/' << dump_allele(a[1])<< '/' 
    << dump_allele(f_allele) << '\t'
    << edge::connection << ' ' << std::flush;

  for(iterator j = edges.begin(); j != edges.end(); ++j)
  {
  #ifdef __TESTING_NODES
    o << (*j)->id();
  #else
    o << (*j)->id()->name() << " ";
  #endif
  }
  o << std::endl;
}

void edge::dump(std::ostream& o)
{
#ifdef __TESTING_NODES
  o << id();
#else
  o << id()->name();
#endif

  o << '\t' << state() << '\t' 
#if 0
  o << '\t';
  if( state() == edge::homozygous )
    o << "homo";
  else if( state() == edge::sorted )
    o << "sort";
  else if( state() == edge::unsorted )
    o << "unso";
  else if( state() == edge::phased )
    o << "phas";
  else
    o << "????";

  o << '\t' 
#endif
    << dump_allele(get_allele(0)) << '/'
    << dump_allele(get_allele(1)) << '\t';

  if(ch != -1) o << ch << '\t'; else o << '\t';

  o << (((int) (bool) nodes[0]) + ((int) (bool) nodes[1])) << '\t';

  if(nodes[1]) o << nodes[0]->id() << ' ' << nodes[1]->id();

  o << std::endl;
}

bool node::fix(int i)
{
  if(!state() || fixed()) return false;

  save(i);

  for(iterator i = edges.begin(); i != edges.end(); ++i)
  {
    if(!(*i)->fixed_propagate(*this))
    {
      for( ; i != edges.begin(); --i) (*i)->undo(*this);

      (*edges.begin())->undo(*this);

      restore();

      return false;
    }
  }

  ch = edge::connection++;

  return true;
}

void node::unfix()
{
  if(!fixed()) return;

  --edge::connection;

  undo( **edges.begin() );

  (*edges.begin())->undo(*this);
}

#ifndef __TESTING_NODES
double node::probability(const SAGE::MLOCUS::inheritance_model& im, bool use)
{
  if(!state()) return 1;
 
  if(!fixed()) return -1;

  if(!use) return 1.0;
  
  if((!my_is_sex_null && f_allele.is_sex_allele()) ||
     (my_is_sex_null && !f_allele.is_sex_allele()))
  {
//    cout << "REQUIRED to be different: " << my_id << " = " << f_allele.name() << endl;
  }

  return im.get_allele(f_allele.name()).frequency();
}

#endif

bool node::add_edge(edge& e, int edge_end)
{
#if 0
  cout << " node::add_edge()..." << endl;
#endif

  // If this is the first edge
  if(edges.empty())
  {
    // We add the edge and set the possible alleles

    edges.push_front(&e);

    a[0] = e.get_allele(0);
    a[1] = e.get_allele(1);

    // Determine procedure based on edge type.  An unsorted edge simply
    // propegates itself.  Sorted and homozygous edges must use allele1,
    // while phased edges specify the edge to be used in edge_end.

    switch(e.state())
    { 
      case edge::unsorted   :
        e.propagate(*this);
        break;

      case edge::homozygous :
      case edge::sorted     :
        edge_end = 1;

      case edge::phased     :
      
        allele choice = e.get_allele(edge_end);
#ifndef __TESTING_NODES
        if((choice.is_sex_allele() && !my_is_sex_null) ||
           (!choice.is_sex_allele() && my_is_sex_null))
        {
          edges.pop_front();
          return false;
        }
#endif
        f_allele = choice;

        ch = edge::connection;
    }
    
    ++edge::connection;

    return true;
  }

  // If we're already fixed, an unsorted edge propegates the information,
  // while the other types of edges check for valid state.  They don't need
  // to propagate any information.

  if(fixed())
  {
    switch(e.state())
    { 
      case edge::unsorted   :
        if(!e.fixed_propagate(*this)) return false;
        break;

      case edge::homozygous :
      case edge::sorted     :
        edge_end = 1;

      case edge::phased     :
        if(e.get_allele(edge_end) != f_allele) return false;
    }

    edges.push_front(&e);

    ++edge::connection;

    return true;
  }

  // At this point, we know we're not fixed, but we have potential alleles.

  // If e is not unsorted, then we can determine which allele it *must* be.

  if(e.state() != edge::unsorted)
  {
    // Phased alleles must be edge_end, while the other types muse be allele 1

    if(e.state() != edge::phased)
      edge_end = 1;

    if     (e.get_allele(edge_end) == a[1]) save(1);
    else if(e.get_allele(edge_end) == a[0]) save(0);
    else                                    return false;

    if(!propagate(e))
    {
      restore();
      return false;
    }

    ch = edge::connection++;

    edges.push_front(&e);

    return true;
  }
  
  // At this point, we and the edge are both unsorted

  // Determine if the alleles for the edge match the alleles for the node. 
  // If they do, we add the edge and do nothing more.

  if(e.get_allele(0) != a[0] && e.get_allele(1) != a[1])
    std::swap(a[0], a[1]);
  
  if(e.get_allele(0) == a[0])
  {
    if(e.get_allele(1) == a[1])
    {
      edges.push_front(&e);
      e.propagate(*this);

      ++edge::connection;

      return true;
    }

    // If only the 0 allele matches, we save it

    save(0);
  }
  else if(e.get_allele(1) == a[1]) save(1);
  else                             return false;

  if(!propagate(e))
  {
    restore();
    return false;
  }

  if(!e.fixed_propagate(*this))
  {
    undo(e);
    return false;
  }

  edges.push_front(&e);

  ch = edge::connection++;

  return true;
}

void node::remove_edge()
{
  if(edges.empty()) return;

  --edge::connection;

  if(ch == edge::connection)
  {
    undo( **edges.begin() );
  }

  (*edges.begin())->undo(*this, true);
    
  edges.pop_front();
}
 
bool node::fix(edge& dont, allele al)
{
  if(fixed())
    if(ch != -1 || al != f_allele) return false;
    else                           return true;

  if     (al == a[0]) save(0);
  else if(al == a[1]) save(1);
  else                return false;

  if(!propagate(dont))
  {
    restore();
    return false;
  }

  ch = edge::connection;

  return true;
}

bool node::propagate(edge& dont)
{
  for(iterator i = edges.begin(); i != edges.end(); ++i)
  {
    if((*i)->id() == dont.id()) continue;

    if(!(*i)->fixed_propagate(*this))
    {
      for( ; i != edges.begin(); --i) (*i)->undo(*this);

      (*edges.begin())->undo(*this);

      return false;
    }
  }

  return true;
}

void node::undo(edge& dont)
{
  if(ch == -1) return;

  ch = -1;

  for(iterator i = edges.begin(); i != edges.end(); ++i)
    if((*i)->id() != dont.id()) (*i)->undo(*this);

  restore();
}

bool edge::propagate(node& n)
{
  if(my_state != unsorted)
  {
    std::cout << "Shouldn't ever see this -- edge propagation" << std::endl;
    return true;
  }

  if(!nodes[0])
  {
    nodes[0] = &n;

    return true;
  }

  if(nodes[1])
  {
    std::cerr << "Shouldn't ever see this -- edge propagation2" << std::endl;
    return true;
  }

  nodes[1] = &n;
  
  return true;
}

bool edge::fixed_propagate(node& n)
{
  if(my_state != unsorted)
    return true;

  if(!nodes[0])
  {
    if(my_genotype.allele1() != n.f_allele)
      if(my_genotype.allele2() != n.f_allele)
        return false;
      else
        my_genotype.flip();

    ch = connection;
    my_state = sorted;

    nodes[0] = &n;

    return true;
  }

  if(!nodes[1])
  {
    if(my_genotype.allele2() != n.f_allele)
      if(my_genotype.allele1() != n.f_allele)
        return false;
      else
        my_genotype.flip();

    nodes[1] = &n;

    ch = connection;
    my_state = sorted;

    if(!nodes[0]->fix(*this, my_genotype.allele1()))
    {
      nodes[1] = NULL;
      ch = -1;
      my_state = unsorted;
      return false;
    }

    return true;
  }

  if(!n.fixed())
  {
    std::cout << "Shouldn't ever see this -- edge propagation4" << std::endl;
    return true;
  }

  if(&n == nodes[0])
  {
    if(my_genotype.allele1() != n.f_allele)
      my_genotype.flip();

    if(!nodes[1]->fix(*this, my_genotype.allele2())) return false;

    my_state = sorted;
    ch = connection;

    return true;
  }

  if(my_genotype.allele2() != n.f_allele)
    my_genotype.flip();

  if(!nodes[0]->fix(*this, my_genotype.allele1())) return false;

  my_state = sorted;
  ch = connection;

  return true;
}

void edge::undo(node& n, bool b)
{
  if(my_state == homozygous || my_state == phased || (ch != connection && !b)) return;

  if(my_state == unsorted && !b)
  {
    std::cerr << "Shouldn't see this -- edge undo" << std::endl;
    return;
  }

  if(&n == nodes[1])
  {
    nodes[0]->undo(*this);
    if(b) nodes[1] = NULL;
  }
  else
  {
    if(b) nodes[0] = NULL;
    else if(nodes[1]) nodes[1]->undo(*this);
  }

  ch = -1;
  my_state = unsorted;
}

}

