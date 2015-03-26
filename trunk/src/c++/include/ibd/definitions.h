#ifndef IBD_DEF_H
#define IBD_DEF_H

//****************************************************************************
//* File:      definitions.h                                                 *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation              - yjs Jun. 11 *
//*                                                                          *
//* Notes:     This header file defines IBD utility structure.               *
//*                                                                          *
//* Copyright (c) 2011 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "globals/SAGEConstants.h"
#include "numerics/corinfo.h"
#include "rped/genome_description.h"
#include "lvec/lvector.h"
#include "lvec/meiosis_map.h"
#include "lvec/mpoint_like.h"
#include "pairs/relpair.h"

namespace SAGE {

const double INF    = std::numeric_limits<double>::infinity();
const double NE_INF = -std::numeric_limits<double>::infinity();
const double EPS    = std::numeric_limits<double>::epsilon();

typedef Likelihood_Vector                        lvector;
typedef RPED::genome_description::region_type    region_type;

typedef FPED::SubpedigreeConstPointer            sped_pointer;
typedef FPED::MemberConstPointer                 mem_pointer;
typedef pair<mem_pointer, mem_pointer>           id_pair;

typedef pair_generator::pair_type                pair_type;
typedef MLOCUS::GenotypeModelType                gmodel_type;

enum pair_sex_type { MM        = 0, MF, FF, UNKNOWN };
enum pair_x_type   { M_M       = 0, M_F, F_F,
                     M_F_M, M_F_F, F_F_F,
                     MFM, MFF, FFM, FFF,
                     M_FM, M_FF, F_FM, F_FF, M_MF, F_MF,
                     FM_MF, MF_FM, MF_FF, FF_FF, FM_FF, FM_FM, INVALID };

struct ibd_option_type
{
  ibd_option_type()
  {
    max_pedigree      = "18";
    scan_type         = "markers";
    allow_loops       = "off";
    ibd_mode          = "multipoint";
    split_pedigrees   = "no";
    use_simulation    = "no";

    exact             = true;
    x_linked          = false;
    old_ibd_format    = false;
    ibd_state_out     = false;
  };

  string title;
  string region;
  string max_pedigree;
  string scan_type;   
  string allow_loops; 
  string ibd_mode;    
  string split_pedigrees;
  string use_simulation; 

  bool   exact;
  bool   x_linked;
  bool   old_ibd_format;
  bool   ibd_state_out;
};

struct ibd_marker_info
{
  ibd_marker_info() : name(""), distance(0.0), type(MLOCUS::AUTOSOMAL) {}
  ibd_marker_info(const string& n, double d, gmodel_type m)
    : name(n), distance(d), type(m) {}

  string        name;
  double        distance;
  gmodel_type   type;
};

struct member_pair : public id_pair
{
  member_pair()
  {
    first = second = NULL;
  }
  
  member_pair(mem_pointer a, mem_pointer b)
  {
    first  = a;
    second = b;
  }
                          
  template<class U, class V>
  member_pair(const pair<U, V> &p)
  {
    first  = p.first;
    second = p.second;
  }
  
  bool operator==(const member_pair& p) const
  {
    return (first == p.first && second == p.second);
  }

  bool operator!=(const member_pair& p)
  {
    return !(*this == p);
  }
  
  bool operator<(const member_pair& p) const
  {
//    if(first < p.first)   return true;
//    if(first > p.first)   return false;
//    if(second < p.second) return true;
//    return false;
    if( first->pedigree()->name() == p.first->pedigree()->name() )
      if( first->name() == p.first->name() )
        return second->name() < p.second->name();
      else
        return first->name() < p.first->name();
    else
      return second->pedigree()->name() < p.second->pedigree()->name();
  }
};

struct member_pair_less : public std::binary_function<member_pair, member_pair, bool>
{
  member_pair_less() {}

  bool operator()(const member_pair& p1, const member_pair& p2) const
  {
    if( p1.first->pedigree()->name() == p2.first->pedigree()->name() )
      if( p1.first->name() == p2.first->name() )
        return p1.second->name() < p2.second->name();
      else
        return p1.first->name() < p2.first->name();
    else
      return p1.second->pedigree()->name() < p2.second->pedigree()->name();
  }
};

typedef map<member_pair, size_t>  rel_map;
//typedef map<member_pair, size_t, member_pair_less >  rel_map;

struct ibd_pair_info
{
  ibd_pair_info() : index((size_t)-1) {}

  ibd_pair_info(mem_pointer i1, mem_pointer i2, pair_type p, size_t n = (size_t)-1)
    : pair(i1,i2), type(p), index(n) {}

  ibd_pair_info(mem_pointer i1, mem_pointer i2, pair_type p, pair_sex_type s, size_t n = (size_t)-1)
    : pair(i1,i2), type(p), index(n), pair_sex(s), x_type(INVALID) {}

  member_pair            pair;
  pair_type              type;

  size_t                 index;

  pair_sex_type          pair_sex;
  pair_x_type            x_type;
};

struct ibd_pair_info_less : public std::binary_function<ibd_pair_info, ibd_pair_info, bool>
{
  ibd_pair_info_less() {}

  bool operator()(const ibd_pair_info& p1, const ibd_pair_info& p2) const
  {
    if( p1.pair.first->pedigree()->name() == p2.pair.first->pedigree()->name() )
      if( p1.pair.first->name() == p2.pair.first->name() )
        return p1.pair.second->name() < p2.pair.second->name();
      else
        return p1.pair.first->name() < p2.pair.first->name();
    else
      return p1.pair.second->pedigree()->name() < p2.pair.second->pedigree()->name();
  }
};

struct ibd_probability_info
{
  ibd_probability_info(size_t n = 0)
  {
    f0s.resize(n);
    f1mp.resize(n);
    f2s.resize(n);
  }

  vector<double>         f0s;
  vector<double>         f1mp;
  vector<double>         f2s;
};

typedef pair< double, vector<size_t> >   ibd_lvec_probability;
typedef vector< ibd_lvec_probability >   a_marker_ibd_state;
typedef vector< a_marker_ibd_state >     markers_ibd_state;

struct ibd_state_info
{
  ibd_state_info() {}

  size_t find_pair_ids(id_pair pid) const
  {
    for( size_t i = 0; i < pair_ids.size(); ++i )
    {
      if(    (pid.first == pair_ids[i].first && pid.second == pair_ids[i].second)
          || (pid.first == pair_ids[i].second && pid.second == pair_ids[i].first) )
        return i;
    }
    return (size_t) -1;
  }

  size_t find_pair_ids(mem_pointer p1, mem_pointer p2) const
  {
    return find_pair_ids(make_pair(p1, p2));
  }

  vector<id_pair>    pair_ids;
  markers_ibd_state  ibd_states;
};


//
// inline functions
//

inline bool is_self(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return true;
  return false;
}

inline bool is_self(const id_pair i)
{
  return is_self(i.first, i.second);
}

inline bool is_sib(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  if( !i1->parent1() || !i1->parent2() || !i2->parent1() || !i2->parent2() )
    return false;

  return (i1->parent1() == i2->parent1() || i1->parent1() == i2->parent2()
       || i1->parent2() == i2->parent1() || i1->parent2() == i2->parent2() );
}

inline bool is_sib(const id_pair i)
{
  return is_sib(i.first, i.second);
}

inline bool is_fsib(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  if( !is_sib(i1,i2) ) return false;

  return ((i1->parent1() == i2->parent1() && i1->parent2() == i2->parent2())
       || (i1->parent1() == i2->parent2() && i1->parent2() == i2->parent1()));
}

inline bool is_fsib(const id_pair i)
{
  return is_fsib(i.first, i.second);
}

inline bool is_hsib(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  return is_sib(i1,i2) && !is_fsib(i1,i2);
}

inline bool is_hsib(const id_pair i)
{
  return is_hsib(i.first, i.second);
}

inline bool is_maternal_hsib(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( !is_hsib(i1,i2) )
    return false;
  
  if(    (i1->parent1() && i1->parent1() == i2->parent1() && i1->parent1()->is_female())
      || (i1->parent1() && i1->parent1() == i2->parent2() && i1->parent1()->is_female())
      || (i1->parent2() && i1->parent2() == i2->parent1() && i1->parent2()->is_female())
      || (i1->parent2() && i1->parent2() == i2->parent2() && i1->parent2()->is_female()) )
    return true;

  return false;
}

inline bool is_maternal_hsib(const id_pair i)
{
  return is_maternal_hsib(i.first, i.second);
}

inline bool is_paternal_hsib(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( !is_hsib(i1,i2) )
    return false;

  return !is_maternal_hsib(i1,i2);
}

inline bool is_paternal_hsib(const id_pair i)
{
  return is_paternal_hsib(i.first, i.second);
}

inline bool is_first_ind_grandp(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  // case 1: i1 is granp of i2
  //
  if( i2->parent1() && i2->parent2() )
  {
    FPED::MateConstIterator m1 = i1->mate_begin();
    for( ; m1 != i1->mate_end(); ++m1 )
    {
      FPED::OffspringConstIterator o1 = i1->offspring_begin(m1);
      for( ; o1 != i1->offspring_end(); ++o1 )
        if( o1->index() == i2->parent1()->index() || o1->index() == i2->parent2()->index() )
          return true;
    }
  }

  return false;
}

inline bool is_grandp(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  return is_first_ind_grandp(i1, i2) || is_first_ind_grandp(i2, i1);
}

inline bool is_grandp(const id_pair i)
{
  return is_grandp(i.first, i.second);
}

inline bool is_m_grandp(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( is_first_ind_grandp(i1, i2) )
    return (  ( i2->parent1() && i2->parent1()->is_female()
                && (    (i2->parent1()->parent1() && i2->parent1()->parent1() == i1)
                     || (i2->parent1()->parent2() && i2->parent1()->parent2() == i1) ) )
           || ( i2->parent2() && i2->parent2()->is_female()
                && (    (i2->parent2()->parent1() && i2->parent2()->parent1() == i1)
                     || (i2->parent2()->parent2() && i2->parent2()->parent2() == i1) ) ) );

  if( is_first_ind_grandp(i2, i1) )
    return (  ( i1->parent1() && i1->parent1()->is_female()
                && (    (i1->parent1()->parent1() && i1->parent1()->parent1() == i2)
                     || (i1->parent1()->parent2() && i1->parent1()->parent2() == i2) ) )
           || ( i1->parent2() && i1->parent2()->is_female()
                && (    (i1->parent2()->parent1() && i1->parent2()->parent1() == i2)
                     || (i1->parent2()->parent2() && i1->parent2()->parent2() == i2) ) ) );


  return false;
}

inline bool is_m_grandp(const id_pair i)
{
  return is_m_grandp(i.first, i.second);
}

inline bool is_p_grandp(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  return is_grandp(i1,i2) && !is_m_grandp(i1,i2);
}

inline bool is_p_grandp(const id_pair i)
{
  return is_p_grandp(i.first, i.second);
}

inline bool is_first_ind_avunc(const mem_pointer i1, const mem_pointer i2)
{
  if( i2->parent1() )
  {
    FPED::SiblingConstIterator s1 = i2->parent1()->sibling_begin();
    for( ; s1 != i2->parent1()->sibling_end(); ++s1 )
      if( s1->index() == i1->index() )
        return true;
  }

  if( i2->parent2() )
  {
    FPED::SiblingConstIterator s2 = i2->parent2()->sibling_begin();
    for( ; s2 != i2->parent2()->sibling_end(); ++s2 )
      if( s2->index() == i1->index() )
        return true;
  }

  return false;
}

inline bool is_avunc(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  return is_first_ind_avunc(i1, i2) || is_first_ind_avunc(i2, i1);
}

inline bool is_avunc(const id_pair i)
{
  return is_avunc(i.first, i.second);
}

inline bool is_m_avunc(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( is_first_ind_avunc(i1, i2) )
  {
    FPED::SiblingConstIterator s1 = i2->parent1()->sibling_begin();
    for( ; s1 != i2->parent1()->sibling_end(); ++s1 )
      if( s1->index() == i1->index() && i2->parent1()->is_female() )
        return true;

    FPED::SiblingConstIterator s2 = i2->parent2()->sibling_begin();
    for( ; s2 != i2->parent2()->sibling_end(); ++s2 )
      if( s2->index() == i1->index() && i2->parent2()->is_female() )
        return true;
  }

  if( is_first_ind_avunc(i2, i1) )
  {
    FPED::SiblingConstIterator s1 = i1->parent1()->sibling_begin();
    for( ; s1 != i1->parent1()->sibling_end(); ++s1 )
      if( s1->index() == i2->index() && i1->parent1()->is_female() )
        return true;

    FPED::SiblingConstIterator s2 = i1->parent2()->sibling_begin();
    for( ; s2 != i1->parent2()->sibling_end(); ++s2 )
      if( s2->index() == i2->index() && i1->parent2()->is_female() )
        return true;
  }

  return false;
}

inline bool is_m_avunc(const id_pair i)
{
  return is_m_avunc(i.first, i.second);
}

inline bool is_cousin(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( !i1->parent1() || !i1->parent2() || !i2->parent1() || !i2->parent2() )
    return false;
  
  mem_pointer i1_p1 = i1->parent1();
  mem_pointer i1_p2 = i1->parent2();

  mem_pointer i2_p1 = i2->parent1();
  mem_pointer i2_p2 = i2->parent2();

  FPED::SiblingConstIterator s1 = i1_p1->sibling_begin();
  for( ; s1 != i1_p1->sibling_end(); ++s1 )
    if( s1->index() == i2_p1->index() || s1->index() == i2_p2->index() )
      return true;

  s1 = i1_p2->sibling_begin();
  for( ; s1 != i1_p2->sibling_end(); ++s1 )
    if( s1->index() == i2_p1->index() || s1->index() == i2_p2->index() )
      return true;

  return false;
}

inline bool is_cousin(const id_pair i)
{
  return is_cousin(i.first, i.second);
}

inline bool is_p_p_cousin(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;

  if( !is_cousin(i1, i2) )
    return false;
  
  mem_pointer i1_p = NULL;
  mem_pointer i2_p = NULL;
  
  if( i1->parent1()->is_male() )
    i1_p = i1->parent2();
  else if( i1->parent2()->is_male() )
    i1_p = i1->parent1();

  if( i2->parent1()->is_male() )
    i2_p = i2->parent2();
  else if( i2->parent2()->is_male() )
    i2_p = i2->parent1();

  if( i1_p == NULL || i2_p == NULL )
    return false;

  FPED::SiblingConstIterator s1 = i1_p->sibling_begin();
  for( ; s1 != i1_p->sibling_end(); ++s1 )
    if( s1->index() == i2_p->index() )
      return true;

  return false;
}

inline bool is_p_p_cousin(const id_pair i)
{
  return is_p_p_cousin(i.first, i.second);
}

inline bool is_m_m_cousin(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  
  if( !is_cousin(i1, i2) )
    return false;

  mem_pointer i1_p = NULL;
  mem_pointer i2_p = NULL;
  
  if( i1->parent1()->is_female() )
    i1_p = i1->parent2();
  else if( i1->parent2()->is_female() )
    i1_p = i1->parent1();

  if( i2->parent1()->is_female() )
    i2_p = i2->parent2();
  else if( i2->parent2()->is_female() )
    i2_p = i2->parent1();

  if( i1_p == NULL || i2_p == NULL )
    return false;

  FPED::SiblingConstIterator s1 = i1_p->sibling_begin();
  for( ; s1 != i1_p->sibling_end(); ++s1 )
    if( s1->index() == i2_p->index() )
      return true;

  return false;
}

inline bool is_m_m_cousin(const id_pair i)
{
  return is_m_m_cousin(i.first, i.second);
}

inline bool is_p_m_cousin(const mem_pointer i1, const mem_pointer i2)
{
  if(i1 == i2) return false;
  
  if( !is_cousin(i1, i2) )
    return false;

  mem_pointer i1_p = NULL;
  mem_pointer i2_p = NULL;
  
  if( i1->parent1()->is_male() )
    i1_p = i1->parent2();
  else if( i1->parent2()->is_male() )
    i1_p = i1->parent1();

  if( i2->parent1()->is_female() )
    i2_p = i2->parent2();
  else if( i2->parent2()->is_female() )
    i2_p = i2->parent1();

  if( i1_p == NULL || i2_p == NULL )
    return false;

  FPED::SiblingConstIterator s1 = i1_p->sibling_begin();
  for( ; s1 != i1_p->sibling_end(); ++s1 )
    if( s1->index() == i2_p->index() )
      return true;

  return false;
}

inline bool is_p_m_cousin(const id_pair i)
{
  return is_p_m_cousin(i.first, i.second);
}

inline bool is_cluster_parent(const mem_pointer p, const set<mem_pointer>& parent_set)
{
  set<mem_pointer>::const_iterator pi = parent_set.begin();
  for( ; pi != parent_set.end(); ++pi )
    if( *pi == p )
      return true;

  return false;
}

inline bool is_unrelated_pair(const id_pair &p1, const id_pair &p2)
{
  set<mem_pointer> parent_set1;
  
  parent_set1.insert(p1.first->parent1());
  parent_set1.insert(p1.first->parent2());
  parent_set1.insert(p1.second->parent1());
  parent_set1.insert(p1.second->parent2());
  
  set<mem_pointer> parent_set2;
  
  parent_set2.insert(p2.first->parent1());
  parent_set2.insert(p2.first->parent2());
  parent_set2.insert(p2.second->parent1());
  parent_set2.insert(p2.second->parent2());
  
  set<mem_pointer>::const_iterator pi = parent_set2.begin();
  for( ; pi != parent_set2.end(); ++pi )
  {
    if( is_cluster_parent(*pi, parent_set1) )
      return false;
  }

  return true;
}

inline int sibs_shared(const id_pair &p1, const id_pair &p2)
{
  if( p1.first == p2.first  && p1.second == p2.second ) return 2;
  if( p1.first == p2.first  || p1.second == p2.second ||
      p1.first == p2.second || p1.second == p2.first ) return 1;

  if( is_unrelated_pair(p1, p2) )  return 3;

  return 0;
}

inline double zlog(double d)
{
  if(d < std::numeric_limits<double>::epsilon())
    return 0.0;
  return log(d);
}

inline bool is_mm_pair(const mem_pointer i1, const mem_pointer i2)
{
  return ( i1->is_male() && i2->is_male() );
}

inline bool is_mm_pair(const id_pair i)
{
  return is_mm_pair(i.first, i.second);
}

inline bool is_ff_pair(const mem_pointer i1, const mem_pointer i2)
{
  return ( i1->is_female() && i2->is_female() );
}

inline bool is_ff_pair(const id_pair i)
{
  return is_ff_pair(i.first, i.second);
}

inline bool is_mf_pair(const mem_pointer i1, const mem_pointer i2)
{
  return (   (i1->is_male() && i2->is_female())
          || (i1->is_female() && i2->is_male()) );
}

inline bool is_mf_pair(const id_pair i)
{
  return is_mf_pair(i.first, i.second);
}

inline bool is_brother_brother(const mem_pointer i1, const mem_pointer i2)
{
  if( !is_fsib(i1, i2) && !is_hsib(i1, i2) )
    return false;

  return is_mm_pair(i1, i2);
}

} // end of namespace SAGE

#endif
