//============================================================================
// File:     stats.cpp
//                                                                          
// Author:   
//                                                                          
// History:  8/00 - Modified and extended to utilize pair_generator class
//                  and to include trait information.  - Dan Baechle
//                                                                          
// Notes:    Implements the following classes -
//              General_Stats
//              Ped_stats
//              MP_stats
//              Base_trait_stats
//              Binary_trait_stats
//              Cont_trait_stats
//              Cmpd_trait_stats
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "pedinfo/stats.h"
#include "error/internal_error.h"
using namespace std;

namespace SAGE {
namespace RPED {

//============================================================================
// IMPLEMENTATION:  General_Stats
//============================================================================

void 
General_Stats::init()
{
  my_total_founders    = my_total_nonfounders  = my_num_gen             = 0;  
  my_total_unknown_sex = my_unconnecteds       = my_total_subpedigrees  = 0;
  my_total_inds        = my_total_females      = my_total_males         = 0;
  my_marriage_loops    = my_non_marriage_loops = my_max_bits            = 0;
  
  for(size_t i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    my_pairs[i] = 0;
  }
  
  my_valid = true;
}

General_Stats& 
General_Stats::operator+=(const General_Stats& gs)
{
  my_total_founders     += gs.founder_count();
  my_total_nonfounders  += gs.nonfounder_count();
  my_total_females      += gs.female_count();
  my_total_males        += gs.male_count();
  my_total_inds         += gs.member_count();
  my_total_unknown_sex  += gs.unknown_sex_count();
  my_marriage_loops     += gs.marriage_loops();
  my_non_marriage_loops += gs.non_marriage_loops();
  my_unconnecteds       += gs.unconnecteds();
  my_total_subpedigrees += gs.total_subpedigrees();

  my_max_bits            = max(my_max_bits, gs.likelihood_bits());

  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    if(!finite(gs.pairs(static_cast<pg::pair_type>(i + 1))))  
      continue;
    my_pairs[i] += gs.pairs(static_cast<pg::pair_type>(i + 1));
  }

  my_sibship += gs.sibship();
  my_sibship_size_freq += gs.sib_size_freq();
  my_multiple_mates.insert(gs.my_multiple_mates.begin(), gs.my_multiple_mates.end());
  my_cons_pairs.insert(gs.my_cons_pairs.begin(), gs.my_cons_pairs.end());
  
  // Add trait data.
  if(my_binary_trait_stats.size() == gs.my_binary_trait_stats.size())
  {
    for(size_t i = 0; i < my_binary_trait_stats.size(); ++i)
    {
      if(my_binary_trait_stats[i].trait() == gs.my_binary_trait_stats[i].trait())
      {
        my_binary_trait_stats[i] += gs.my_binary_trait_stats[i];
      }
    }
  }
  
  if(my_cont_trait_stats.size() == gs.my_cont_trait_stats.size())
  {
    for(size_t i = 0; i < my_cont_trait_stats.size(); ++i)
    {
      if(my_cont_trait_stats[i].trait() == gs.my_cont_trait_stats[i].trait())
      {
        my_cont_trait_stats[i] += gs.my_cont_trait_stats[i];
      }
    }
  }
  
  if(my_cmpd_trait_stats.size() == gs.my_cmpd_trait_stats.size())
  {
    for(size_t i = 0; i < my_cmpd_trait_stats.size(); ++i)
    {
      if(my_cmpd_trait_stats[i].traits() == gs.my_cmpd_trait_stats[i].traits())
      {
        my_cmpd_trait_stats[i] += gs.my_cmpd_trait_stats[i];
      }
    }
  }

  return *this;
}

// Print trait information to std out for testing.
//
void
General_Stats::print_traits()
{
  for(size_t i = 0; i < my_binary_trait_stats.size(); ++i)
  {
    my_binary_trait_stats[i].print();
  }
  
  for(size_t i = 0; i < my_cont_trait_stats.size(); ++i)
  {
    my_cont_trait_stats[i].print();
  }
}

//============================================================================
// IMPLEMENTATION:  Ped_stats
//============================================================================
//
// Without trait.
//

Ped_stats::Ped_stats(const RefPedigree* p, cerrormultistream &e) 
      : General_Stats(e), my_pedigree(NULL)
{
  invalidate();
  compute(p);
}

// With trait.
//
Ped_stats::Ped_stats(const RefPedigree* p, size_t trait, cerrormultistream &e) 
      : General_Stats(e), my_pedigree(NULL)
{
  invalidate();
  compute(p, trait);
}

// With compound trait.
//
Ped_stats::Ped_stats(const RefPedigree* p, std::vector<size_t> traits, cerrormultistream &e) 
      : General_Stats(e), my_pedigree(NULL)
{
  invalidate();
  compute(p, traits);
}

Ped_stats::Ped_stats(cerrormultistream &e) 
      : General_Stats(e), my_pedigree(NULL)
{
  invalidate();
}

Ped_stats::~Ped_stats()  
{}

// Without trait.
//
bool 
Ped_stats::compute(const RefPedigree* pedig)
{
  if(!pedig) 
  {
    return false;
  }

  if( pedig==ped() )
    return true;
    
  non_trait_computations(pedig);

  my_valid = true;
  return true;
}

// With simple trait.
//
bool 
Ped_stats::compute(const RefPedigree* pedig, size_t trait)
{
  if(!pedig) 
  {
    return false;
  }

  if( pedig==ped() )
    return true;
    
  non_trait_computations(pedig);
  compute_traits(trait);

  my_valid = true;
  return true;
}

// With compound trait.
//
bool 
Ped_stats::compute(const RefPedigree* pedig, std::vector<size_t> traits)
{
  if(!pedig) 
  {
    return false;
  }

  if( pedig==ped() )
    return true;
    
  non_trait_computations(pedig);
  compute_traits(traits);

  my_valid = true;
  return true;
}

// Assumes pedig is non-zero and not equal to my_pedigree.
// This is checked in compute().
//
void
Ped_stats::non_trait_computations(const RefPedigree* pedig)
{
  init(); 

  //size_t male_sib, female_sib;
  size_t cur_sib;
  list<RefPedigree::member_const_pointer>  founder_list; 

  RefPedigree::family_const_iterator     family ;
  RefPedigree::member_const_iterator     ind;

  my_pedigree = pedig; 

  my_total_inds         = pedig->member_count();
  my_unconnecteds       = pedig->unconnected_count();
  my_total_subpedigrees = pedig->subpedigree_count();

// This is _not_ the current definition
//my_total_founders    += pedig->unconnected_count();

  // obtaining sex counts info, founder/nonfounder info, ind's w. multiple mates
  for(ind = pedig->member_begin(); ind != pedig->member_end(); ++ind)
  {
    if     ( ind->is_female() )  ++my_total_females;
    else if( ind->is_male()   )  ++my_total_males;
    else                         ++my_total_unknown_sex;
    
    if((!ind->parent1() || !ind->parent2()) && ind->offspring_count()) 
     { ++my_total_founders; founder_list.push_back(&(*ind)); }
    
    if( ind->parent1() && ind->parent2() )
      ++my_total_nonfounders;
      
    // - Individuals w. multiple mates.  6-18-3, djb
    //
    if(ind->mate_count() > 1)
    {
      full_name  ind_name(ind->pedigree()->name(), ind->name());
      vector<string>  mates;
      
      RefMultiPedigree::mate_const_iterator  m_iter;
      for(m_iter = ind->mate_begin(); m_iter != ind->mate_end(); ++m_iter)
      {
        mates.push_back(m_iter->mate().name());
      }
      
      assert(my_multiple_mates.insert(mate_map::value_type(ind_name, mates)).second);
    }
  }

  pp_vector  pairs;
  for(family = pedig->family_begin(); family != pedig->family_end(); ++family) 
  {
    cur_sib = family->offspring_count();
    my_sibship += cur_sib;
    
    // - Consanguineous mating pairs.  6-19-3, djb
    //
    if(consanguineous(family))
    {
      pairs.push_back(parental_pair(family->parent1()->name(), family->parent2()->name()));
    }
  }
  
  my_cons_pairs.insert(pp_map::value_type(pedig->name(), pairs));

  compute_pairs(pedig);
  compute_loops(pedig);

  if(non_marriage_loops() == 0 )
    get_generation(founder_list);

  compute_bits(pedig);
}

// - Are the parents of the given family a consanguineous mating pair?
//   See .../c++/pedinfo/cons_matings.txt for algorithm description.
//
bool
Ped_stats::consanguineous(RefPedigree::family_const_iterator& family)
{
  set<string>  set1;
  set<string>  set2;
  set<string>  intersection;
  
  fill_set(set1, family->parent1());
  fill_set(set2, family->parent2());
  
  set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), 
                   inserter(intersection, intersection.begin())       );
  
  return  ! intersection.empty();
}

// - Fill set w. names of ind, his parent, his parent's parents ...
//   This is a binary tree traversal w. ind as the root node.
//
void 
Ped_stats::fill_set(std::set<string>& s, RefPedigree::member_const_pointer ind)
{
  if(ind)
  {
    s.insert(ind->name());
    fill_set(s, ind->parent1());
    fill_set(s, ind->parent2());
  }
}

void 
Ped_stats::get_generation(list<RefPedigree::member_const_pointer> &fl)
{
  // - Pedigrees consisting soley of unconnected individuals will
  //   be said to have 1 generation instead of 0. -djb 3/25/2
  //
  my_num_gen = 1;

  list<RefPedigree::member_const_pointer>::iterator i;
  deque<RefPedigree::member_const_pointer> qu;

  size_t generation;
  RefPedigree::member_const_pointer current_person;
  RefPedigree::offspring_const_iterator child;

  for(i=fl.begin(); i != fl.end(); ++i)
  {
    generation = 1;

    qu.push_back( *i );
    qu.push_back(0);
  
    while(true)  
    {
      current_person = qu.front();
      qu.pop_front();
      if( qu.empty() ) break;

      if(current_person == 0)
      {
         qu.push_back(0);
         current_person = qu.front();
         qu.pop_front();
         ++generation;
      }
      
      if( current_person == 0 || qu.empty() ) break;

      if( current_person->mate_count() > 0 )
      {  
        RefPedigree::mate_const_iterator mate = current_person->mate_begin();
        for( ; mate != current_person->mate_end(); ++mate )
        {
          child = current_person->offspring_begin(mate);
          for( ; child != current_person->offspring_end(); ++child)
            qu.push_back( &(*child) );
        }
      }
    }

    if(generation > my_num_gen ) my_num_gen = generation; 
  } 
}

// - Generate data for all traits.
//
void
Ped_stats::compute_traits()
{
  if(my_pedigree != 0)
  {
    RefMPedInfo            multi = my_pedigree->multipedigree()->info();
    RefTraitInfo           rti;
    RefTraitInfo::trait_t  type;
    for(size_t i = 0; i < multi.trait_count(); ++i)
    {
      rti   = multi.trait_info(i);
      type  = rti.type();
      if(type == RefTraitInfo::binary_trait)
      {
        Binary_trait_stats bts(my_pedigree, i, errors);    
        my_binary_trait_stats.push_back(bts);
      }
      else if(type == RefTraitInfo::continuous_trait ||
              type == RefTraitInfo::discrete_trait)
      {
        Cont_trait_stats cts(my_pedigree, i, errors);
        my_cont_trait_stats.push_back(cts);
      }
    }
  }
}

// - Generate data for specified trait.
//
void
Ped_stats::compute_traits(size_t trait)
{
  if(my_pedigree != 0)
  {
    RefMPedInfo            multi = my_pedigree->multipedigree()->info();
    RefTraitInfo           rti;
    RefTraitInfo::trait_t  type;
    size_t                 t_count = multi.trait_count();
    if(trait < t_count)
    {
      rti   = multi.trait_info(trait);
      type  = rti.type();
      if(type == RefTraitInfo::binary_trait)
      {
        Binary_trait_stats bts(my_pedigree, trait, errors);    
        my_binary_trait_stats.push_back(bts);
      }
      else if(type == RefTraitInfo::continuous_trait ||
              type == RefTraitInfo::discrete_trait)
      {
        Cont_trait_stats cts(my_pedigree, trait, errors);
        my_cont_trait_stats.push_back(cts);
      }
    }
  }
}

// - Generate data for a compound trait.
//
void
Ped_stats::compute_traits(std::vector<size_t> traits)
{
  if(my_pedigree != 0)
  {
    Cmpd_trait_stats cmpd(my_pedigree, traits, errors);
    my_cmpd_trait_stats.push_back(cmpd);
  }
}

// - Generate pair counts.  Assumes calling function insures pedig is non-null.
//
void 
Ped_stats::compute_pairs(const RefPedigree* pedig)
{
  pg::pair_type t;
  pair_generator p_gen(const_cast<RefPedigree*>(pedig), pg::ALL_TYPES);    // ****** address const issue
  for(pair_generator::iterator iter = p_gen.begin(); iter != p_gen.end(); ++iter)
  {
    t = iter->type();
    if(t <= TOTAL_PAIR_TYPES)
    {
      ++my_pairs[t - 1];
    }
  }
}

// - Count loops in the pedigree.  Assumes calling function insures pedig is non-null.
//
void 
Ped_stats::compute_loops(const RefPedigree* pedig)
{
  LoopChecker  lck;

  if(lck.check_pedigree(pedig) && lck.loops())
  {
    my_marriage_loops     = lck.marriage_loops();
    my_non_marriage_loops = lck.non_marriage_loops();
  }
}

// - likelihood bits is maximum over _each_ subpedigree.  Assumes calling 
//   function insures pedig is non-null.
//
void
Ped_stats::compute_bits(const RefPedigree* pedig)
{
  if(pedig->subpedigree_count() > 1)
  {
    size_t sfounder;
    RefPedigree::subpedigree_const_iterator sp;
    RefSubpedigree::member_const_iterator mp;

    for(sp = pedig->subpedigree_begin(); sp != pedig->subpedigree_end(); ++sp)
    {
      sfounder = 0; 
      for(mp= sp->member_begin(); mp != sp->member_end(); ++mp)
        if( !mp->parent1() && !mp->parent2() )  
          ++sfounder;
    
      my_max_bits = max(my_max_bits, 2*(signed)sp->member_count()-3*sfounder);
    }
  }
  else
    my_max_bits = 2*my_total_nonfounders - my_total_founders;
}



//============================================================================
// IMPLEMENTATION:  MP_stats
//============================================================================
//
// Without trait.
//
MP_stats::MP_stats(const RefMultiPedigree* mp, cerrormultistream &e) 
      : General_Stats(e), my_multipedigree(NULL)
{
  invalidate();
  my_stats_type = MULTIPEDIGREE_STATS;
  compute(mp);
}

// With simple trait.
//
MP_stats::MP_stats(const RefMultiPedigree* mp, size_t trait, cerrormultistream &e) 
      : General_Stats(e), my_multipedigree(NULL)
{
  invalidate();
  my_stats_type = MULTIPEDIGREE_STATS;
  compute(mp, trait);
}

// With compound trait.
//
MP_stats::MP_stats(const RefMultiPedigree* mp, std::vector<size_t> traits, cerrormultistream &e) 
      : General_Stats(e), my_multipedigree(NULL)
{
  invalidate();
  my_stats_type = MULTIPEDIGREE_STATS;
  compute(mp, traits);
}

// Computation must be invoked after construction for this constructor.
//
MP_stats::MP_stats(cerrormultistream &e) 
      : General_Stats(e), my_multipedigree(NULL)
{ 
  invalidate(); 
  my_stats_type = MULTIPEDIGREE_STATS;
}

MP_stats::~MP_stats()  
{}

void
MP_stats::init()
{
  my_ped_stats.clear();
  General_Stats::init();
}

bool
MP_stats::trait_valid(size_t trait)
{
  if(my_multipedigree != 0)
  {
    RefMPedInfo multi = my_multipedigree->info();
    size_t t_count = multi.trait_count();
    if(trait < t_count)
    {
      RefTraitInfo            rti = multi.trait_info(trait);
      RefTraitInfo::trait_t   type = rti.type();
      if(type == RefTraitInfo::binary_trait)
      {
        // Create 'empty' Binary_trait_stats object in the multi-pedigree.
        my_binary_trait_stats.push_back(Binary_trait_stats(trait, errors, rti.name(), true));
        return true;
      }
      else if(type == RefTraitInfo::continuous_trait ||
              type == RefTraitInfo::discrete_trait)
      {
        // Create 'empty' Cont_trait_stats object in the multi-pedigree.
        my_cont_trait_stats.push_back(Cont_trait_stats(trait, errors, rti.name(), true));
        return true;
      }
    }
  }
  return false;
}

bool
MP_stats::cmpd_trait_valid(std::vector<size_t> traits)
{
  if(!traits.empty())
  {
    if(my_multipedigree != 0)
    {
      RefMPedInfo multi = my_multipedigree->info();
      size_t t_count = multi.trait_count();
      for(size_t t = 0; t < traits.size(); ++t)
      {
        // Need at least one valid trait for compound trait to be valid.     
        if(t < t_count)
        {
          // Create 'empty' Cmpd_trait_stats object in the multi-pedigree.
          std::string name = compute_traits_name(traits);
          my_cmpd_trait_stats.push_back(Cmpd_trait_stats(traits, errors, name, true));
          return true;
        }
      }
    }
  }
  return false;
}

// Without trait.
//
bool 
MP_stats::compute(const RefMultiPedigree* mp)
{
  if(mp == mped())      return true;

  my_multipedigree = mp;
  init();
  compute_pedigrees();
  if ( pedigree_count() <= 0 )
    return false;
    
  my_valid = true;
  return true;
}

// Without traits.
//
void
MP_stats::compute_pedigrees()
{
  if(my_multipedigree != 0)
  {
    RefMultiPedigree::pedigree_const_iterator ped;

    for( ped = my_multipedigree->pedigree_begin(); ped != my_multipedigree->pedigree_end(); ++ped)
    {
      my_ped_stats.push_back( Ped_stats(&(*ped), errors) );
      Ped_stats& ps = my_ped_stats.back();    

      if(!ps.valid())
      {
        my_ped_stats.pop_back();
        continue;
      }
     
      // Add the new pedigree General_Stats
      *this += ps;

      my_pedigree_size_info += ps.member_count();
      my_family_count_info  += ps.family_count();

      // Count 2n-f Frequencies
      my_likelihood_bits_freq.add(ps.likelihood_bits());
    
      // Count Family Size
      my_family_count_freq.add(ps.family_count());    

      // Count Generations
      my_generations_freq.add(ps.generations());
    }
  }
}

// With simple trait.
//
bool 
MP_stats::compute(const RefMultiPedigree* mp, size_t trait)
{
  if(mp == mped())      return true;

  my_multipedigree = mp;
  init();
  
  if(trait_valid(trait))
  {
    compute_pedigrees(trait);
  }
  else
  {
    compute_pedigrees();
  }
  
  if ( pedigree_count() <= 0 )
    return false;
    
  my_valid = true;
  return true;
}

// With simple trait.
//
void
MP_stats::compute_pedigrees(size_t trait)
{
  if(my_multipedigree != 0)
  {
    RefMultiPedigree::pedigree_const_iterator ped;

    for( ped = my_multipedigree->pedigree_begin(); ped != my_multipedigree->pedigree_end(); ++ped)
    {
      my_ped_stats.push_back( Ped_stats(&(*ped), trait, errors) );
      Ped_stats& ps = my_ped_stats.back();    

      if(!ps.valid())
      {
        my_ped_stats.pop_back();
        continue;
      }
     
      // Add the new pedigree General_Stats
      *this += ps;

      my_pedigree_size_info += ps.member_count();
      my_family_count_info  += ps.family_count();

      // Count 2n-f Frequencies
      my_likelihood_bits_freq.add(ps.likelihood_bits());
    
      // Count Family Size
      my_family_count_freq.add(ps.family_count());    

      // Count Generations
      my_generations_freq.add(ps.generations());
    }
  }
}

// With compound trait.
//
bool 
MP_stats::compute(const RefMultiPedigree* mp, std::vector<size_t> traits)
{
  if(mp == mped())      return true;

  my_multipedigree = mp;
  init();
  
  if(cmpd_trait_valid(traits))
  {
    compute_pedigrees(traits);
  }
  else
  {
    compute_pedigrees();
  }
  
  if ( pedigree_count() <= 0 )
    return false;
    
  my_valid = true;
  return true;
}

// With compound trait.
//
void
MP_stats::compute_pedigrees(std::vector<size_t> traits)
{
  if(my_multipedigree != 0)
  {
    RefMultiPedigree::pedigree_const_iterator ped;

    for( ped = my_multipedigree->pedigree_begin(); ped != my_multipedigree->pedigree_end(); ++ped)
    {
      my_ped_stats.push_back( Ped_stats(&(*ped), traits, errors) );
      Ped_stats& ps = my_ped_stats.back();    

      if(!ps.valid())
      {
        my_ped_stats.pop_back();
        continue;
      }
     
      // Add the new pedigree General_Stats
      *this += ps;

      my_pedigree_size_info += ps.member_count();
      my_family_count_info  += ps.family_count();

      // Count 2n-f Frequencies
      my_likelihood_bits_freq.add(ps.likelihood_bits());
    
      // Count Family Size
      my_family_count_freq.add(ps.family_count());    

      // Count Generations
      my_generations_freq.add(ps.generations());
    }
  }
}

std::string
MP_stats::compute_traits_name(std::vector<size_t> traits) const
{
  std::string temp = "";
  if(my_multipedigree != 0)
  {
    bool first = true;
    RefMPedInfo multi = my_multipedigree->info();
    size_t trait_cnt = multi.trait_count();
    for(size_t i = 0; i < traits.size(); ++i)
    {
      if(traits[i] < trait_cnt)
      { 
        RefTraitInfo rti = multi.trait_info(traits[i]);
        RefTraitInfo::trait_t t = rti.type();
        if(t == RefTraitInfo::continuous_trait ||
           t == RefTraitInfo::discrete_trait   ||
           t == RefTraitInfo::binary_trait        )
        {      
          if(!first)
          {
            temp += ", ";
          }
        
          temp += multi.trait_info(traits[i]).name();
          first = false;
          continue;
        }
      }
      
      if(!first)
      {
        temp += ", ";
      }
      
      temp += "???";
      first = false;
    } 
  }
  return temp;
}

//============================================================================
// IMPLEMENTATION:  Base_trait_stats
//============================================================================
//
void 
Base_trait_stats::init()
{
  for(int i = 0; i < 3; ++i)
  {
    my_sibships[i] = 0;
  }
  
  my_pedigree_size_info.clear();
  my_sibship_size_info.clear();
  my_family_count_freq.clear();
}

// - Tabulate sibships by the number of informative parents they have.
//
void 
Base_trait_stats::compute_sibships(const ind_filter& i_filter)
{
  if(my_pedigree == 0)
  {
    return;
  }
  
  int count;
  RefPedigree::family_iterator iter;
  for(iter = const_cast<RefPedigree*>(my_pedigree)->family_begin(); 
      iter != const_cast<RefPedigree*>(my_pedigree)->family_end(); ++iter)     // ***** const
  {
    count = 0;
    if(i_filter.informative(iter->parent1()))
    {
      ++count;
    }
    
    if(i_filter.informative(iter->parent2()))
    {
      ++count;
    }
    
    ++my_sibships[count];
  }
}

// - Determine size of pedigree counting only those people who have information
//   for the trait.
//
void
Base_trait_stats::compute_pedigree_size(const ind_filter& i_filter)
{
  assert(my_pedigree != 0);
  
  size_t  count = 0;
  RefPedigree::member_iterator iter;
  for(iter = const_cast<RefPedigree*>(my_pedigree)->member_begin();
      iter != const_cast<RefPedigree*>(my_pedigree)->member_end();            // ***** const
      ++iter)
  {
    if(i_filter.informative(&(*iter)))
    {
      ++count;
    }
  }
  
  my_pedigree_size_info.add(static_cast<double>(count)); 
}

// - determine number and sizes of the sibships in the pedigree counting
//   only those children who are informative for the trait in question.
//
void
Base_trait_stats::compute_sibship_size(const ind_filter& i_filter)
{
  RefMultiPedigree::family_iterator  f_iter = const_cast<RefPedigree*>(my_pedigree)->family_begin();
  for(; f_iter != const_cast<RefPedigree*>(my_pedigree)->family_end(); ++f_iter)
  {
    size_t  count = 0;
    RefMultiPedigree::offspring_iterator  c_iter = f_iter->offspring_begin();
    for(; c_iter != f_iter->offspring_end(); ++c_iter)
    {
      if(i_filter.informative(&(*c_iter)))
      {
        ++count;
      }
    }
    
    my_sibship_size_info.add(static_cast<double>(count));
  }
}

// - Count the number of nuclear families in the member pedigree pedigree.  
//   In this context a family is counted if at least one parent and at least 
//   one child is informative for the trait in question.
//
void                 
Base_trait_stats::compute_family_count(const ind_filter& i_filter)       // ***** const
{
  size_t  count = 0;
  RefMultiPedigree::family_iterator  f_iter = const_cast<RefPedigree*>(my_pedigree)->family_begin();
  for(; f_iter != const_cast<RefPedigree*>(my_pedigree)->family_end(); ++f_iter)
  {
    if(! (i_filter.informative(f_iter->parent1()) || 
          i_filter.informative(f_iter->parent2())   ))
    {
      continue;
    }
    else
    {
      RefMultiPedigree::offspring_iterator  c_iter = f_iter->offspring_begin();
      for(; c_iter != f_iter->offspring_end(); ++c_iter)
      {
        if(i_filter.informative(&(*c_iter)))
        {
          ++count;
          break;
        }
      }
    }
  }
  
  my_family_count_freq.add(static_cast<double>(count));
}

void
Base_trait_stats::print() const
{
  cout << "\nSibships w. no inf. parents  " << my_sibships[0];
  cout << "\nSibships w. one inf. parent   " << my_sibships[1];
  cout << "\nSibships w. two inf. parents  " << my_sibships[2] << endl;
}

Base_trait_stats& 
Base_trait_stats::operator+=(const Base_trait_stats& other)
{
    for(int i = 0; i < 3; ++i)
    {
      my_sibships[i] += other.my_sibships[i];
    }
    
    my_pedigree_size_info += other.my_pedigree_size_info;
    my_sibship_size_info += other.my_sibship_size_info;
    my_family_count_freq += other.my_family_count_freq;
    
  return *this;
}

// Convert from SAGE::sex to the enumeration, gender.
//
Base_trait_stats::gender
Base_trait_stats::determine_gender(MPED::SexCode ind_sex) const
{
  gender g;

  if(MPED::is_female(ind_sex))
  {  
    g = FEMALE;
  }
  else if(MPED::is_male(ind_sex))
  {
    g = MALE;
  }
  else                                           
  {
    g = UNKNOWN;
  }
  
  return g;
}

// - Determine if a member is a founder, non-founder, or unconnected.
//
Base_trait_stats::founder_status      
Base_trait_stats::determine_founder_status(member_const_pointer P1,
                                           member_const_pointer P2,
                                           size_t offspring_count) const
{
  founder_status f;

  if((!P1 || !P2) && offspring_count) 
  {
    f = FOUNDER;
  }
  else if(P1 && P2)
  {
    f = NON_FOUNDER; 
  }
  else
  {
    f = UNCONNECTED;
  }
  
  return f;   
}

//============================================================================
// IMPLEMENTATION:  Binary_trait_stats
//============================================================================
//
void                                       
Binary_trait_stats::compute()
{
  if(my_pedigree != 0)
  {
    RefMPedInfo multi = my_pedigree->multipedigree()->info();
    size_t trait_cnt = multi.trait_count();
    if(my_trait < trait_cnt)
    {
      RefTraitInfo rti = multi.trait_info(my_trait);
      if(rti.type() == RefTraitInfo::binary_trait)
      {
        my_trait_name = rti.name();
        compute_sibships();
        compute_pedigree_size();
        compute_sibship_size();
        compute_family_count();
        compute_inds_gender();
        compute_inds_founder();
        compute_pairs();
        my_valid = true;
      }
    }
  }
}

Binary_trait_stats& 
Binary_trait_stats::operator+=(const Binary_trait_stats& other)
{
  if(other.trait() == my_trait)
  {
    // Add sibship data.
    Base_trait_stats::operator+=(other);
  
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        my_ind_counts_gender[i][j] += other.my_ind_counts_gender[i][j];
      }
    }
    
    for(int i = 0; i < 3; ++i)
    {
      for(int j = 0; j < 3; ++j)
      {
        my_ind_counts_founder[i][j] += other.my_ind_counts_founder[i][j];
      }
    }
    
    for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
    {
      for(int j = 0; j < 4; ++j)
      {
        my_pair_counts[i][j] += other.my_pair_counts[i][j];
      }
    }

    for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
    {
      my_corr_infos[i] += other.my_corr_infos[i];
    }
  }
  return *this;
}
  
void          
Binary_trait_stats::init()
{
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      my_ind_counts_gender[i][j] = 0;
    }
  }
  
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < 3; ++j)
    {
      my_ind_counts_founder[i][j] = 0;
    }
  }
  
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    for(int j = 0; j < 4; ++j)
    {
      my_pair_counts[i][j] = 0;
    }
  }
}

// - Tabulate sibships according to number of informative parents: 0, 1, or 2.
//
void
Binary_trait_stats::compute_sibships()
{
  if(my_pedigree == 0)
  {
    return;
  }
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_sibships(filter);
}

// - Compute pedigree size counting only people with information for the trait.
//
void
Binary_trait_stats::compute_pedigree_size()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_pedigree_size(filter);
}

// - Determine number and size of the sibships in the pedigree counting only 
//   people with information for the trait in question.
//
void
Binary_trait_stats::compute_sibship_size()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_sibship_size(filter);
}

// - Count families in pedigree.  In this context a family must have at least
//   one informative parent and one informative child for the trait in question
//   to be counted.
//
void
Binary_trait_stats::compute_family_count()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_family_count(filter);
}

// - Tabulate individuals by gender and affection status.
//
void
Binary_trait_stats::compute_inds_gender()
{
  if(my_pedigree == 0)
  {
    return;
  }

  double                 v;             // Trait value.
  ind_affection_status   s;
  gender                 g;
  

  // For each individual, determine their state.

  RefPedigree:: member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    g = determine_gender(ind->get_effective_sex());

    // Determine if the trait is available.

    if(!my_pedigree->info().trait_missing(ind->index(), my_trait))
    {
      // If trait available, determine if affected or not.  We include a
      // third error state which should never occur, as other non-missing
      // values should have been removed at file read time.

      v = my_pedigree->info().trait(ind->index(), my_trait);

      if(v == 1.0)
      {
        s = AFF;
      }
      else if(v == 0.0)
      {
        s = UNAFF; 
      }
      else
      {
        SAGE_internal_error(); // Should never happen!
      }
    }
    else // Trait is missing
    {
      s = UNINF;
    }
    
    ++my_ind_counts_gender[g][s];
  }    
}

// - Tabulate individuals by founder and affection status.
//
void
Binary_trait_stats::compute_inds_founder()
{
  if(my_pedigree == 0)
  {
    return;
  }

  double                 v;             // Trait value.
  ind_affection_status   s;
  founder_status         f;
  
  // For each individual, determine their state.

  RefPedigree:: member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    f = determine_founder_status(ind->parent1(), ind->parent2(), ind->offspring_count());
    
    // Determine if the trait is available.

    if(!my_pedigree->info().trait_missing(ind->index(), my_trait))
    {
      // If trait available, determine if affected or not.  We include a
      // third error state which should never occur, as other non-missing
      // values should have been removed at file read time.

      v = my_pedigree->info().trait(ind->index(), my_trait);

      if(v == 1.0)
      {
        s = AFF;
      }
      else if(v == 0.0)
      {
        s = UNAFF; 
      }
      else
      {
        SAGE_internal_error(); // Should never happen!
      }
    }
    else // Trait is missing
    {
      s = UNINF;
    }
    
    ++my_ind_counts_founder[f][s];
  }    
}

// - Tabulate pairs by pair type and affection status and correlations
//   for each pair type.
//
inline void
Binary_trait_stats::compute_pairs()
{ 
  if(my_pedigree == 0)
  {
    return;
  }

  pg::pair_type             p_type;
  pair_filter_trait         pf_trait(my_trait);
  pair_filter               filter;
  pair_generator            gen(const_cast<RefPedigree*>(my_pedigree));  // ****** const
  filtering_pair_generator  fpg;
  
  fpg.set_pair_generator(gen);
  
  // - Pair counting.
  //
  for(int i = 1; i < 5; ++i)     // For each affection status.
  {
    pf_trait.set_status(static_cast<pft::affection_status>(i));  // Set to filter for one status at a time.
    filter.clear_traits();
    filter.add_trait(pf_trait);
    fpg.set_filter(filter);
    for(filtering_pair_generator::iterator iter = fpg.begin(); iter != fpg.end(); ++iter)
    {
      p_type = iter->type();
      ++my_pair_counts[p_type - 1][i - 1];
      
      // - Correlations.
      //
      if(static_cast<pft::affection_status>(i) != pft::UNINFORM)
      {
        RefPedigree::member_pointer member1 = iter->member_one();
        RefPedigree::member_pointer member2 = iter->member_two();
        
        compute_correlations(member1, member2, p_type);
      }
    }
  }
}

// - Add trait data to appropriate corinfo object.
//
inline void 
Binary_trait_stats::compute_correlations(RefPedigree::member_pointer member1, 
                                         RefPedigree::member_pointer member2, pg::pair_type p_type)
{
  RefPedInfo  ped_info  = my_pedigree->info();
  size_t      index1    = member1->index();
  size_t      index2    = member2->index();
  double      value1    = ped_info.trait(index1, my_trait);
  double      value2    = ped_info.trait(index2, my_trait);

  // Intraclass pair types.
  if(p_type == pg::SIBSIB  || p_type == pg::SISSIS || p_type == pg::BROBRO ||
     p_type == pg::HALFSIB || p_type == pg::COUSIN)
  {
    my_corr_infos[p_type - 1].add(value1, value2);
    my_corr_infos[p_type - 1].add(value2, value1);
  }
  // Interclass pair types.
  else if(p_type == pg::PARENTAL)
  {
    my_corr_infos[p_type - 1].add(index1 < index2 ? value1 : value2,
                                  index1 < index2 ? value2 : value1);
  }
  else if(p_type == pg::BROSIS)
  {
    if(member1->is_male())
    {
      my_corr_infos[p_type - 1].add(value1, value2);
    }
    else
    {
      my_corr_infos[p_type - 1].add(value2, value1);
    }
  }
  else if(p_type == pg::GRANDP)
  {
    my_corr_infos[p_type - 1].add(index1 < index2 ? value1 : value2,
                                  index1 < index2 ? value2 : value1);
  }
  else if(p_type == pg::AVUNC)
  {
    if(member1_is_uncle(index1, member2))
    {
      my_corr_infos[p_type - 1].add(value1, value2);
    }
    else
    {
      my_corr_infos[p_type - 1].add(value2, value1);
    }
  }
}

inline bool
Binary_trait_stats::member1_is_uncle(size_t index1, 
                                     RefPedigree::member_pointer member2)
{
  RefPedigree::member_pointer P1 = member2->parent1();
  RefPedigree::member_pointer P2 = member2->parent2();
  RefPedigree::sibling_iterator iter;
  for(iter = P1->sibling_begin(); iter != P1->sibling_end(); ++iter)
  {
    if(iter->index() == index1)
    {
      return true;
    }
  }
  
  for(iter = P2->sibling_begin(); iter != P2->sibling_end(); ++iter)
  {
    if(iter->index() == index1)
    {
      return true;
    }
  }
  
  return false;
}

// - Write object contents to standard output for testing purposes.
//
void
Binary_trait_stats::print() const
{
  Base_trait_stats::print();

  size_t total;

  // Individual gender data.
  cout << endl << endl;
  cout << setw(20) << left << "Trait:  ";
  cout << my_trait_name << endl;
  cout << setw(20) << left << "Males:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_gender[MALE][i];
    total += my_ind_counts_gender[MALE][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  cout << setw(20) << left << "Females:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_gender[FEMALE][i];
    total += my_ind_counts_gender[FEMALE][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  cout << setw(20) << left << "Unknowns:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_gender[UNKNOWN][i];
    total += my_ind_counts_gender[UNKNOWN][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  
  // Individual founder data.
  cout << endl << endl;
  cout << setw(20) << left << "Trait:  ";
  cout << my_trait_name << endl;
  cout << setw(20) << left << "Founders:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_founder[FOUNDER][i];
    total += my_ind_counts_founder[FOUNDER][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  cout << setw(20) << left << "Non-founder:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_founder[NON_FOUNDER][i];
    total += my_ind_counts_founder[NON_FOUNDER][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  cout << setw(20) << left << "Unconnecteds:";
  total = 0;
  for(int i = 0; i < 3; ++i)
  {
    cout << setw(10) << right << my_ind_counts_founder[UNCONNECTED][i];
    total += my_ind_counts_founder[UNCONNECTED][i];
  }
  cout << setw(10) << right << total;
  cout << endl;
  
  // Pair data.
  cout << endl;
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    cout << setw(20) << left << pg::pair_type_to_string(static_cast<pg::pair_type>(i + 1));
    total = 0;
    for(int j = 0; j < 4; ++j)
    {
      cout << setw(10) << right << my_pair_counts[i][j];
      total += my_pair_counts[i][j];
    }
    cout << setw(10) << right << total;
    cout << setw(20) << right << my_corr_infos[i].correlation() << "  ";
    cout << endl;
  }
  cout << endl << endl;
}

//============================================================================
// IMPLEMENTATION:  Cont_trait_stats
//============================================================================
//
void                                       
Cont_trait_stats::compute()
{
  if(my_pedigree != 0)
  {
    RefMPedInfo multi = my_pedigree->multipedigree()->info();
    size_t trait_cnt = multi.trait_count();
    if(my_trait < trait_cnt)
    {
      RefTraitInfo rti = multi.trait_info(my_trait);
      if(rti.type() == RefTraitInfo::continuous_trait ||
         rti.type() == RefTraitInfo::discrete_trait)
      {
        my_trait_name = rti.name();
        compute_sibships();
        compute_pedigree_size();
        compute_sibship_size();
        compute_family_count();
        compute_inds_gender();
        compute_inds_founder();
        compute_pairs();
        my_valid = true;
      }
    }
  }
}

Cont_trait_stats& 
Cont_trait_stats::operator+=(const Cont_trait_stats& other)
{
  if(other.trait() == my_trait)
  {
    // Add sibship data.
    Base_trait_stats::operator+=(other);
  
    for(int i = 0; i < 4; ++i)
    {
      my_inds_gender[i] += other.my_inds_gender[i];
      my_inds_founder[i] += other.my_inds_founder[i];
    }
    
    for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
    {
      my_pairs[i] += other.my_pairs[i];
      my_corr_infos[i] += other.my_corr_infos[i];
    }
  }
  return *this;
}

void          
Cont_trait_stats::init()
{
  for(int i = 0; i < 4; ++i)
  {
    my_inds_gender[i].clear();
    my_inds_founder[i].clear();
  }

  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    my_pairs[i].clear();  
  }
}

// - Tabulate sibships according to number of informative parents: 0, 1, or 2.
//
void
Cont_trait_stats::compute_sibships()
{
  if(my_pedigree == 0)
  {
    return;
  }
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_sibships(filter);
}

// - Compute pedigree size counting only people with information for the trait.
//
void
Cont_trait_stats::compute_pedigree_size()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_pedigree_size(filter);
}

// - Determine number and size of the sibships in the pedigree counting only 
//   people with information for the trait in question.
//
void
Cont_trait_stats::compute_sibship_size()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_sibship_size(filter);
}

// - Count families in pedigree.  In this context a family must have at least
//   one informative parent and one informative child for the trait in question
//   to be counted.
//
void
Cont_trait_stats::compute_family_count()
{
  assert(my_pedigree != 0);
  
  ind_filter_trait   if_trait(my_trait);
  ind_filter         filter(if_trait);
  Base_trait_stats::compute_family_count(filter);
}

// - Calculate individual trait data for each gender.
//
void
Cont_trait_stats::compute_inds_gender()
{
  if(my_pedigree == 0)
  {
    return;
  }

  double                 v;             // Trait value.
  gender                 g;
  
  RefPedigree:: member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    g = determine_gender(ind->get_effective_sex());
    
    v = my_pedigree->info().trait(ind->index(), my_trait);
    my_inds_gender[g].add(v);
    my_inds_gender[ALL].add(v);
  }    
}

// - Calculate individual trait data for each founder status.
//
void
Cont_trait_stats::compute_inds_founder()
{
  if(my_pedigree == 0)
  {
    return;
  }

  double                 v;             // Trait value.
  founder_status         f;
  
  RefPedigree:: member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    f = determine_founder_status(ind->parent1(), ind->parent2(), ind->offspring_count());
    
    v = my_pedigree->info().trait(ind->index(), my_trait);
    my_inds_founder[f].add(v);
    my_inds_founder[ALL].add(v);
  }    
}

// - Calculate pair trait data.
//
inline void
Cont_trait_stats::compute_pairs()
{ 
  if(my_pedigree == 0)
  {
    return;
  }

  pg::pair_type             p_type;
  pair_filter_trait         pf_trait(my_trait);   // Default filters on basis of informativity.
  
  pair_filter               p_filter(pf_trait);
  pair_generator            gen(const_cast<RefPedigree*>(my_pedigree));  // ****** const
  filtering_pair_generator  fpg(gen, p_filter);
  
  RefPedigree::member_pointer   member1;
  RefPedigree::member_pointer   member2;
  double                        v1;
  double                        v2;
  size_t                        index1;
  size_t                        index2;
  for(filtering_pair_generator::iterator iter = fpg.begin(); iter != fpg.end(); ++iter)
  {
    p_type     = iter->type();
    member1    = iter->member_one();
    member2    = iter->member_two();
    index1     = member1->index();
    index2     = member2->index();
    v1         = my_pedigree->info().trait(index1, my_trait);
    v2         = my_pedigree->info().trait(index2, my_trait);
    
    my_pairs[p_type - 1].add(v1);
    my_pairs[p_type - 1].add(v2);
    compute_correlations(member1, member2, index1, index2, v1, v2, p_type);
  }
}

// - Add trait data to appropriate corinfo object.
//
inline void 
Cont_trait_stats::compute_correlations(RefPedigree::member_pointer member1, 
                                       RefPedigree::member_pointer member2, 
                                       size_t index1, size_t index2,
                                       double value1, double value2, pg::pair_type p_type)
{
  // Intraclass pair types.
  if(p_type == pg::SIBSIB  || p_type == pg::SISSIS || p_type == pg::BROBRO ||
     p_type == pg::HALFSIB || p_type == pg::COUSIN)
  {
    my_corr_infos[p_type - 1].add(value1, value2);
    my_corr_infos[p_type - 1].add(value2, value1);
  }
  // Interclass pair types.
  else if(p_type == pg::PARENTAL)
  {
    my_corr_infos[p_type - 1].add(index1 < index2 ? value1 : value2,
                                  index1 < index2 ? value2 : value1);
  }
  else if(p_type == pg::BROSIS)
  {
    if(member1->is_male())
    {
      my_corr_infos[p_type - 1].add(value1, value2);
    }
    else
    {
      my_corr_infos[p_type - 1].add(value2, value1);
    }
  }
  else if(p_type == pg::GRANDP)
  {
    my_corr_infos[p_type - 1].add(index1 < index2 ? value1 : value2,
                                  index1 < index2 ? value2 : value1);
  }
  else if(p_type == pg::AVUNC)
  {
    if(member1_is_uncle(index1, member2))
    {
      my_corr_infos[p_type - 1].add(value1, value2);
    }
    else
    {
      my_corr_infos[p_type - 1].add(value2, value1);
    }
  }
}

inline bool
Cont_trait_stats::member1_is_uncle(size_t index1, 
                                     RefPedigree::member_pointer member2)
{
  RefPedigree::member_pointer P1 = member2->parent1();
  RefPedigree::member_pointer P2 = member2->parent2();
  RefPedigree::sibling_iterator iter;
  for(iter = P1->sibling_begin(); iter != P1->sibling_end(); ++iter)
  {
    if(iter->index() == index1)
    {
      return true;
    }
  }
  
  for(iter = P2->sibling_begin(); iter != P2->sibling_end(); ++iter)
  {
    if(iter->index() == index1)
    {
      return true;
    }
  }
  
  return false;
}

// - Write object contents to standard output for testing purposes.
//
void
Cont_trait_stats::print() const
{
  Base_trait_stats::print();

  // Individual gender data.
  cout << endl << endl;
  cout << setw(20) << left << "Trait:  ";
  cout << my_trait_name << endl;
  cout << setw(20) << left << "Males:";
  cout << setw(10) << right << ind_gender_count(MALE);
  cout << setw(10) << right << ind_gender_mean(MALE);
  cout << " +/- ";
  cout << setw(10) << right << ind_gender_std_dev(MALE);
  cout << setw(10) << right << ind_gender_min(MALE);
  cout << setw(10) << right << ind_gender_max(MALE);
  cout << endl;
  
  cout << setw(20) << left << "Females:";
  cout << setw(10) << right << ind_gender_count(FEMALE);
  cout << setw(10) << right << ind_gender_mean(FEMALE);
  cout << " +/- ";
  cout << setw(10) << right << ind_gender_std_dev(FEMALE);
  cout << setw(10) << right << ind_gender_min(FEMALE);
  cout << setw(10) << right << ind_gender_max(FEMALE);
  cout << endl;
  
  cout << setw(20) << left << "Unknown";
  cout << setw(10) << right << ind_gender_count(UNKNOWN);
  cout << setw(10) << right << ind_gender_mean(UNKNOWN);
  cout << " +/- ";
  cout << setw(10) << right << ind_gender_std_dev(UNKNOWN);
  cout << setw(10) << right << ind_gender_min(UNKNOWN);
  cout << setw(10) << right << ind_gender_max(UNKNOWN);
  cout << endl;
  
  cout << setw(20) << left << "All";
  cout << setw(10) << right << ind_gender_count(ALL);
  cout << setw(10) << right << ind_gender_mean(ALL);
  cout << " +/- ";
  cout << setw(10) << right << ind_gender_std_dev(ALL);
  cout << setw(10) << right << ind_gender_min(ALL);
  cout << setw(10) << right << ind_gender_max(ALL);
  cout << endl;
  
  cout << endl;
  
  // Individual founder data.
  cout << endl << endl;
  cout << setw(20) << left << "Trait:  ";
  cout << my_trait_name << endl;
  cout << setw(20) << left << "Founders:";
  cout << setw(10) << right << ind_founder_count(FOUNDER);
  cout << setw(10) << right << ind_founder_mean(FOUNDER);
  cout << " +/- ";
  cout << setw(10) << right << ind_founder_std_dev(FOUNDER);
  cout << setw(10) << right << ind_founder_min(FOUNDER);
  cout << setw(10) << right << ind_founder_max(FOUNDER);
  cout << endl;
  
  cout << setw(20) << left << "Non-founders:";
  cout << setw(10) << right << ind_founder_count(NON_FOUNDER);
  cout << setw(10) << right << ind_founder_mean(NON_FOUNDER);
  cout << " +/- ";
  cout << setw(10) << right << ind_founder_std_dev(NON_FOUNDER);
  cout << setw(10) << right << ind_founder_min(NON_FOUNDER);
  cout << setw(10) << right << ind_founder_max(NON_FOUNDER);
  cout << endl;
  
  cout << setw(20) << left << "Unconnected";
  cout << setw(10) << right << ind_founder_count(UNCONNECTED);
  cout << setw(10) << right << ind_founder_mean(UNCONNECTED);
  cout << " +/- ";
  cout << setw(10) << right << ind_founder_std_dev(UNCONNECTED);
  cout << setw(10) << right << ind_founder_min(UNCONNECTED);
  cout << setw(10) << right << ind_founder_max(UNCONNECTED);
  cout << endl;
  
  cout << setw(20) << left << "All";
  cout << setw(10) << right << ind_founder_count(ALL);
  cout << setw(10) << right << ind_founder_mean(ALL);
  cout << " +/- ";
  cout << setw(10) << right << ind_founder_std_dev(ALL);
  cout << setw(10) << right << ind_founder_min(ALL);
  cout << setw(10) << right << ind_founder_max(ALL);
  cout << endl;
  
  cout << endl;
  
  // Pair data.
  cout << endl;
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    cout << setw(20) << left << pg::pair_type_to_string(static_cast<pg::pair_type>(i + 1));
    cout << setw(10) << right << pair_count(static_cast<pg::pair_type>(i + 1));
    cout << setw(10) << right << pair_mean(static_cast<pg::pair_type>(i + 1));
    cout << " +/- ";
    cout << setw(10) << right << pair_std_dev(static_cast<pg::pair_type>(i + 1));
    cout << setw(10) << right << pair_min(static_cast<pg::pair_type>(i + 1));
    cout << setw(10) << right << pair_max(static_cast<pg::pair_type>(i + 1));
    
    cout << setw(20) << right << my_corr_infos[i].correlation() << "  ";
    cout << endl;
  }
  cout << endl << endl;
}

//============================================================================
// IMPLEMENTATION:  Cmpd_trait_stats
//============================================================================
//
Cmpd_trait_stats& 
Cmpd_trait_stats::operator+=(const Cmpd_trait_stats& other)
{
  if(my_traits == other.my_traits)
  {
    // Add sibship data.
    Base_trait_stats::operator+=(other);  
  
    for(int i = 0; i < 3; ++i)
    {
      my_ind_counts_gender[i] += other.my_ind_counts_gender[i];
      my_ind_counts_founder[i] += other.my_ind_counts_founder[i]; 
    }
    
    for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
    {
      my_pair_counts[i] += other.my_pair_counts[i];
    }
  }
  return *this;
}

void
Cmpd_trait_stats::compute()
{
  if(my_pedigree != 0)
  {
    compute_traits_name();
    compute_sibships();
    compute_pedigree_size();
    compute_sibship_size();
    compute_family_count();
    compute_inds_gender();
    compute_inds_founder();
    compute_pairs();
    my_valid = true;
  }
}

void
Cmpd_trait_stats::compute_traits_name()
{
  if(my_pedigree != 0)
  {
    bool first = true;
    my_traits_name = "";
    RefMPedInfo multi = my_pedigree->multipedigree()->info();
    size_t trait_cnt = multi.trait_count();
    for(size_t i = 0; i < my_traits.size(); ++i)
    {
      if(my_traits[i] < trait_cnt)
      { 
        RefTraitInfo rti = multi.trait_info(my_traits[i]);
        RefTraitInfo::trait_t t = rti.type();
        if(t == RefTraitInfo::continuous_trait ||
           t == RefTraitInfo::discrete_trait   ||
           t == RefTraitInfo::binary_trait        )
        {      
          if(!first)
          {
            my_traits_name += ", ";
          }
        
          my_traits_name += multi.trait_info(my_traits[i]).name();
          first = false;
          continue;
        }
      }
      
      if(!first)
      {
        my_traits_name += ", ";
      }
      
      my_traits_name += "???";
      first = false;
    } 
  }
}

void          
Cmpd_trait_stats::init()
{
  for(int i = 0; i < 3; ++i)
  {
    my_ind_counts_gender[i] = 0;
    my_ind_counts_founder[i] = 0;
  }
  
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    my_pair_counts[i] = 0;
  }
  
  make_ind_filter();
}

void
Cmpd_trait_stats::make_ind_filter()
{
  ind_filter_trait    if_trait;
  my_ind_filter.clear_traits();
  for(size_t t = 0; t < my_traits.size(); ++t)
  {
    if_trait.set_trait(my_traits[t]);
    my_ind_filter.add_trait(if_trait);
  }
}

// - Tabulate sibships according to number of informative parents: 0, 1, or 2.
//
void
Cmpd_trait_stats::compute_sibships()
{
  if(my_pedigree == 0)
  {
    return;
  }
  
  Base_trait_stats::compute_sibships(my_ind_filter);
}

// - Compute pedigree size counting only people with information for the trait.
//
void
Cmpd_trait_stats::compute_pedigree_size()
{
  assert(my_pedigree != 0);
  
  Base_trait_stats::compute_pedigree_size(my_ind_filter);
}

// - Determine number and size of the sibships in the pedigree counting only 
//   people with information for the trait in question.
//
void
Cmpd_trait_stats::compute_sibship_size()
{
  assert(my_pedigree != 0);
  
  Base_trait_stats::compute_sibship_size(my_ind_filter);
}

// - Count families in pedigree.  In this context a family must have at least
//   one informative parent and one informative child for the trait in question
//   to be counted.
//
void
Cmpd_trait_stats::compute_family_count()
{
  assert(my_pedigree != 0);
  
  Base_trait_stats::compute_family_count(my_ind_filter);
}

// - Tabulate informative individuals by gender.
//
void
Cmpd_trait_stats::compute_inds_gender()
{
  if(my_pedigree == 0)
  {
    return;
  }
  
  gender       g;
  
  RefPedigree::member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    g = determine_gender(ind->get_effective_sex());
    
    if(my_ind_filter.informative(const_cast<RefPedigree::member_pointer>(&(*ind))))    // const ******
    
    
    {
      ++my_ind_counts_gender[g];
    }
  }    
}

// - Tabulate informative individuals by founder/nonfounder/unconnected catagories.
//
void
Cmpd_trait_stats::compute_inds_founder()
{
  if(my_pedigree == 0)
  {
    return;
  }

  founder_status       f;
  
  RefPedigree::member_const_iterator  ind;

  for(ind = my_pedigree->member_begin(); ind != my_pedigree->member_end(); ++ind)
  {
    f = determine_founder_status(ind->parent1(), ind->parent2(), ind->offspring_count());
      
    if(my_ind_filter.informative(const_cast<RefPedigree::member_pointer>(&(*ind))))       // ****** const
    {
      ++my_ind_counts_founder[f];
    }
  }    
}

// Tabulate informative pairs by type.
//
void 
Cmpd_trait_stats::compute_pairs()
{
  if(my_pedigree == 0)
  {
    return;
  }

  // Make pair filter.
  pair_filter_trait      pf_trait;
  pair_filter   filter;
  filter.clear_traits();
  for(size_t t = 0; t < my_traits.size(); ++t)
  {
    pf_trait.set_trait(my_traits[t]);
    filter.add_trait(pf_trait);
  }

  // Count pairs.
  pg::pair_type             p_type;
  pair_generator            gen(const_cast<RefPedigree*>(my_pedigree));  // ****** const
  filtering_pair_generator  fpg(gen, filter);
  for(filtering_pair_generator::iterator iter = fpg.begin(); iter != fpg.end(); ++iter)
  {
    p_type = iter->type();
    ++my_pair_counts[p_type - 1];
  }
}

// - Write object contents to standard output for testing purposes.
//
void
Cmpd_trait_stats::print() const
{
  Base_trait_stats::print();

  // Individual data.
  cout << endl << endl;
  cout << setw(20) << left << "Trait:  ";
  cout << my_traits_name << endl;
  
  cout << setw(20) << left << "Males:";
  cout << setw(10) << right << my_ind_counts_gender[MALE];
  cout << endl;
  
  cout << setw(20) << left << "Females:";
  cout << setw(10) << right << my_ind_counts_gender[FEMALE];
  cout << endl;
  
  cout << setw(20) << left << "Unknowns:";
  cout << setw(10) << right << my_ind_counts_gender[UNKNOWN];
  cout << endl;
  
  cout << endl;
  cout << setw(20) << left << "Founders:";
  cout << setw(10) << right << my_ind_counts_founder[FOUNDER];
  cout << endl;
  
  cout << setw(20) << left << "Non Founder:";
  cout << setw(10) << right << my_ind_counts_founder[NON_FOUNDER];
  cout << endl;
  
  cout << setw(20) << left << "Unconnected:";
  cout << setw(10) << right << my_ind_counts_founder[UNCONNECTED];
  cout << endl;
  
  // Pair data.
  cout << endl;
  for(int i = 0; i < TOTAL_PAIR_TYPES; ++i)
  {
    cout << setw(20) << left << pg::pair_type_to_string(static_cast<pg::pair_type>(i + 1));
    cout << setw(10) << right << my_pair_counts[i];
    cout << endl;
  }
  cout << endl << endl;
}


} // End namespace RPED
} // End namespace SAGE

