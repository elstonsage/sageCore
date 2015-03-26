//============================================================================
// File:      relpair.cpp                      
//                                                                          
// Author:    Dan Baechle & Kevin Jacobs                                    
//                                                                          
// History:   7/00   created.  - djb
//            4/4/01 modified so that relative_pair class includes connecting 
//                   members.  -djb
//                                                                          
// Notes:     Non-inline implementation for the following classes -    
//              pair_generator (and iterator classes) 
//              ind_filtering_trait
//              pair_filtering_trait
//              ind_filter
//              pair_filter
//              filtering_pair_generator (and iterator classes)
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "pairs/relpair.h"

namespace SAGE {

//============================================================================
// IMPLEMENTATION:  pair_generator::base_iterator      
//============================================================================
//
// - Initialize state variables for a family.
//
void
pair_generator::base_iterator::init_family()
{
  if(my_family == my_generator->my_p->family_end())
  {
    at_end = true;
    return;
  }
  else
  {
    my_parent             = my_family->parent1();
    my_grandp             = my_parent->parent1();
    my_children           = my_family->offspring_begin();
    my_fullsibs           = my_family->offspring_begin();
    ++my_fullsibs;
    my_parents_mates      = my_parent->mate_begin();
    my_avuncs             = my_parent->sibling_begin();
    halfsibs_init         = false;
    cousins_init          = false;                        
    
    my_parent_number      = P1;
    my_grandp_number      = GP1;
  }
}

// - my_subped just initialized or incremented.
//
void
pair_generator::base_iterator::init_every()
{
  if(my_subped == my_generator->my_p->subpedigree_end())
  {
    at_end = true;
    return;
  }
  else
  {
    my_member_one = my_subped->member_begin();
    my_member_two = my_member_one;
    ++my_member_two;
  }
}

void
pair_generator::base_iterator::seek()
{
  my_current_type = my_generator->next_type(my_current_type);
  if(my_current_type != my_generator->my_first_type)           // Pair type.
  {
    my_children      = my_family->offspring_begin();
    my_fullsibs      = my_family->offspring_begin();
    ++my_fullsibs;
    my_avuncs        = my_parent->sibling_begin();
    my_parents_mates = my_parent->mate_begin();
  }  
  else if(my_parent_number == P1)                              // Parent.
  {
    init_family();
    my_parent         = my_family->parent2();
    my_parent_number  = P2;
    my_grandp         = my_parent->parent1();
    my_parents_mates  = my_parent->mate_begin();       
    my_avuncs         = my_parent->sibling_begin();    
  }
  else                                                         // Family.
  {
    ++(my_family);
    init_family();
  }
}

// - Called when member_two comes to an end.
//
void
pair_generator::base_iterator::seek_every()
{
  ++my_member_one;
  if(my_member_one != my_subped->member_end())
  {
    my_member_two = my_member_one;
    ++my_member_two;
    if(my_member_two != my_subped->member_end())
    {
      return;
    }
    else
    {
      seek_every();
    }
  }
  else
  {
    ++my_subped;
    init_every();
  }
}

void
pair_generator::base_iterator::set_every_pair()
{
  my_pair.my_member_one     = &(*my_member_one);
  my_pair.my_member_two     = &(*my_member_two);
  my_pair.my_connector_one  = 0;
  my_pair.my_connector_two  = 0;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

void 
pair_generator::base_iterator::every()
{
  if(my_member_two != my_subped->member_end())
  {
    set_every_pair();
    
    // Point to next member two
    ++(my_member_two);
  }
  else
  {
    seek_every();  
  }
}

// - Assign values to parental pair. 
//  
void
pair_generator::base_iterator::set_parental_pair()
{
  my_pair.my_member_one     = my_parent;
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = 0;
  my_pair.my_connector_two  = 0;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

// - Find parental pairs.
//
void 
pair_generator::base_iterator::parental()
{
  if(my_children != my_family->offspring_end())
  {
    set_parental_pair();
    
    // Point to next child in parental relationship
    ++(my_children);
  }
  else
  {
    seek();  
  }
}

// - Increment my_fullsibs and or my_children so that they point
//   to the next candidate pair.
//
void
pair_generator::base_iterator::incr_sibs()
{
  if(my_fullsibs != my_family->offspring_end())
  {
    ++my_fullsibs;
  }
  else
  {
    ++my_children;
    if(my_children != my_family->offspring_end())
    {
      // Keep my_fullsibs "one ahead" of my_children to avoid duplicate pairs.
      my_fullsibs = my_children;
      ++my_fullsibs;
    }
  }
}

// - Assign values to pair formed from a single sibship.
//  
void
pair_generator::base_iterator::set_sibling_pair()
{
  my_pair.my_member_one     = &(*my_fullsibs);
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = my_family->parent1();
  my_pair.my_connector_two  = my_family->parent2();

  my_pair.my_type = my_current_type;
  pair_valid = true;
}


// - Find sibling - sibling pairs.
//
void 
pair_generator::base_iterator::sibsib()
{
  if(my_parent_number == P2)                            // All sibsib pairs found for this family.
  {
    seek();
  }
  else if(my_fullsibs != my_family->offspring_end())    // Valid pair.
  {
    set_sibling_pair();  
    incr_sibs();
  }
  else if(my_children != my_family->offspring_end())    // More pairs to be found for this family.
  {
    incr_sibs();
    sibsib();
  }
  else
  {
    seek();
  }
}

// - Find sister - sister pairs.
//
void 
pair_generator::base_iterator::sissis()
{
  if(my_parent_number == P2)            // All sissis pairs found for this family.
  {
    seek();
  }
  else if(my_fullsibs != my_family->offspring_end())
  {
    if(!(my_fullsibs->is_female() &&  my_children->is_female()))
    { 
      incr_sibs();
      sissis();
    }
    else
    {
      set_sibling_pair();
      incr_sibs();
    }
  }
  else if(my_children != my_family->offspring_end())
  {
    incr_sibs();
    sissis();
  }
  else
  {
    seek();
  }
}

// - Find brother - brother pairs.
//
void 
pair_generator::base_iterator::brobro()
{
  if(my_parent_number == P2)    // All brobro pairs have been found for this family.
  {
    seek();
  }
  else if(my_fullsibs != my_family->offspring_end())
  {
    if(!(my_fullsibs->is_male() &&  my_children->is_male()))
    { 
      incr_sibs();
      brobro();
    }
    else
    {
      set_sibling_pair();
      incr_sibs();
    }
  }
  else if(my_children != my_family->offspring_end())
  {
    incr_sibs();
    brobro();
  }
  else
  {
    // Look for next pair.
    seek();
  }
}

// - Find brother - sister pairs.
//
void 
pair_generator::base_iterator::brosis()
{
  if(my_parent_number == P2)    // All brosis pairs have been found for this family.
  {
    seek();
  }
  else if(my_fullsibs != my_family->offspring_end())
  {
    if(!( (my_fullsibs->is_male()   &&   my_children->is_female() ) ||
          (my_fullsibs->is_female() &&   my_children->is_male()   )    ))
    { 
      incr_sibs();
      brosis();
    }
    else
    {
      set_sibling_pair();
      incr_sibs();
    }
  }
  else if(my_children != my_family->offspring_end())
  {
    incr_sibs();
    brosis();
  }
  else
  {
    // Look for next pair.
    seek();
  }
}

// - Assign values to grandparental pairs.
//
void
pair_generator::base_iterator::set_grandp_pair()
{
  my_pair.my_member_one     = my_grandp;
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = my_parent;
  my_pair.my_connector_two  = 0;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

// - Increment grandparent or children and grandparent so that they point
//   to the next candidate pair.
//
void
pair_generator::base_iterator::incr_grandps()
{
  if(my_grandp_number == GP1)
  {
    my_grandp = my_parent->parent2();
    my_grandp_number = GP2;
  }
  else
  {
    ++my_children;
    if(my_children != my_family->offspring_end())
    {
      my_grandp = my_parent->parent1();
      my_grandp_number = GP1;
    }
  }
}

// - Find grandparental pairs.
//
void 
pair_generator::base_iterator::grandp()
{
  if(my_children != my_family->offspring_end())
  { 
    if(my_grandp)
    {
      set_grandp_pair();  
      incr_grandps();
    }
    else if(my_grandp_number == GP1)
    {
      incr_grandps();
    }
    else
    {
      seek();
    }
  }
  else
  {
    seek();
  }
}

// - increment uncles and aunts, or children and uncles and aunts so that they
//   point to the next candidate pair. 
//
void
pair_generator::base_iterator::incr_avuncs()
{
  if(my_avuncs != my_parent->sibling_end())
  {
    ++my_avuncs;
  }
  else
  {
    ++my_children;
    if(my_children != my_family->offspring_end())
    {
      my_avuncs = my_parent->sibling_begin();
    }
  }
}

// - Assign values to avuncular pair.
//  
void
pair_generator::base_iterator::set_avunc_pair()
{
  my_pair.my_member_one     = &(*my_avuncs);
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = my_parent;
  my_pair.my_connector_two  = 0;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

// - Find avuncular pairs.
//
void 
pair_generator::base_iterator::avunc()
{
  if(my_avuncs != my_parent->sibling_end())
  {
    set_avunc_pair();  
    incr_avuncs();
  }
  else if(my_children != my_family->offspring_end())
  {
    incr_avuncs();
    avunc();
  }
  else
  {
    seek();
  }
}

// - Assign values to half sibling - half sibling pair.
//  
void
pair_generator::base_iterator::set_halfsib_pair()
{
  my_pair.my_member_one     = &(*my_halfsibs);
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = my_parent;
  my_pair.my_connector_two  = 0;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

// - Increment parents mates and reset halfsibs and children.
//  
void
pair_generator::base_iterator::incr_parents_mates_hs()
{
  if(my_parents_mates == my_parent->mate_end())
  {
    return;
  }
  else
  {
    ++my_parents_mates;
    if(my_parents_mates == my_parent->mate_end())
    {
      return;
    }
    else if(other_parent(&(my_parents_mates->mate())))
    {
      incr_parents_mates_hs();
    }
    else
    {
      my_halfsibs = my_parent->offspring_begin(my_parents_mates);
      my_children = my_family->offspring_begin();
    }
  }
}

// - Increment children and reset halfsibs.
//  
void
pair_generator::base_iterator::incr_children_hs()
{
  if(my_children == my_family->offspring_end())
  {
    incr_parents_mates_hs();
  }
  else
  {
    ++my_children;
    if(my_children == my_family->offspring_end())
    {
      incr_parents_mates_hs();
    }
    else
    {
      my_halfsibs = my_parent->offspring_begin(my_parents_mates);
    }
  }
}

// - Advance iterators to next candidate pair if one exists. 
//
void
pair_generator::base_iterator::incr_halfsibs_hs()
{
  if(my_halfsibs == my_parent->offspring_end())
  {
    incr_children_hs();
  }
  else
  {
    ++my_halfsibs;
  }
}
  
// - Find initial canditate pair if one exists.
//  
void
pair_generator::base_iterator::init_halfsibs()
{
  if(my_parents_mates == my_parent->mate_end())       // No more halfsibs to be found for this parent.
  {
    halfsibs_init = true;
  }
  else if(other_parent(&(my_parents_mates->mate())))  // Mate is other parent in current nuclear family.
  {
    ++my_parents_mates;
    init_halfsibs();
  }
  else
  {
    my_halfsibs = my_parent->offspring_begin(my_parents_mates);
    halfsibs_init = true;
  }
}

// - Is parent's mate the other parent in the current nuclear family?
//  
bool
pair_generator::base_iterator::other_parent(RPED::RefPedigree::member_pointer mate)
{
  return mate == my_family->parent1() || mate == my_family->parent2();
}

// - Find half sibling - half sibling pairs.
//
void 
pair_generator::base_iterator::halfsib()
{
  if(!halfsibs_init)
  {
    init_halfsibs();
  }
  
  if(my_parents_mates == my_parent->mate_end())        // Done searching for halfsibs for this parent.
  {
    halfsibs_init = false;
    seek();
  }
  else if(my_halfsibs == my_parent->offspring_end())   // Valid halfsib not found.
  {
    incr_halfsibs_hs();
    halfsib();
  }
  else if(&(*my_family) > my_halfsibs->family())       // Halfsib pair not found elsewhere.
  {
    set_halfsib_pair();
    incr_halfsibs_hs();
  }
  else
  {
    incr_halfsibs_hs();
    halfsib();
  }
}

// - Assign values to cousin - cousin pair.
//  
void
pair_generator::base_iterator::set_cousin_pair()
{
  my_pair.my_member_one     = &(*my_cousins);
  my_pair.my_member_two     = &(*my_children);
  my_pair.my_connector_one  = &(*my_avuncs);
  my_pair.my_connector_two  = my_parent;
  
  my_pair.my_type = my_current_type;
  pair_valid = true;
}

// - Increment uncles/aunts and reset mates, cousins and children.
//  
void
pair_generator::base_iterator::incr_avuncs_c()
{
  if(my_avuncs == my_parent->sibling_end())
  {
    return;
  }
  else
  {
    ++my_avuncs;
    if(my_avuncs == my_parent->sibling_end())
    { 
      return;
    }
    else
    {
      my_avuncs_mates = my_avuncs->mate_begin();
      if(my_avuncs_mates == my_avuncs->mate_end())    // This uncle/aunt has no mates.
      {
        incr_avuncs_c();
      }
      else
      {
        my_cousins = my_avuncs->offspring_begin(my_avuncs_mates);
        my_children = my_family->offspring_begin();
      }
    }
  }
}

// - Increment uncles/aunts mates and reset cousins and children.
//  
void
pair_generator::base_iterator::incr_avuncs_mates_c()
{
  if(my_avuncs_mates == my_avuncs->mate_end())
  {
    incr_avuncs_c();
  }
  else
  {
    ++my_avuncs_mates;
    if(my_avuncs_mates == my_avuncs->mate_end())
    {
      incr_avuncs_c();
    }
    else
    {
      my_cousins = my_avuncs->offspring_begin(my_avuncs_mates);
      my_children = my_family->offspring_begin();
    }
  }
}

// - Inrement children and reset cousins.
//  
void
pair_generator::base_iterator::incr_children_c()
{
  if(my_children == my_family->offspring_end())
  {
    incr_avuncs_mates_c();
  }
  else
  {
    ++my_children;
    if(my_children == my_family->offspring_end())
    {
      incr_avuncs_mates_c();
    }
    else
    {
      my_cousins = my_avuncs->offspring_begin(my_avuncs_mates);
    }
  }
}

// - Advance iterators to next candidate cousin pair if one exists.
//  
void
pair_generator::base_iterator::incr_cousins_c()
{
  if(my_cousins == my_avuncs->offspring_end())   
  {
    incr_children_c();
  }
  else
  {
    ++my_cousins;
  }
}

// - Advance iterators to first candidate cousin pair if one exists. 
//  
void
pair_generator::base_iterator::init_cousins()
{
  if(my_avuncs == my_parent->sibling_end())         // No cousins thru this parent.
  {
    cousins_init = true;
  }
  else
  {
    my_avuncs_mates = my_avuncs->mate_begin();
    if(my_avuncs_mates == my_avuncs->mate_end())    // This uncle/aunt has no mates.
    {
      ++my_avuncs;
      init_cousins();
    }
    else
    {
      my_cousins = my_avuncs->offspring_begin(my_avuncs_mates);
      cousins_init = true;
    }
  }
}

// - Find cousin - cousin pairs.
//
void 
pair_generator::base_iterator::cousin()
{
  if(!cousins_init)
  {
    init_cousins();
  }
  
  if(my_avuncs == my_parent->sibling_end())            // Cousin search done for this parent.
  {
    seek();
    cousins_init = false;
  }
  else if(my_cousins == my_avuncs->offspring_end())    // Valid cousin pair not found.
  {
    incr_cousins_c();
    cousin();
  }
  else if(&(*my_family) > my_cousins->family())        // Cousin pair not found elsewhere.
  {
    set_cousin_pair();
    incr_cousins_c();
  }
  else
  {
    incr_cousins_c();
    cousin();
  }
}

// - Invalid type or null type.  Set iterator to end.
//  
void
pair_generator::base_iterator::null_type()
{
  at_end = true;
}

//============================================================================
// IMPLEMENTATION:  pair_generator::iterator
//============================================================================
//
pair_generator::iterator::iterator(const pair_generator* generator,
                                   RPED::RefPedigree::family_iterator family) 
      : base_iterator(generator, family)
{
  my_current_type = my_generator->my_first_type;

  if(my_family == my_generator->my_p->family_end() ||
     my_current_type == pair_generator::NULL_TYPE    )
  {
    at_end = true;
    return;
  }
  else
  { 
    at_end = false;
    if(my_current_type == pair_generator::EVERY)
    {
      init_every();
    }
    else
    {
      init_family();
    }
    
    ++(*this);
  }
}

// - Traverse pedigree by family looking for pairs of the specified types in
//   a predefined order.  When a pair is found, set data members of the 
//   a temporary pair, and return.
//
pair_generator::iterator*
pair_generator::iterator::operator ++()
{
  pair_valid = false;
  while(!pair_valid && !at_end)
  {
    switch(my_current_type) 
    {
      case EVERY:
        every();
        break;
      case PARENTAL:
        parental();
        break;
      case SIBSIB:
        sibsib();
        break;
      case SISSIS:
        sissis();
        break;
      case BROBRO:
        brobro();
        break;
      case BROSIS:
        brosis();
        break;
      case GRANDP:
        grandp();
        break;  
      case AVUNC:
        avunc();
        break;
      case HALFSIB:
        halfsib();
        break;
      case COUSIN:
        cousin();
        break;
      default:
        null_type(); 
    } 
  }
  return this;
}

//============================================================================
// IMPLEMENTATION:  pair_generator::const_iterator
//============================================================================
//
pair_generator::const_iterator::const_iterator(const pair_generator* generator,
                                   RPED::RefPedigree::family_iterator family) 
      : base_iterator(generator, family)
{
  if(my_family == my_generator->my_p->family_end())
  {
    at_end = true;
    return;
  }
  else
  { 
    at_end = false;
    my_current_type = my_generator->my_first_type;
    init_family();
    ++(*this);
  }
}

// - Traverse pedigree by family looking for pairs of the specified types in
//   a predefined order.  When a pair is found, set data members of the 
//   a temporary pair, and return.
//
pair_generator::const_iterator*
pair_generator::const_iterator::operator ++()
{
  pair_valid = false;
  while(!pair_valid && !at_end)
  {
    switch(my_current_type) 
    {
      case PARENTAL:
        parental();
        break;
      case SIBSIB:
        sibsib();
        break;
      case SISSIS:
        sissis();
        break;
      case BROBRO:
        brobro();
        break;
      case BROSIS:
        brosis();
        break;
      case GRANDP:
        grandp();
        break;  
      case AVUNC:
        avunc();
        break;
      case HALFSIB:
        halfsib();
        break;
      case COUSIN:
        cousin();
        break;
      default:
        null_type(); 
    } 
  }
  return this;
}

//============================================================================
// IMPLEMENTATION:  ind_filter_trait
//============================================================================
//
// - Determine whether a pedigree member is informative for this trait. 
//  
bool 
ind_filter_trait::informative(const RPED::RefPedigree::member_pointer member) const
{
  if(member != 0)
  { 
    RPED::RefMPedInfo rmpi = member->multipedigree()->info();
    size_t t_count = rmpi.trait_count();
    if(my_trait < t_count)
    {
      RPED::RefTraitInfo::trait_t trait_type = rmpi.trait_info(my_trait).type();

      if(   trait_type == RPED::RefTraitInfo::binary_trait 
         || trait_type == RPED::RefTraitInfo::continuous_trait
         || trait_type == RPED::RefTraitInfo::discrete_trait
         || trait_type == RPED::RefTraitInfo::categorical_trait)
      {
        return !SAGE::isnan(member->pedigree()->info().trait(member->index(), my_trait));
      }
    }
  }
  return false;
}

// - Determine whether trait value for a pedigree member is between min and max.
//  
bool
ind_filter_trait::in_range(const RPED::RefPedigree::member_pointer member) const
{
  if(informative(member))
  {
    double value = member->pedigree()->info().trait(member->index(), my_trait);
    if(value >= my_min && value <= my_max)                            
    {
      return true;
    }
  }
  return false;
}

// - Determine whether a pedigree member is affected for this trait.
//  
bool
ind_filter_trait::affected(const RPED::RefPedigree::member_pointer member) const
{
  if(in_range(member))
  {
    RPED::RefTraitInfo::trait_t trait_type = member->multipedigree()->info().trait_info(my_trait).type();
    
    // Binary.
    if(trait_type == RPED::RefTraitInfo::binary_trait)
    {
      if(member->pedigree()->info().trait(member->index(), my_trait) == 1.0)
      {
        return true;
      }
    }
    else
    { 
      // Continuous or discrete.
      double trait_value = member->pedigree()->info().trait(member->index(), my_trait);
      if(trait_value > my_affected_min && trait_value < my_affected_max)    
      {
        return true;
      }
    }
  } 
  return false; 
}

// - Determine whether a pedigree member is unaffected for this trait.
//  
bool
ind_filter_trait::unaffected(const RPED::RefPedigree::member_pointer member) const
{
  if(in_range(member))
  {
    RPED::RefTraitInfo::trait_t trait_type = member->multipedigree()->info().trait_info(my_trait).type();
    
    // Binary.
    if(trait_type == RPED::RefTraitInfo::binary_trait)
    {
      if(member->pedigree()->info().trait(member->index(), my_trait) == 0.0)
      {
        return true;
      }
    }
    else
    { 
      // Continuous or discrete.
      double trait_value = member->pedigree()->info().trait(member->index(), my_trait);
      if(trait_value >= my_unaffected_min && trait_value <= my_unaffected_max)   
      {
        return true;
      }
    }
  } 
  return false; 
}

//============================================================================
// IMPLEMENTATION:  pair_filter_trait
//============================================================================
//
// - Determine validity of a pair for a single trait.
//
bool 
pair_filter_trait::valid(const pair_generator::relative_pair& pair) const
{
  // Don't filter if no valid affection status is specified.
  if(!valid_status)
  {
    return true;
  }

  if(my_status & UNINFORM_MASK)
  {
    if(!informative(pair.member_one()) || !informative(pair.member_two()))
    {
      return true;
    }
  }

  if(my_status & CONCORD_UNAFF_MASK)
  {
    if(unaffected(pair.member_one()) && unaffected(pair.member_two()))
    {
      return true;
    }
  }
  
  if(my_status & DISCORD_MASK)
  {
    if(    (unaffected(pair.member_one()) &&   affected(pair.member_two()))
        || (  affected(pair.member_one()) && unaffected(pair.member_two())) )
    {
      return true;
    }
  }
  
  if(my_status & CONCORD_AFF_MASK)
  {
    if(affected(pair.member_one()) && affected(pair.member_two()))
    {
      return true;
    }
  }
  
  return false;
}

//============================================================================
// IMPLEMENTATION:  ind_filter
//============================================================================
//
// - To be informative, an individual must be informative for all 
//   ind_filter_traits in the ind_filter.
//
bool
ind_filter::informative(const RPED::RefPedigree::member_pointer member) const
{
  std::list<ind_filter_trait>::const_iterator trait;
  for(trait = my_traits.begin(); trait != my_traits.end(); ++trait)
  {
    if(!trait->informative(member))
    {
      return false;
    }
  }
  return true;
}

// - To be in range, an individual must be in range for all 
//   ind_filter_traits in the ind_filter.
//
bool
ind_filter::in_range(const RPED::RefPedigree::member_pointer member) const
{
  std::list<ind_filter_trait>::const_iterator trait;
  for(trait = my_traits.begin(); trait != my_traits.end(); ++trait)
  {
    if(!trait->in_range(member))
    {
      return false;
    }
  }
  return true;
}

// - To be affected, an individual must be affected for all 
//   ind_filter_traits in the ind_filter.
//
bool
ind_filter::affected(const RPED::RefPedigree::member_pointer member) const
{
  std::list<ind_filter_trait>::const_iterator trait;
  for(trait = my_traits.begin(); trait != my_traits.end(); ++trait)
  {
    if(!trait->affected(member))
    {
      return false;
    }
  }
  return true;
}

// - To be unaffected, an individual must be unaffected for all 
//   ind_filter_traits in the ind_filter.
//
bool
ind_filter::unaffected(const RPED::RefPedigree::member_pointer member) const
{
  std::list<ind_filter_trait>::const_iterator trait;
  for(trait = my_traits.begin(); trait != my_traits.end(); ++trait)
  {
    if(!trait->unaffected(member))
    {
      return false;
    }
  }
  return true;
}


//============================================================================
// IMPLEMENTATION:  pair_filter
//============================================================================
//
// - To be valid, a pair must be valid for all pair_filter_traits in the pair_filter.
//
bool
pair_filter::valid(const pair_generator::relative_pair& pair) const
{
  std::list<pair_filter_trait>::const_iterator trait;
  for(trait = my_traits.begin(); trait != my_traits.end(); ++trait)
  {
    if(!trait->valid(pair))
    {
      return false;
    }
  }
  return true;
}

//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep::iterator
//============================================================================
//
// - Find next pair that is valid for the filter.
//
filtering_pair_generator_rep::iterator*
filtering_pair_generator_rep::iterator::operator ++()
{
  while(my_iterator != *(my_rep->get_end()))
  {
    // Don't increment underlying iterator the first time thru.
    if(initializing)
    {
      initializing = false;
    }
    else
    { 
      ++my_iterator;
    }
    
    if(my_rep->get_filter()->valid(*my_iterator))
    {
      return this;
    }
  }
  return this;
}

//============================================================================
// IMPLEMENTATION:  filtering_pair_generator_rep::const_iterator
//============================================================================
//
// - Find next pair that is valid for the filter.
//
filtering_pair_generator_rep::const_iterator*
filtering_pair_generator_rep::const_iterator::operator ++()
{
  while(my_iterator != *(my_rep->get_end()))
  {
    // Don't increment underlying iterator the first time thru.
    if(initializing)
    {
      initializing = false;
    }
    else
    { 
      ++my_iterator;
    }
    
    if(my_rep->get_filter()->valid(*my_iterator))
    {
      return this;
    }
  }
  return this;
}

} // End namespace
