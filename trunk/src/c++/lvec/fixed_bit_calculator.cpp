//==========================================================================
//  File:    fixed_bit_calculator.cpp
//
//  Author:  Geoff Wedig
//
//  History: 0.1 Initial Implementation
//           1.0 Updated to new libraries                        yjs Sep. 04
//
//  Notes:
//
//  Copyright (c) 1998 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "lvec/fixed_bit_calculator.h"

namespace SAGE
{

ostream& operator<<(ostream& o, const vector<bool>& v)
{
  vector<bool>::const_iterator i = v.begin();

  for( ; i != v.end(); ++i)
  {
    if(*i) o << 'X';
    else   o << '.';
  }

  return o;
}

// ====================
// fixed_bit_container
// ====================

void
fixed_bit_container::dump(ostream& o) const
{
  double bits_saved = 0;
  
  size_t bit_count = 0;
  
  OUTPUT::Section info;

  // Output member info
  
  OUTPUT::Table   member_info;
  
  member_info << OUTPUT::TableColumn("Name")
              << OUTPUT::TableColumn("M. Fixed")
              << OUTPUT::TableColumn("F. Fixed")
              << OUTPUT::TableColumn("M. Don't Care")
              << OUTPUT::TableColumn("F. Don't Care");
  
  for(size_t i = 0; i < my_ped->member_count(); ++i)
  {
    const FPED::Member& mem = my_ped->member_index(i);
    
    if(!mem.is_founder())
    {
      bit_count += 2;
      
      OUTPUT::TableRow trow;
      
      trow << mem.name()
           << mother_fixed(mem)
           << father_fixed(mem)
           << mother_dont_care(mem)
           << father_dont_care(mem);

      member_info << trow;
    }
  }
  
  info << member_info;
  
  // Output Family Info
  OUTPUT::Section  family_info;
  
  for(FPED::FamilyConstIterator fam = my_ped->family_begin();
      fam != my_ped->family_end(); ++fam)
  {
    OUTPUT::Table family_table(fam->parent1()->name() + " x " + fam->parent2()->name());
    
    family_table << OUTPUT::TableColumn("");
    
    // Create names for header
    for(FPED::OffspringConstIterator off1 = fam->offspring_begin();
        off1 != fam->offspring_end(); ++off1)
      family_table << OUTPUT::TableColumn(off1->name());

    family_table << OUTPUT::TableColumn("|");
    
    for(FPED::OffspringConstIterator off1 = fam->offspring_begin();
        off1 != fam->offspring_end(); ++off1)
      family_table << OUTPUT::TableColumn(off1->name());

    // Create row for each child
    for(FPED::OffspringConstIterator off1 = fam->offspring_begin();
        off1 != fam->offspring_end(); ++off1)
    {
      double m_sync = 0;
      double f_sync = 0;
      
      OUTPUT::TableRow row;
      
      row << off1->name();

      for(FPED::OffspringConstIterator off2 = fam->offspring_begin();
          off2 != fam->offspring_end(); ++off2)
      {
        row << mother_synchronized(*off1, *off2);

        if(mother_synchronized(*off1, *off2)) ++m_sync;
      }
      row << "|";

      for(FPED::OffspringConstIterator off2 = fam->offspring_begin();
          off2 != fam->offspring_end(); ++off2)
      {
        row << father_synchronized(*off1, *off2);

        if(father_synchronized(*off1, *off2)) ++f_sync;
      }
      
      family_table << row;

      if(mother_fixed(*off1) || mother_dont_care(*off1)) ++bits_saved;
      else                                               bits_saved += (m_sync) / (m_sync+1);

      if(father_fixed(*off1) || father_dont_care(*off1)) ++bits_saved;
      else                                               bits_saved += (f_sync) / (f_sync+1);

    }
    
    // If the family is big enough, add the output.
    if(fam->offspring_count() > 1)
      family_info << family_table;
  }
  
  info << family_info;
  
  // Output Info
  
  o << info;
  
  o << "Bits saved: " << bits_saved << '/' << bit_count << endl;
}

// ====================
// fixed_bit_calculator
// ====================

fixed_bit_calculator::fixed_bit_calculator
    (const FPED::Subpedigree& p,
     const MLOCUS::inheritance_model& imodel)
  : my_data(p)
{
  for(size_t i = 0; i < p.family_count(); ++i)
  {
    const FPED::Family& fam = p.family_index(i);

    run_nuclear_family(fam, imodel);
  }
}

void
fixed_bit_calculator::dump(ostream& o) const
{
  my_data.dump(o);
}

void
fixed_bit_calculator::run_nuclear_family
    (const FPED::Family& fam,
     const MLOCUS::inheritance_model& imodel)
{
#if 0
  cout << "fixed_bit_calculator::run_nuclear_family(fam = "
       << fam.name() << ")..." << endl;

  cout << "offspring count = " << fam.offspring_count() << endl;
  cout << "parents : " << endl;
  cout << " " << fam.parent1()->name() << endl
       << " " << fam.parent2()->name() << endl;

  FPED::OffspringConstIterator oi = fam.offspring_begin();
  cout << "offsprings : " << endl;
  for( ; oi != fam.offspring_end(); ++oi )
    cout << " " << oi->name();
  cout << endl;
#endif

  my_current_family = &fam;

  test_parent(*fam.parent1(), imodel);
  test_parent(*fam.parent2(), imodel);

  if(my_potential_fixed_bit.size() || my_potential_synchronized_bits.size())
  {
    test_family(imodel);

    potential_iterator j = my_potential_fixed_bit.begin();
    for( ; j != my_potential_fixed_bit.end(); ++j)
      my_data.set_parent_bit_type(j->is_father, *j->id1, fixed_bit_container::fixed);
      
    j = my_potential_synchronized_bits.begin();
    for( ; j != my_potential_synchronized_bits.end(); ++j)
      my_data.synchronize_parent(j->is_father, *j->id1, *j->id2);
  }

  my_potential_fixed_bit         = potential_list();
  my_potential_synchronized_bits = potential_list();
}

void
fixed_bit_calculator::test_parent
    (const FPED::Member& par,
     const MLOCUS::inheritance_model& imodel)
{
  // Check parent for easy things

  vector<int> allele_sources(imodel.allele_count(), 0);

  size_t parent = par.subindex();
  
  // The type of the parent determined how much we can predetermine before
  // expensive visiting of all possibilities.
  //   0 - all children no_rel bits (parent has both homo and hetero genotypes)
  //   1 - children may be fixed.  Any given parental allele has only one gparental source
  //   2 - children are all don't care (parent has all genotypes homozygous)
  //   3 - No fixed bits.  At least one
  //       parental allele can come from either grandparent.

  // Start with 1, then disprove it
  int type = 1;

  // Keep track of finding homo and heterozygous genotypes
  bool found_heterozygous = false;
  bool found_homozygous   = false;

  MLOCUS::penetrance_model::phased_penetrance_iterator g =
      imodel.phased_penetrance_begin(parent+1);

  for( ; g != imodel.phased_penetrance_end(parent+1); ++g)
  {
    if(g.phased_geno().homozygous())
    {
      found_homozygous = true;

      if(found_heterozygous)
      {
        // Found both, there's no relationships that we can detect
        
        type = 0; // All children no_rel
        break;
      } 
    }
    else
    {
      found_heterozygous = true;

      if(found_homozygous)
      {
        // Found both, there's no relationships that we can detect
        
        type = 0; // All children no_rel
        break;
      }

      // Store sources of the particular alleles of the parentally phased genotype
      allele_sources[g.phased_geno().allele1().id()] |= 1;
      allele_sources[g.phased_geno().allele2().id()] |= 2;
      
      if(allele_sources[g.phased_geno().allele1().id()] > 2 ||
         allele_sources[g.phased_geno().allele2().id()] > 2)
      {
        type = 3;
        // We don't discontinue testing at this point, because there may be
        // synchronization, though there are no fixed bits.  IF we find a 
        // homozygous genotype, though, we're no_rel all the way.
      }
    }
  }
  
  // If we never found heterozygous genotypes, we're type 2
  if(!found_heterozygous)
    type = 2;
  
  bool is_father = !par.is_female();
  
  // Depending on type, do the appropriate thing with the children.
  switch(type)
  {
    case 0: set_all_no_rel(is_father);                               break;

    case 1: do_basic_fixed_check(is_father, imodel, allele_sources); break;

    case 2: set_all_homozygous(is_father);                           break;

    case 3: create_all_potential_sync(is_father);                    break;
  }
}

void
fixed_bit_calculator::set_all_no_rel(bool is_father)
{
  // Set all children to no_rel for parent p
  for(FPED::OffspringConstIterator i = my_current_family->offspring_begin();
      i != my_current_family->offspring_end(); ++i)
  {
    my_data.set_parent_bit_type(is_father, *i, fixed_bit_container::no_rel);
  }

  create_all_potential_sync(is_father);
}

// Check basic fixed status with regards to a single parent.  This is
// not exhaustive, sinced fixed status may be determined by both parents
// states, but this does the easy level first, and stores away the
// more difficult stuff for exhaustive testing.
void
fixed_bit_calculator::do_basic_fixed_check
  (bool                             is_father,
   const MLOCUS::inheritance_model& imodel,
   const vector<int>&               allele_sources)
{
  // Look for fixed bits in children of the particular parent.
  
  // Create a vector to store the children's fixed status.  Initially, we assume
  // they are fixed and disprove it.
  vector<bool> fixed_children(my_current_family->offspring_count(), true);

  for(FPED::OffspringConstIterator i = my_current_family->offspring_begin();
      i != my_current_family->offspring_end(); ++i)
  {
    int allele_source = 0;
    
    size_t sibling_index = my_data.my_ind_sib_indices[i->subindex()];

    // Iterate over phased genotypes of the i'th child
    MLOCUS::penetrance_model::phased_penetrance_iterator g =
        imodel.phased_penetrance_begin(i->subindex()+1);

    for( ; fixed_children[sibling_index] == true &&
           g != imodel.phased_penetrance_end(i->subindex()+1); ++g)
    {
      const MLOCUS::phased_genotype& gt = g.phased_geno();
      
      // Get the allele ids, sorted so the first allele comes from the 
      // parent in question
      int al;

      if(!is_father) // mother
      {
        al = gt.allele1().id();
      }
      else       // father
      {
        al = gt.allele2().id();
      }

      allele_source |= allele_sources[al];
      
      if(allele_source > 2)
        fixed_children[sibling_index] = false;
    }
  }

  // Determine fixed and potential sync bits based upon fixed information we've
  // collected.
  for(FPED::OffspringConstIterator i = my_current_family->offspring_begin();
      i != my_current_family->offspring_end(); ++i)
  {
    size_t sibling_index = my_data.my_ind_sib_indices[i->subindex()];
    if(fixed_children[sibling_index])
    {
      // We know this child is fixed, so label them.
      my_data.set_parent_bit_type(is_father, *i, fixed_bit_container::fixed);
      
      FPED::OffspringConstIterator j = i;
      ++j;
      for( ; j != my_current_family->offspring_end(); ++j)
      {
        size_t sibling_index1 = my_data.my_ind_sib_indices[j->subindex()];

        if(fixed_children[sibling_index1])
        {
          // If both are fixed, then they're synced by definition
          my_data.synchronize_parent(is_father, *i, *j);
        }
        else
        {
          // Otherwise, they could still be synced, but we need both parents
          // to figure it out.
          my_potential_synchronized_bits.push_back(potential_data(is_father, &*i, &*j));
        }
      }
    }
    else
    {
      // Bit is still potentially fixed, but it will require both parents
      // to determine, so store for later.  
      my_potential_fixed_bit.push_back(potential_data(is_father, &*i));
      
      // Label all other children as potentially synced with this one.
      FPED::OffspringConstIterator j = i;
      ++j;
      for( ; j != my_current_family->offspring_end(); ++j)
      {
        my_potential_synchronized_bits.push_back(potential_data(is_father, &*i, &*j));
      }
    }
  }
}

void
fixed_bit_calculator::set_all_homozygous(bool is_father)
{
  for(FPED::OffspringConstIterator i = my_current_family->offspring_begin();
      i != my_current_family->offspring_end(); ++i)
  {
    my_data.set_parent_bit_type(is_father, *i, fixed_bit_container::dont_care);
  }
}

void
fixed_bit_calculator::create_all_potential_sync(bool is_father)
{
  for(FPED::OffspringConstIterator i = my_current_family->offspring_begin();
      i != my_current_family->offspring_end(); ++i)
  {
    FPED::OffspringConstIterator j = i;
    ++j;
    for( ; j != my_current_family->offspring_end(); ++j)
    {
      my_potential_synchronized_bits.push_back(potential_data(is_father, &*i, &*j));
    }
  }
}

void
fixed_bit_calculator::test_family(const MLOCUS::inheritance_model& imodel)
{
  FPED::MemberConstPointer mother = my_current_family->parent1();
  FPED::MemberConstPointer father = my_current_family->parent2();
  
  // Make sure mother and father are sorted
  if( mother->is_male() )
    std::swap(mother, father);
  
  size_t moth_id = mother->subindex();
  size_t fath_id = father->subindex();

  bool b = my_potential_fixed_bit.size() || my_potential_synchronized_bits.size();

  MLOCUS::penetrance_model::phased_penetrance_iterator m_genotype =
      imodel.phased_penetrance_begin(moth_id+1);
  for( ; b && m_genotype != imodel.phased_penetrance_end(moth_id+1); ++m_genotype)
  {
    MLOCUS::penetrance_model::phased_penetrance_iterator f_genotype =
        imodel.phased_penetrance_begin(fath_id+1);
    for( ; b && f_genotype != imodel.phased_penetrance_end(fath_id+1); ++f_genotype)
    {
      test_children(imodel, m_genotype.phased_geno(),
                            f_genotype.phased_geno());
      
      b = my_potential_fixed_bit.size() || my_potential_synchronized_bits.size();
    }
  }
}

// Tests the children against a particular parental genoset
void
fixed_bit_calculator::test_children(const MLOCUS::inheritance_model& imodel,
                                    const MLOCUS::phased_genotype& mg,
                                    const MLOCUS::phased_genotype& fg)
{
  // Build the possible child genotypes.
  MLOCUS::child_genotype_set child_genos(mg, fg);
  
  // Determine the possible iv bits for each child.

  bool valid = true;

  vector<pair<int, int> > c_flags(my_current_family->offspring_count(), make_pair(0,0));
  
  // For each child, determine what their possible bit patterns are given the
  // parental genotypes.  If ever a child is invalid, then we know the parents
  // cannot have those genotypes
  FPED::OffspringConstIterator oi = my_current_family->offspring_begin();
  for( size_t i = 0; valid & oi != my_current_family->offspring_end(); ++oi, ++i )
  {
    MLOCUS::penetrance_model::phased_penetrance_iterator cg =
        imodel.phased_penetrance_begin(oi->subindex()+1);
        
    pair<int, int>& cf = c_flags[i];
    
    // Determine which bit patterns fit for the child's possible genotypes
    for( ; cg != imodel.phased_penetrance_end(oi->subindex()+1); ++cg )
    {
      const MLOCUS::phased_genotype& cgg = cg.phased_geno();
      
      if( cgg == child_genos[0])
      {
        cf.first |= 1; cf.second |= 1;
      }
      if( cgg == child_genos[1])
      {
        cf.first |= 1; cf.second |= 2;
      }
      if( cgg == child_genos[2])
      {
        cf.first |= 2; cf.second |= 1;
      }
      if( cgg == child_genos[3])
      {
        cf.first |= 2; cf.second |= 2;
      }
    }
    
    if(!cf.first || !cf.second) valid = false;
  }

  if(!valid) return;
  
  // We know the parents genotypes result in valid children at this point
  
  // Test the fixed bits
  
  potential_iterator i = my_potential_fixed_bit.begin();
  
  for( ; i != my_potential_fixed_bit.end(); )
  {
    // Initially assume the bit may be fixed
    bool bit_fixed = true;

    size_t sib_index = my_data.my_ind_sib_indices[i->id1->subindex()];
    
    if(!i->is_father)  // Mother Check
    {
      // If the child could have had either allele from the parent, we're not
      // fixed
           if(c_flags[sib_index].first > 2)                 bit_fixed = false;

      // If we haven't identified any alleles from the parent yet store 
      // the pattern for later comparisons
      else if(i->bit_pattern == -1)                      i->bit_pattern = c_flags[sib_index].first;

      // If the bits from the parent that we've previously identified don't match
      // the present one, we're not fixed
      else if(i->bit_pattern != c_flags[sib_index].first)   bit_fixed = false;
    }
    else
    {
      // If the child could have had either allele from the parent, we're not
      // fixed
           if(c_flags[sib_index].second > 2)                bit_fixed = false;

      // If we haven't identified any alleles from the parent yet store 
      // the pattern for later comparisons
      else if(i->bit_pattern == -1)                      i->bit_pattern = c_flags[sib_index].second;

      // If the bits from the parent that we've previously identified don't match
      // the present one, we're not fixed
      else if(i->bit_pattern != c_flags[sib_index].second)  bit_fixed = false;
    }

    if(!bit_fixed)
    {
      // Erase this potential iterator
      potential_iterator j = i;
      ++i;
      
      my_potential_fixed_bit.erase(j);
    }
    else
      ++i;
  }

  // Test the synchronization bits

  i = my_potential_synchronized_bits.begin();
  
  for( ; i != my_potential_synchronized_bits.end(); )
  {
    // Initially assume they're synced.
    bool bits_synced = true;
    
    size_t sib_index1 = my_data.my_ind_sib_indices[i->id1->subindex()];
    size_t sib_index2 = my_data.my_ind_sib_indices[i->id1->subindex()];

    if(!i->is_father) // Test mother
    {
      // Get the valid bit patterns for the mother for both children
      const int& c1 = c_flags[sib_index1].first;
      const int& c2 = c_flags[sib_index2].first;
      
      // If either child could get either allele from the mother, they're
      // not synced
      if(c1 > 2 || c2 > 2)                 bits_synced = false;
      
      // Otherwise, if we've not seen any patterns before, store the pattern.
      // In this case, the pattern indicates whether the two children were the same
      // or different by xor.
      else if(i->bit_pattern == -1)        i->bit_pattern = c1 ^ c2;
      
      // If the present pattern doesn't match previous patterns, then they're
      // not synced
      else if(i->bit_pattern != (c1 ^ c2)) bits_synced = false;
    }
    else
    {
      // Get the valid bit patterns for the father for both children
      const int& c1 = c_flags[sib_index1].second;
      const int& c2 = c_flags[sib_index2].second;

      // If either child could get either allele from the mother, they're
      // not synced
      if(c1 > 2 || c2 > 2)                 bits_synced = false;
      
      // Otherwise, if we've not seen any patterns before, store the pattern.
      // In this case, the pattern indicates whether the two children were the same
      // or different by xor.
      else if(i->bit_pattern == -1)        i->bit_pattern = c1 ^ c2;
      
      // If the present pattern doesn't match previous patterns, then they're
      // not synced
      else if(i->bit_pattern != (c1 ^ c2)) bits_synced = false;
    }

    if(!bits_synced)
    {
      potential_iterator j = i;
      ++i;
      
      my_potential_synchronized_bits.erase(j);
    }
    else
      ++i;
  }
}

}

