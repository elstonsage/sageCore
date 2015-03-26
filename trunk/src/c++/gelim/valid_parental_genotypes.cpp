#include "gelim/valid_parental_genotypes.h"
#include "containers/bitfield.h"

namespace SAGE
{

struct unphased_parental_data
{ 
  //lint --e{1023}
  //lint --e{1703}
  //lint --e{1712}
  unphased_parental_data(const MLOCUS::inheritance_model& model, size_t id)
   : index(id),
     begin(model.unphased_penetrance_begin(id+1)),
     end  (model.unphased_penetrance_end  (id+1)),
     bits (model.unphased_penetrance_count(id+1), false)
  { }

  size_t index;

  MLOCUS::inheritance_model::unphased_penetrance_iterator begin;
  MLOCUS::inheritance_model::unphased_penetrance_iterator end;

  bit_field bits;
};

struct phased_parental_data
{ 
  //lint --e{1023}
  //lint --e{1703}
  //lint --e{1712}
  phased_parental_data(const MLOCUS::inheritance_model& model, size_t id)
   : index(id),
     begin(model.phased_penetrance_begin(id+1)),
     end  (model.phased_penetrance_end  (id+1)),
     bits (model.phased_penetrance_count(id+1), false)
  { }

  size_t index;

  MLOCUS::inheritance_model::phased_penetrance_iterator begin;
  MLOCUS::inheritance_model::phased_penetrance_iterator end;

  bit_field bits;
};

bool valid_parental_genotypes::generate_valid_parental_genotypes 
  (const family_type& fam, const imodel& model, bool phased)
{
  if(phased) return generate_valid_phased_parental_genotypes(fam,model);
  else       return generate_valid_unphased_parental_genotypes(fam,model);
}

bool valid_parental_genotypes::generate_valid_unphased_parental_genotypes 
  (const family_type& fam, const imodel& model)
{
  typedef imodel::unphased_penetrance_iterator unphased_iterator;
  typedef imodel::phased_penetrance_iterator   phased_iterator;

  // Clean out the parental genotypes

  clean_internals();

  // Construct mother and father data

  FPED::MemberConstPointer mother_id = fam.get_mother();
  FPED::MemberConstPointer father_id = fam.get_father();

  if(!mother_id)
  {
    mother_id = fam.parent1();
    father_id = fam.parent2();
  }

  unphased_parental_data mother(model, mother_id->subindex());
  unphased_parental_data father(model, father_id->subindex());

  // Find the indices of the sibship

  family_type::offspring_const_iterator begin_child = fam.offspring_begin();
  family_type::offspring_const_iterator end_child   = fam.offspring_end();

  // Create a list to store the results

  list<parental_genotype_pair> genotype_pairs;

  // For each parental pair, look for child genotypes that are valid.
  // If there is at least one valid genotype for each child, the parental
  // pair is valid

  unphased_iterator   gm = mother.begin;
  bit_field::iterator bm = mother.bits.begin();

  for( ; gm != mother.end; ++bm, ++gm)
  {
    unphased_iterator   gf = father.begin;
    bit_field::iterator bf = father.bits.begin();

    for( ; gf != father.end; ++bf, ++gf)
    {
      // Assume the set is valid until proven otherwise
      bool valid_set = true;

      // Set child genotypes - We use the reduced set since we don't care
      // about grand parental origins.  Similarly, we can reduce the set
      MLOCUS::child_genotype_set cg(gm.unphased_geno(), gf.unphased_geno(), true);
     
      // Iterate over children.  If ever a child sets valid_set to false,
      // the set isn't any good.
      family_type::offspring_const_iterator current_child  = begin_child; 
      for( ; valid_set && current_child != end_child; ++current_child)
      {
        size_t child_index = current_child->subindex();

        phased_iterator gc     = model.phased_penetrance_begin(child_index+1);
        phased_iterator gc_end = model.phased_penetrance_end  (child_index+1);

        // We only evaluate children that have genotypes.  Children that do
        // not would make everyone else invalid by contagion.
        if(gc != gc_end)
        {
          // The child is only valid if shown to be by this loop, so false
          // initially
          bool child_valid = false;

          // We don't have to check all the genotypes.  We only have to find
          // one that sets child_valid to true.
          for( ; !child_valid && gc != gc_end; ++gc)
          {
            if(cg.contains(gc.phased_geno()))
            {
              child_valid = true;
            }
          }

          // If no version is valid for the child, then this parental set
          // is not consistent.
          if(!child_valid) valid_set = false;
        }
      }

      // If all the children are valid, we validate the parental pair.
      if(valid_set)
      {
        genotype_pairs.push_back(make_pair(gm.geno_id(), gf.geno_id()));

        *bm = *bf = 1;
      }
    }
  }

  // Copy our list of valid elements into our vector

  my_genotype_pairs.resize(genotype_pairs.size());

  copy(genotype_pairs.begin(), genotype_pairs.end(), my_genotype_pairs.begin());

  // Construct the valid/invalid vectors for the mother

  gm = mother.begin;
  bm = mother.bits.begin();

  for( ; gm != mother.end; ++bm, ++gm)
  {
    if(*bm) my_mother_valid.push_back(gm.geno_id());
    else    my_mother_invalid.push_back(gm.geno_id());
  }

  // Construct the valid/invalid vectors for the father

  unphased_iterator   gf = father.begin;
  bit_field::iterator bf = father.bits.begin();

  for( ; gf != father.end; ++bf, ++gf)
  {
    if(*bf) my_father_valid.push_back(gf.geno_id());
    else    my_father_invalid.push_back(gf.geno_id());
  }

  return my_genotype_pairs.size() > 0;
}

bool valid_parental_genotypes::generate_valid_phased_parental_genotypes 
  (const family_type& fam, const imodel& model)
{
  typedef imodel::phased_penetrance_iterator   phased_iterator;

  // Clean out the parental genotypes

  clean_internals();

  // Construct mother and father data
  FPED::MemberConstPointer mother_id = fam.get_mother();
  FPED::MemberConstPointer father_id = fam.get_father();

  // If parents are unsexed, it doesn't really matter.
  if(!mother_id)
  {
    mother_id = fam.parent1();
    father_id = fam.parent2();
  }

  phased_parental_data mother(model, mother_id->subindex());
  phased_parental_data father(model, father_id->subindex());

  // Find the indices of the sibship

  family_type::offspring_const_iterator begin_child = fam.offspring_begin();
  family_type::offspring_const_iterator end_child   = fam.offspring_end();

  // For each parental pair, look for child genotypes that are valid.
  // If there is at least one valid genotype for each child, the parental
  // pair is valid

  phased_iterator     gm = mother.begin;
  bit_field::iterator bm = mother.bits.begin();

  for( ; gm != mother.end; ++bm, ++gm)
  {
    phased_iterator     gf = father.begin;
    bit_field::iterator bf = father.bits.begin();

    for( ; gf != father.end; ++bf, ++gf)
    {
      // Assume the set is valid until proven otherwise
      bool valid_set = true;

      // Set child genotypes - We use the reduced set since we don't care
      // about grand parental origins.
      MLOCUS::child_genotype_set cg(gm.phased_geno(), gf.phased_geno(), true);
     
      // Iterate over children.  If ever a child sets valid_set to false,
      // the set isn't any good.
      family_type::offspring_const_iterator current_child  = begin_child; 
      for( ; valid_set && current_child != end_child; ++current_child)
      {
        size_t child_index = current_child->subindex();

        phased_iterator gc     = model.phased_penetrance_begin(child_index+1);
        phased_iterator gc_end = model.phased_penetrance_end  (child_index+1);

        // We only evaluate children that have genotypes.  Children that do
        // not would make everyone else invalid by contagion.
        if(gc != gc_end)
        {
          // The child is only valid if shown to be by this loop, so false
          // initially
          bool child_valid = false;

          // We don't have to check all the genotypes.  We only have to find
          // one that sets child_valid to true.
          for( ; !child_valid && gc != gc_end; ++gc)
          {
            if(cg.contains(gc.phased_geno()))
            {
              child_valid = true;
            }
          }

          // If no version is valid for the child, then this parental set
          // is not consistent.
          if(!child_valid) valid_set = false;
        }
      }

      // If all the children are valid, we validate the parental pair.
      if(valid_set)
      {
        my_genotype_pairs.push_back(make_pair(gm.geno_id(), gf.geno_id()));

        *bm = *bf = 1;
      }
    }
  }

  // Construct the valid/invalid vectors for the mother

  gm = mother.begin;
  bm = mother.bits.begin();

  for( ; gm != mother.end; ++bm, ++gm)
  {
    if(*bm) my_mother_valid.push_back(gm.geno_id());
    else    my_mother_invalid.push_back(gm.geno_id());
  }

  // Construct the valid/invalid vectors for the father

  phased_iterator     gf = father.begin;
  bit_field::iterator bf = father.bits.begin();

  for( ; gf != father.end; ++bf, ++gf)
  {
    if(*bf) my_father_valid.push_back(gf.geno_id());
    else    my_father_invalid.push_back(gf.geno_id());
  }

  return my_genotype_pairs.size() > 0;
}


} // end namespace SAGE
