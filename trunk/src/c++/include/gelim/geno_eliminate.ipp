#include "gelim/geno_eliminate.h"

inline size_t
genotype_eliminator::build_family_data(imodel& model, const family_type& fam)
{
  // Verify that the family listed is valid (in the current subpedigree)

  if(fam.subpedigree() != subpedigree)
    return 1;

  // Find the parents

  mother_ptr = fam.parent1();
  father_ptr = fam.parent2();

  if(mother_ptr->is_male() || father_ptr->is_female())
    std::swap(mother_ptr, father_ptr);

  // Test the parents for validity.  If they're already empty, we return an
  // error state.

  if(model.unphased_penetrance_count(mother_ptr->subindex()+1) == 0 ||
     model.unphased_penetrance_count(father_ptr->subindex()+1) == 0 )
    return 2;

  // Find the first and last child

  family_type::offspring_const_iterator begin_child = fam.offspring_begin();
  family_type::offspring_const_iterator end_child   = fam.offspring_end();

  // Now define the children.

  my_child_info.resize(fam.offspring_count());
  
  bool valid_children = false;

  // Initialize children.

  size_t index = 0;

  for(family_type::offspring_const_iterator i  = begin_child; 
                                            i != end_child;
                                          ++i, ++index)
  {
    build_child_data(model, *i, my_child_info[index]);

    if(my_child_info[index].data.size() > 0)
      valid_children = true;
  }    

  // If there's no valid children, it means there's a problem.
  if(!valid_children) return 1;

  return 0;
}

inline void
genotype_eliminator::build_child_data
    (imodel& model, const individual_type& ind, child_info_type& child)
{
  child.ptr   = &ind;
  child.begin = model.phased_penetrance_begin(ind.subindex()+1);
  child.end   = model.phased_penetrance_end  (ind.subindex()+1);

  child.data.resize(model.phased_penetrance_count(ind.subindex()+1), false);

  child.data.clear();
}

inline void
genotype_eliminator::mark_uninformative(imodel& model, size_t marker)
{
  err.mark_info(subpedigree->pedigree(), marker, false, false);
}
