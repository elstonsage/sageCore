//============================================================================
// File:      phenoset.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/28/2 - created.                                djb
//                                                                          
// Notes:     inline implementation for class, phenoset.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  phenoset
//============================================================================
//
inline
phenoset::phenoset(const subped_type& sp, size_t trait, size_t marker,
                   const member_type& ind)
      : my_trait(trait), my_marker(marker), my_ind(ind),
        my_trait_pm(ge_models::get_model(sp, trait)),
        my_marker_pm(ge_models::get_model(sp, marker)),
        my_trait_phenotype(ind.subindex() + 1),
        my_marker_phenotype(ind.subindex() + 1)
{}

inline phenoset::phenoset_iterator
phenoset::begin() const
{
  return phenoset_iterator(trait_iter_begin(), marker_iter_begin(), this);
}

inline phenoset::phenoset_iterator
phenoset::end() const
{
  return phenoset_iterator(trait_iter_end(), marker_iter_end(), this);
}

inline phenoset::phased_penetrance_iterator  
phenoset::trait_iter_begin() const
{
  return my_trait_pm.phased_penetrance_begin(my_trait_phenotype);
}

inline phenoset::phased_penetrance_iterator  
phenoset::marker_iter_begin() const
{
  return my_marker_pm.phased_penetrance_begin(my_marker_phenotype);
}

inline phenoset::phased_penetrance_iterator  
phenoset::trait_iter_end() const
{
  return my_trait_pm.phased_penetrance_end(my_trait_phenotype);
}

inline phenoset::phased_penetrance_iterator  
phenoset::marker_iter_end() const
{
  return my_marker_pm.phased_penetrance_end(my_marker_phenotype);
}


//============================================================================
// IMPLEMENTATION:  phenoset::phenoset_iterator
//============================================================================
//
inline
phenoset::phenoset_iterator::phenoset_iterator
          (const phased_penetrance_iterator& trait_iter,
           const phased_penetrance_iterator& marker_iter,
           const phenoset* ph_set)
      : my_trait_iter(trait_iter), my_marker_iter(marker_iter),
        my_phenoset(ph_set)
{}

inline joint_pen_iter
phenoset::phenoset_iterator::operator*() const
{
  return joint_pen_iter(my_trait_iter, my_marker_iter);
}

inline phenoset::phenoset_iterator&
phenoset::phenoset_iterator::operator++()
{
  if(my_marker_iter != my_phenoset->marker_iter_end())
  {
    ++my_marker_iter;
    if(my_marker_iter == my_phenoset->marker_iter_end())
    {
      ++my_trait_iter;
      if(my_trait_iter != my_phenoset->trait_iter_end())
      {
        my_marker_iter = my_phenoset->marker_iter_begin();
      }
    }
  }
    
  return *this;
}

inline bool
phenoset::phenoset_iterator::operator==(const phenoset_iterator& other) const
{
  return my_trait_iter  == other.my_trait_iter  &&
         my_marker_iter == other.my_marker_iter &&
         my_phenoset    == other.my_phenoset      ;
}

inline bool
phenoset::phenoset_iterator::operator!=(const phenoset_iterator& other) const
{
  return  ! operator==(other);
}
                   


