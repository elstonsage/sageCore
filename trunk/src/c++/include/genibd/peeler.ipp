//============================================================================
// INLINE IMPLEMENTATION:  peeler
//============================================================================
//
inline const peeler::subped_type&
peeler::subpedigree() const
{
  return my_subpedigree;
}

inline const inheritance_model&
peeler::marker() const
{
  return my_model;
}

// - Sum likelihood of and individual and his posterior over his unphased_genotype.
//
inline log_double
peeler::sum_ind_and_posterior(const ind_genotype& mg, const ind_genotype& fg, const member_type& ind)
{
  log_double  r(0);

  size_t ind_index = ind.subindex() + 1;

  phased_pen_iter iter = my_model.phased_penetrance_begin(ind_index);
  
  for( ; iter != my_model.phased_penetrance_end(ind_index); ++iter )
  {
    if(    my_model.is_x_linked()
        && ind.is_male()
        && !is_Y_genotype(iter.phased_geno()) )
      continue;

    MLOCUS::penetrance_model::phased_penetrance_iterator  ipi(iter);
    ind_genotype  ind_geno(ipi);

    double  trans = my_tcalc.transition(mg, fg, ind_geno);
    if( trans == 0 )
      continue;

    r += (*ipi) * posterior(ind, ipi) * trans;
  }
  
  return r;
}
