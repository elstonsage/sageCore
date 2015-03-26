//============================================================================
// IMPLEMENTATION:  trans_calculator
//============================================================================
//
inline  
trans_calculator::trans_calculator()
{}

// - Kronecker delta.
//
inline bool
trans_calculator::kd(const allele& one, const allele& two)
{
  return (one == two) ? true : false;
}

inline double
trans_calculator::transition(const ind_genotype& mom,
                             const ind_genotype& dad,
                             const ind_genotype& kid ) const
{
  double  trans1 = transmission(mom, MPED::SEX_FEMALE, kid.mg.allele1());
  double  trans2 = transmission(dad, MPED::SEX_MALE,   kid.mg.allele2());

  return trans1 * trans2;
}

inline double
trans_calculator::transmission(const ind_genotype& p, MPED::SexCode s, const allele& ka) const
{
  const allele& p_allele1 = p.mg.allele1();
  const allele& p_allele2 = p.mg.allele2();

  if( p_allele1 != ka && p_allele2 != ka )
    return 0.;

  double  term1 = kd(p_allele1, ka);
  double  term2 = kd(p_allele2, ka);

  return (term1 + term2) / 2.0;
}
