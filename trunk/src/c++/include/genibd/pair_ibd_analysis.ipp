// ======================================
// Inline functions of pair_ibd_analysis
// ======================================
  
inline
cerrorstream
pair_ibd_analysis::get_errors() const
{
  return errors;
}

inline
cerrorstream
pair_ibd_analysis::set_errors(SAGE::cerrorstream& s)
{
  errors = s;
  return errors;
}

inline bool
pair_ibd_analysis::built() const
{
  return my_built;
}

inline bool
pair_ibd_analysis::valid() const
{
  return my_valid;
}

inline SAGE::IBD*
pair_ibd_analysis::ibd_adaptor() const
{
  return my_ibds;
}

inline
double
pair_ibd_analysis::cond_tran (const phased_genotype& g, const conditional_genotype&  cg) const
{
  return 0.25 * ( (cg.b[0] && cg.child_geno[0] == g) +
                  (cg.b[1] && cg.child_geno[1] == g) +
                  (cg.b[2] && cg.child_geno[2] == g) +
                  (cg.b[3] && cg.child_geno[3] == g) );
}

inline
size_t
pair_ibd_analysis::has_allele(const phased_genotype& g, const allele& al) const
{
  return ( (al == g.allele1()) + (al == g.allele2()) );
}
