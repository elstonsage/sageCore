inline pedigree_imodel_generator::pedigree_imodel_generator()
 : prior_remap(true), post_remap(true), geno_elim(true)
{ }

inline bool pedigree_imodel_generator::do_prior_remap() const 
{ return prior_remap; }

inline bool pedigree_imodel_generator::do_post_remap()  const 
{ return post_remap;  }

inline bool pedigree_imodel_generator::do_genotype_elimination() const 
{ return geno_elim; }

inline void pedigree_imodel_generator::set_prior_remap(bool b)
{
  prior_remap = b;
}

inline void pedigree_imodel_generator::set_post_remap(bool b)
{
  post_remap = b;
}

inline void pedigree_imodel_generator::set_genotype_elimination (bool b)
{
  geno_elim = b; 
}

inline bool pedigree_imodel_generator::inconsistent() const
{
  return last_model_incon;
}

inline bool pedigree_imodel_generator::informative() const
{
  return last_model_inform;
}
