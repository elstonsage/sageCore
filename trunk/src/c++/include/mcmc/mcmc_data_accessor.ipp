//===========================================
// Inline Implemantation: mcmc_data_accessor 
//===========================================

inline
bool
mcmc_data_accessor::is_valid_locus(size_t m) const
{
  return my_valid_loci[m];
}

inline
bool 
mcmc_data_accessor::is_geno_miss(size_t id, size_t m)  const
{
  assert(id < individual_count());
  assert(m  < locus_count());

//  size_t ped_i = my_mcmc_meiosis_map[id]->index();

  const FPED::Member& mem = my_mcmc_meiosis_map.get_subpedigree().member_index(id);

  const MLOCUS::inheritance_model& model = my_mcmc_meiosis_map.get_multipedigree().info().marker_info(m);

  return mem.info().phenotype_missing(m, model);
}

inline
int
mcmc_data_accessor::mother_bit(size_t id, size_t m) const
{
  assert(id < individual_count());
  assert(m  < locus_count());

  const FPED::Member& mem = my_mcmc_meiosis_map.get_subpedigree().member_index(id);
  
  return my_indicator[m][my_mcmc_meiosis_map.get_mother_meiosis(mem)];
}

inline
int
mcmc_data_accessor::father_bit(size_t id, size_t m) const
{
  assert(id < individual_count());
  assert(m  < locus_count());

  const FPED::Member& mem = my_mcmc_meiosis_map.get_subpedigree().member_index(id);
  
  return my_indicator[m][my_mcmc_meiosis_map.get_father_meiosis(mem)];
}

inline
size_t
mcmc_data_accessor::individual_count() const
{
  return my_mcmc_meiosis_map.get_subpedigree().member_count();
}

inline
size_t
mcmc_data_accessor::locus_count() const
{
  return my_locus_count;
}

inline
void
mcmc_data_accessor::set_valid_locus(size_t m, bool b)
{
  my_valid_loci[m] = b;
}

inline
mcmc_data_accessor::indicator_type&
mcmc_data_accessor::get_indicator()
{
  return my_indicator;
}

inline
const mcmc_data_accessor::indicator_type&
mcmc_data_accessor::get_indicator() const
{
  return my_indicator;
}

inline
bit_field&
mcmc_data_accessor::get_indicator(size_t m)
{
  return my_indicator[m];
}

inline
const bit_field&
mcmc_data_accessor::get_indicator(size_t m) const
{
  return my_indicator[m];
}

inline
const McmcMeiosisMap&
mcmc_data_accessor::get_mcmc_meiosis_map() const
{
  return my_mcmc_meiosis_map;
}
