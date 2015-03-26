// ===========================================
// Inline Implementation of sim_mcmc_simulator
// ===========================================

inline
vector<sim_relative_pair>*
ibd_mcmc_simulator::set_pairs(vector<sim_relative_pair>* rp)
{
  return my_relative_pairs = rp;
}

inline
vector<sim_relative_pair>*
ibd_mcmc_simulator::get_pairs() const
{
  return my_relative_pairs;
}

inline
void
ibd_mcmc_simulator::singlepoint_ibd_sharing(int offset)
{
  ibd_sharing(offset, my_current_marker);
}

inline
void
ibd_mcmc_simulator::multipoint_ibd_sharing(int offset)
{
  for( size_t j=0; j < my_total_loci; j++ )
    if( my_markers_used[j] )
      ibd_sharing(offset, j);
}
