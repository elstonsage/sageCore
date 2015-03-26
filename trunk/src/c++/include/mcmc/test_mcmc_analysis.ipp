//---------------------------------------------------------------------------
// Inline Implementation of test_mcmc_analysis
//---------------------------------------------------------------------------

inline const test_mcmc_parameters*
test_mcmc_analysis::get_parameters() const
{
  return my_parameters;
}

inline const RPED::RefMultiPedigree*
test_mcmc_analysis::get_multipedigree() const
{
  return my_multipedigree;
}
