////////////////////////////////////////////////////////////////////////////
//             Implementation of meantest.h (Inline)                      //
////////////////////////////////////////////////////////////////////////////

inline
void
SibMeanTest::set_parameters(const meantest_parameters& p)
{
  my_parameters = p;
}

inline
meantest_parameters&
SibMeanTest::parameters()
{
  return my_parameters;
}

inline
const meantest_parameters&
SibMeanTest::parameters() const
{
  return my_parameters;
}

inline void
SibMeanTest::set_use_pairs(const pair<bool, bool>& up)
{
  my_parameters.set_use_full_sibs(up.first);
  my_parameters.set_use_half_sibs(up.second);
}

inline pair<bool, bool>
SibMeanTest::get_use_pairs() const
{
  return make_pair(my_parameters.get_use_full_sibs(), my_parameters.get_use_half_sibs());
}
