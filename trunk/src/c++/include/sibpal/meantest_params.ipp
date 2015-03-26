////////////////////////////////////////////////////////////////////////////
//             Implementation of meantest_params.h (Inline)               //
////////////////////////////////////////////////////////////////////////////

inline
void mean_estimate::set_pi(double p)
{
  my_pi = p;
}

inline
void mean_estimate::set_f1(double f)
{
  my_f1 = f;
}

inline
void mean_estimate::set_cov(const matrix& c)
{
  my_ss = c;
}

inline
void mean_estimate::set_w(double w)
{
  my_w = w;
}

inline
double mean_estimate::pi() const
{
  return my_pi;
}

inline
double mean_estimate::f1() const
{
  return my_f1;
}

inline
const matrix& mean_estimate::cov() const
{
  return my_ss;
}

inline
double mean_estimate::w() const
{
  return my_w;
}

////////////////////////////////////////////////////////////////////////////

inline
string marker_parameter::name(const relative_pairs& pairs) const
{
  return pairs.marker_name(marker);
}

inline
bool marker_parameter::operator==(const marker_parameter& p) const
{
  return marker == p.marker;
}

inline
bool marker_parameter::operator!=(const marker_parameter& p) const
{
  return !((*this) == p);
}

inline
bool marker_parameter::operator< (const marker_parameter& p) const
{
  if( marker < p.marker )
    return true;

  return false;
}

//////////////////////////////////////////////////////////////////////////

inline
string trait_parameter::name(const relative_pairs& pairs) const
{
  return pairs.trait_name(trait);
}

inline
size_t trait_parameter::affected_count() const
{
  if( affected_types() == 1 )
    for(int i=0; i < 3; ++i)
      if(affection[i])
        return i;
  return (size_t)-1;
}

inline
size_t trait_parameter::affected_types() const
{
  size_t a = 0;

  for(int i=0; i < 3; ++i)
    if(affection[i])
      ++a;
  return a;
}

inline
bool trait_parameter::is_set() const
{
  return (trait != (size_t)-1);
}

inline
bool trait_parameter::operator==(const trait_parameter& p) const
{
  return !((*this) != p);
}

inline
bool trait_parameter::operator!=(const trait_parameter& p) const
{
  if( trait != p.trait )
    return true;

  for(int i=0; i < 3; ++i)
    if(affection[i] != p.affection[i])
      return true;
  return false;
}

inline
bool trait_parameter::operator< (const trait_parameter& p) const
{
  if( trait < p.trait )
    return true;
  return false;
}

//////////////////////////////////////////////////////////////////////////

inline
const marker_parameter&
meantest_parameters::beta(size_t i) const
{
  return my_betas[i];
}

inline
const marker_parameter&
meantest_parameters::operator[](size_t i) const
{
  return beta(i);
}

inline
const trait_parameter&
meantest_parameters::subsets(size_t i) const
{
  return my_subsets[i];
}

inline
void
meantest_parameters::set_affection_status(int a)
{
  my_trait.set_affection(a);
}

inline
void
meantest_parameters::set_trait(size_t t)
{
  my_trait = trait_parameter(t);
}

inline
void
meantest_parameters::set_trait(size_t t, int a)
{
  my_trait = trait_parameter(t, a);
}

inline
void
meantest_parameters::validate()
{
  my_valid_regression = true;
}

inline
void
meantest_parameters::invalidate()
{
  my_valid_regression = false;
}

inline
bool
meantest_parameters::valid() const
{
  return my_valid_regression;
}

inline
int
meantest_parameters::affection_status() const
{
  return my_trait.affected_count();
}

inline
const trait_parameter&
meantest_parameters::trait() const
{
  return my_trait;
}

inline
bool
meantest_parameters::trait_set() const
{
  return my_trait.is_set();
}

inline
size_t
meantest_parameters::marker_count() const
{
  return my_betas.size();
}

inline
size_t
meantest_parameters::subset_count() const
{
  return my_subsets.size();
}

inline
void
meantest_parameters::set_pvalues_scientific_notation(bool a)
{
  my_pval_s_notation = a;
}

inline
bool
meantest_parameters::get_pvalues_scientific_notation() const
{
  return my_pval_s_notation;
}

inline
void
meantest_parameters::set_export_output(bool a)
{
  my_export_output = a;
}

inline
bool
meantest_parameters::get_export_output() const
{
  return my_export_output;
}

inline
marker_parameter&
meantest_parameters::beta(size_t i)
{
  return my_betas[i];
}

inline
marker_parameter&
meantest_parameters::operator[](size_t i)
{
  return beta(i);
}

inline
trait_parameter&
meantest_parameters::subsets(size_t i)
{
  return my_subsets[i];
}

inline
void
meantest_parameters::set_use_full_sibs(bool s)
{
  my_use_full_sibs = s;
}

inline
void
meantest_parameters::set_use_half_sibs(bool s)
{
  my_use_half_sibs = s;
}

inline
bool
meantest_parameters::get_use_full_sibs() const
{
  return my_use_full_sibs;
}

inline
bool
meantest_parameters::get_use_half_sibs() const
{
  return my_use_half_sibs;
}

inline
void
meantest_parameters::set_w(double w)
{
  my_w = w;
}

inline
double
meantest_parameters::get_w() const
{
  return my_w;
}
