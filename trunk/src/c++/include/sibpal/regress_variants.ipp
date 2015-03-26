////////////////////////////////////////////////////////////////////////////
//             Implementation of regressvariants.h (Inline)               //
////////////////////////////////////////////////////////////////////////////

//
//------------------------------------------------------------------------
//
/*
inline
double
VariantImplementation::p0_correlation(const trait_parameter& trait, bool empirical) const
{
  double p0 = fixed_p0_correlation();

  if( empirical && !finite(p0) )
    p0 = trait.p_empirical_correlation[0];

  if( !finite(p0) )
    // per Dr. Elston, always start with identity matrix
    //p0 = trait.p_correlation[0];
    p0 = 0.0;

  return p0;
}

inline
double
VariantImplementation::p1_correlation(const trait_parameter& trait, bool empirical) const
{
  double p1 = fixed_p1_correlation();

  if( empirical && !finite(p1) )
    p1 = trait.p_empirical_correlation[1];

  if( !finite(p1) )
    // per Dr. Elston, always start with identity matrix
    //p1 = trait.p_correlation[1];
    p1 = 0.0;

  return p1;
}
*/
inline
string
VariantImplementation::description() const
{
  return "No description provided";
}

//
//------------------------------------------------------------------------
//

inline
string
SumSquared::description() const
{
  return "Mean-corrected squared trait sum";
}
/*
inline
double
SumSquared::p0_correlation(const trait_parameter& trait, bool use_empirical_correlations) const
{
  double p0 = fixed_p0_correlation();

  if( !finite(p0) && use_empirical_correlations )
      p0 = trait.p_empirical_correlation[0];

  if( !finite(p0) )
      p0 = trait.p_correlation[0];

  if( !finite(p0) )
    p0 = trait.p_sum_correlation[0];

  return p0;
}

inline
double
SumSquared::p1_correlation(const trait_parameter& trait, bool use_empirical_correlations) const
{
  double p1 = fixed_p1_correlation();

  if( !finite(p1) && use_empirical_correlations )
      p1 = trait.p_empirical_correlation[1];

  if( !finite(p1) )
      p1 = trait.p_correlation[1];

  if( !finite(p1) )
    p1 = trait.p_sum_correlation[1];

  return p1;
}
*/
//
//------------------------------------------------------------------------
//

inline
string
DifferenceSquared::description() const
{
  return "Squared trait difference";
}
/*
inline
double
DifferenceSquared::p0_correlation(const trait_parameter& trait, bool use_empirical_correlations) const
{
  double p0 = fixed_p0_correlation();

  if( !finite(p0) && use_empirical_correlations )
      p0 = trait.p_empirical_correlation[0];

  if( !finite(p0) )
      p0 = trait.p_correlation[0];

  if( !finite(p0) )
    p0 = trait.p_diff_correlation[0];

  return p0;
}

inline
double
DifferenceSquared::p1_correlation(const trait_parameter& trait, bool use_empirical_correlations) const
{
  double p1 = fixed_p1_correlation();

  if( !finite(p1) && use_empirical_correlations )
      p1 = trait.p_empirical_correlation[1];

  if( !finite(p1) )
      p1 = trait.p_correlation[1];

  if( !finite(p1) )
    p1 = trait.p_diff_correlation[1];

  return p1;
}
*/
//
//------------------------------------------------------------------------
//

inline
string
CrossProduct::description() const
{
  return "Mean-corrected trait cross-product";
}

//
//------------------------------------------------------------------------
//

inline
string
WeightedVariance::description() const
{
  return "Weighted square trait sum and difference (W2)";
}

//
//------------------------------------------------------------------------
//

inline
string
WeightedCorrelatedVariance::description() const
{
  return "Weighted square trait sum and difference (W3)";
}

inline
bool
WeightedCorrelatedVariance::apply_weighting() const
{
  return false;
}

//
//------------------------------------------------------------------------
//

inline
string
WeightedCorrelatedVarianceAndTraits::description() const
{
  return "Weighted square trait sum and difference (W4)";
}

