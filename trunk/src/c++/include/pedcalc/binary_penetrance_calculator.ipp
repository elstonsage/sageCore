#ifndef BINARY_PENETRANCE_CALCULATOR_H
#include "pedcalc/binary_penetrance_calculator.h"
#endif
namespace SAGE
{
namespace PED_CALC
{

/// Basic constructor.
///
/// \param susc_func The susceptibility functor.  Must have interface compatible with
///                  double (MEMBER_TYPE, GENOTYPE)
/// \param aff_func  The affection functor.  Must have interface compatible with
///                  bool (MEMBER_TYPE)
template <typename MEMBER_TYPE, typename GENOTYPE>
template <typename SUSC, typename AFF>
BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>::BinaryPenetranceCalculator
    (const SUSC& susc_func,
     const AFF&  aff_func)
  : my_susc_func (susc_func),
    my_aff_func  (aff_func)
{ }

/// Copy Constructor
///
/// \param s The object to copy
template <typename MEMBER_TYPE, typename GENOTYPE>
BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>::BinaryPenetranceCalculator
    (const SelfType& s)
  : my_susc_func (s.my_susc_func),
    my_aff_func  (s.my_aff_func)
{ }
    
/// Destructor
///
template <typename MEMBER_TYPE, typename GENOTYPE>
BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>::~BinaryPenetranceCalculator()
{ }
    
/// Copy Operator
///
/// \param s The object to copy
template <typename MEMBER_TYPE, typename GENOTYPE>
BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>&
  BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>::operator=
    (const SelfType& s)
{
  if(this != &s)
  {
    my_susc_func = s.my_susc_func;
    
    my_aff_func  = s.my_aff_func;
  }
  
  return *this;
}
    
/// The calculation operator.  Calculates:
///
/// \f[
///    Pen(m, u) = \frac{e^{theta_u(m)}y_m}{1+e^{theta_u(m)}}
/// \f]
///
/// where \f$theta_u(m)\f$ is the susceptibility and 
/// \f$y_m\f$ is the affection status for member \i m.
///
/// \param m The member
/// \param g The genotype
template <typename MEMBER_TYPE, typename GENOTYPE>
double BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE>::operator()
    (const MEMBER_TYPE& m,
     const GENOTYPE&    g) const
{
  double susceptibility = my_susc_func (m, g);
  bool   affection      = my_aff_func  (m);
  
  return (exp(susceptibility * affection)) / (1.0 + exp(susceptibility));
}

}
}

