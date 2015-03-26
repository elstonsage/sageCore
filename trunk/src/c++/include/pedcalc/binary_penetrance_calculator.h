#ifndef BINARY_PENETRANCE_CALCULATOR_H
#define BINARY_PENETRANCE_CALCULATOR_H

#include "boost/function.hpp"
#include <functional>
#include <math.h>

namespace SAGE
{
namespace PED_CALC
{
  
/// \brief Calculates the penetrance of a member for a genotype given susceptibility and affection functors
///
/// The BinaryPenetranceCalculator calculates penetrances of a binary trait
/// given fuctors which return the susceptibility and the affection status.
///
/// The equation used is
///
/// \f[
///    Pen(m, u) = \frac{e^{\theta_u(m)}y_m}{1+e^{\theta_u(m)}}
/// \f]
///
/// where \f$\theta_u(m)\f$ is the susceptibility and 
/// \f$y_m\f$ is the affection status for member \i m.
///
/// The functors which are given must meet the signature (or be convertible to
/// the signature.  See boost function references):
///
///  - susceptibility: double (const MEMBER_TYPE&, const GENOTYPE&)
///  - affection: bool (const MEMBER_TYPE&);
template <typename MEMBER_TYPE, typename GENOTYPE>
class BinaryPenetranceCalculator : public std::binary_function<double, const MEMBER_TYPE&, const GENOTYPE&>
{
  public:
  
    typedef BinaryPenetranceCalculator<MEMBER_TYPE, GENOTYPE> SelfType;
  
  /// \name Object Management
  //@{
    template <typename SUSC, typename AFF>
    BinaryPenetranceCalculator(const SUSC& susc_func,
                               const AFF&  aff_func);
                               
    BinaryPenetranceCalculator(const SelfType& s);
    
    ~BinaryPenetranceCalculator();
    
    SelfType& operator=(const SelfType& s);
  //@}
    
    double operator() (const MEMBER_TYPE& m, const GENOTYPE& g) const;

  private:
  
    typedef boost::function<double (const MEMBER_TYPE&, const GENOTYPE&)> SuscFuncType;
    typedef boost::function<double (const MEMBER_TYPE&)>                  AffFuncType;

    SuscFuncType my_susc_func; ///< Susceptibility Functor
    AffFuncType  my_aff_func;  ///< Affection Functor
};

}
}

#include "pedcalc/binary_penetrance_calculator.ipp"

#endif

