//======================================================================
//
//  File:  Calculator.cpp
//
//  Author:  Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================

#include "ageon/Calculator.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
//  Calculator() CONSTRUCTOR
//
//======================================================================
Calculator::Calculator(
  Model                                       & mod, 
  const SAMPLING::PartitionedMemberDataSample & sample,
  int                                           t)
  :
  my_model         (mod),
  my_kernel        (mod, sample, t),
  my_analysis_type (t)
{ }

//======================================================================
//
//  evaluate(...)
//
//======================================================================
double
Calculator::evaluate(vector<double> & params)
{
  nfe++;

  double lh = my_kernel.get_sample_likelihood();

  return lh;
}

//======================================================================
//
//  update_bounds(...)
//
//======================================================================
int
Calculator::update_bounds(vector<double> & params)
{
  if(my_model  .update (params, my_analysis_type)) return 1;
  if(my_kernel .update ())                         return 1;

  return 0;
}

//======================================================================
//
//  get_mcc()
//
//======================================================================
const MemberCovariateCalculator &
Calculator::get_mcc() const
{
  return my_kernel.get_mcc();
}

} // End namespace AO
} // End namespace SAGE
