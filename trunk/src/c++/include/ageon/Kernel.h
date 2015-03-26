#ifndef AO_KERNEL_H
#define AO_KERNEL_H
//======================================================================
//
//  File:	Kernel.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//
//  A sample likelihood is calculated as sum of individual log likelihoods.
//  Accordingly, the get_sample_likelihood(...) function sets in motion
//  a series of individual likelihood calculations. Each of the three
//  individual likelihood equations are relatively simple. They all fetch
//  the individual's pre-calculated mean, variance, susceptibility, age of
//  onset, age at exam, lambda1, and lambda2.
//
//======================================================================


#include <string>
#include <list>
#include <vector>
#include <cmath>
#include "numerics/cephes.h"
#include "numerics/functions.h"
#include "rped/rped.h"
#include "sampling/sampling.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "ageon/Datatypes.h"
#include "ageon/MemberCovariateCalculator.h"
#include "ageon/Model.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
//  class Kernel
//
//======================================================================
class Kernel
{
  public:
    //==================================================================
    // Constructor:
    //==================================================================

    Kernel(const Model & mod, const SAMPLING::PartitionedMemberDataSample & sample, int t);

    //==================================================================
    // Public utility functions:
    //==================================================================

    int update(); // update() is called by the ao_calculator before
                  // get_sample_likelihood() is fetched. update() takes
                  // subsequently calls all the other update() functions
                  // in the other data members.

    //==================================================================
    // Public accessors:
    //==================================================================

    double get_sample_likelihood(); // get_sample_likelihood() is the
                                    // most important public function available
                                    // in the Kernel. It calculates, and
                                    // returns, the overall sample likelihood.

    const MemberCovariateCalculator & get_mcc() const;

  private:
    //==================================================================
    // Private utility functions:
    //==================================================================

    double std_nor_density (double x, double stdev);

    double sign(double x);

    log_double get_indiv_likelihood                     (size_t indiv);
    log_double get_indiv_likelihood_affected_AO_known   (size_t indiv);
    log_double get_indiv_likelihood_affected_AO_unknown (size_t indiv);
    log_double get_indiv_likelihood_unaffected          (size_t indiv);

    void update_truncation_denominator(double new_denominator);

    double transform(double);

    //==================================================================
    // Data members:
    //==================================================================

          int                                     my_analysis_type;
          MemberCovariateCalculator               my_mcc;
    const Model                                 & my_model;
    const SAMPLING::PartitionedMemberDataSample & my_sample;
};

//======================================================================
//
//  std_nor_density(...)
//
//======================================================================
inline double
Kernel::std_nor_density(double x, double stdev)
{
  return ONE_OVER_TWO_PI * exp(-0.5 * x * x / stdev);
}

//======================================================================
//
//  update_truncation_denominator(...)
//
//======================================================================
inline void 
Kernel::update_truncation_denominator(double new_denominator)
{
  if(new_denominator < my_model.get_min_denominator())
    my_model.min_denominator() = new_denominator;
}

//======================================================================
//
//  sign(...)
//
//======================================================================
inline double 
Kernel::sign(double x)
{
  return x >= 0 ? 1.0 : -1.0;
}

//======================================================================
//
//  transform(...)
//
//======================================================================
inline double 
Kernel::transform(double x)
{
  double y = x;

  my_mcc.transform(y);

  return y;
}

}} // End namespace

#endif
