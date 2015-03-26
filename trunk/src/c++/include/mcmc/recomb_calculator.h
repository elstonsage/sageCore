#ifndef RECOMB_CALCULATOR_H
#define RECOMB_CALCULATOR_H

//==========================================================================
//  File:     recomb_calculator.h
//
//  Author:   Geoff Wedig
//
//  History:  0.1 Initial Implementation
//            1.0 Updated to new libraries                       yjs May. 04
//
//  Notes:    Calculates the inheritance probabilities of intervals
//            between markers based upon Lander-Green Inheritance patterns
//            observed at those markers.  It can also claculate the change
//            in ratio when those patterns change.
//
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/mcmc_data_accessor.h"


namespace SAGE
{

namespace MCMC
{

class recombination_calculator
{
  public:

    recombination_calculator(const pedigree_region& pr, const mcmc_data_accessor& data);

    // Calculation funtions
    
    double operator()       (size_t marker) const; // log inheritance prob. between m and m+1
    double log_recombination(size_t marker) const; // log inheritance prob. between m and m+1

    double operator()       ()              const; // log inheritance prob. for all markers
    double log_recombination()              const; // log inheritance prob. for all markers

    double operator()       (size_t m, const bit_field& b1, const bit_field& b2) const;
    double log_recombination(size_t m, const bit_field& b1, const bit_field& b2) const;

    // Ratios - Ratio of before and after the bits in b1 and b2 are applied.
    //          Useful in calculating the change between successive steps in
    //          the mcmc simulation.

    double log_recombination_ratio(size_t m, const bit_field& b1, const bit_field& b2) const;

    void   dump_recombination_calculator(ostream& o) const;

  protected:

    typedef vector<double> theta_ratio_type;

    const pedigree_region&       my_ped_region;
    const mcmc_data_accessor&    my_data;

    theta_ratio_type             my_log_thetas;
    theta_ratio_type             my_log_one_minus_thetas;
    theta_ratio_type             my_theta_ratios;
};


#include "mcmc/recomb_calculator.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
