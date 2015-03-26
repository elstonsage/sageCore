#ifndef POLYGENIC_TRANSITION_CALC_H
#define POLYGENIC_TRANSITION_CALC_H
//===========================================================================
//
//  File:	polygenic_transition_calculator.h
//
//  Author:	Stephen Gross
//
//  Purpose:	Calculate polygenic transition values.
//
//  Copyright (c) 2001, R. C. Elston
//
//===========================================================================


#include "numerics/binomial_dist.h"
#include "segreg/sub_model_base.h"
#include "segreg/model.h"
#include "segreg/segreg_datatypes.h"

namespace SAGE {
namespace SEGREG {

//============================================================================ 
// Class:	polygenic_transition_calculator
//============================================================================ 
class polygenic_transition_calculator 
{
  //==========================================================================
  // Constructor
  //==========================================================================
  public:
    polygenic_transition_calculator(const model & mod_param);

  //==========================================================================
  // Basic public get_transition function: prob(...)
  //==========================================================================
    double prob(size_t polygenotype_indiv,  size_t polygenotype_mother, 
                size_t polygenotype_father, size_t num_of_locii);

  //==========================================================================
  // Functions to calculate probability
  //==========================================================================
  private:
    double int_prob(size_t polygenotype_indiv,  size_t polygenotype_mother, 
                    size_t polygenotype_father, size_t num_of_locii);

    double calculate_prob(size_t polygenotype_indiv,  size_t polygenotype_mother, 
                          size_t polygenotype_father, size_t num_of_locii);

    double bin_trans(size_t genotype,         size_t new_polygenotype, 
                     size_t old_polygenotype, size_t num_of_locii);

  //==========================================================================
  // 3-d array to hold transition probabilities
  //==========================================================================
    vector<vector<vector<vector<double> > > > probs;

    const model & mod;

    size_t get_num_of_locii();

    //lint --e{1712}
};

}} // End namespace

#endif

