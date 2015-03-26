#ifndef  STARTING_STATE_H
#define  STARTING_STATE_H

//==========================================================================
//  File:     starting_state.h
//
//  Author:   Geoff Wedig                                                   
//
//  History:  0.1 Initial Implementation
//            1.0 Updated to new libraries                       yjs Jun. 04
//
//  Notes:
//                                                                          
//  Copyright (c) 2004 R. C. Elston
//  All Rights Reserved
//==========================================================================

#include "numerics/mt.h"
#include "mcmc/marker_likelihoods.h"
#include "mcmc/recomb_calculator.h"


namespace SAGE
{

namespace MCMC
{


class starting_state_generator
{
  public:

    starting_state_generator(const McmcMeiosisMap&               mmap,
                             const pedigree_region&              ped_r,
                             const marker_likelihood_calculator& m_like,
                             const recombination_calculator&     re_cal,
                             MersenneTwister&                    ra_gen);

    ~starting_state_generator();

    int set_indicator_count(int=10);       // Number of indicators generated for each marker
    int get_indicator_count() const;

    int set_random_trial_count(int=5);     // Number of random trials to attempt when looking for 
    int get_random_trial_count() const;    // a valid indicator for a given marker.  Only useful
                                           // with pedigrees with loops which might not be genoset
                                           // minimized.

    double generate_starting_state(mcmc_data_accessor::indicator_type& state);

    bool   is_valid() const;

  private:

    typedef vector<pair<CommonAlleleSet::AlleleID, CommonAlleleSet::AlleleID> > allele_vector;

    struct indicator_data
    {
      bit_field     my_indicator;
      double        log_likelihood;
      vector<int>   previous_markers;
    };

    typedef vector<indicator_data>    marker_indicators;
    typedef vector<marker_indicators> region_indicators;

    // Create a valid marker state for marker m, set log_likelihood to marker likelihood
    //
    void create_marker_state(indicator_data& id, int m, genotype_eliminator& ge);

    // Check previous marker states and choose one, set previous markers and recalc log like.
    //
    void choose_previous_marker(indicator_data& id, int m);

    // Given a vector of log likelihoods, choose one based on their relative
    // probabilities.  Return the index.
    //
    int choose_one(vector<double>& like) const;

    //member functions
    //
    void assign_bit        (const allele_vector& av, indicator_data& id);
    void assign_random_bit (indicator_data& id);

    // Reduce the genoset of a given pedigree to a single genotype per individual.
    // Randomization of which individuals are reduced first and which genotype of each
    // individual is chosen.  A genotype elimination step is used between each reduction
    // step to increate the probability of a successful genoset resulting.
    //
    bool random_genotype_reduction(allele_vector& av, const MLOCUS::inheritance_model& imodel, int m, genotype_eliminator& ge);

    //data members
    //
    const pedigree_region&               my_ped_region;
    const McmcMeiosisMap&                my_pedigree;
    const marker_likelihood_calculator&  my_marker_like;
    const recombination_calculator&      my_recomb_calculator;

    int                                  my_indicator_count;
    int                                  my_trial_count;

    region_indicators                    my_indicators;

    MersenneTwister&                     my_random_generator;

    bool                                 my_valid;
};

#include "mcmc/starting_state.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
