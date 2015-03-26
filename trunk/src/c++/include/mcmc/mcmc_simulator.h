#ifndef  MCMC_SIMULATOR_H
#define  MCMC_SIMULATOR_H

//==========================================================================
//  File:       mcmc_simulator.h
//
//  Author:     Qing Sun
//
//  History:    Version 0.90
//              Updated to new libraries                       - yjs Sep. 04
//
//  Notes:      This header defines MCMC simulator which walks on the space
//              based on inheritance vetors and data structures defined in
//              MCMC_DataStruct.h file.(See Soble and Lang, Am.,J. Hum
//              Genet. 58:1323-1337,1996)
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//==========================================================================

// The MCMC Simulator simulates descent state at each of several markers,
//   either single- or multi-point.  This is done by generating a series of
//   transitions and randomly walking over the available states.  A
//   secondary class must be derived from the MCMC_simulator class, which controls
//   program flow.  See IBD_MCMC_Simulator.h for an example of this.

// In general, processing of a step is done as follows:

//    1. Initialize by calling start_step().
//    2. Generate a transition.  This can be done with generate_transtition()
//           or with some other function.  The transtion bit pattern should
//           be stored in the state_tables[markers].local_bits_changed.
//    3. Apply the transition with apply_transition().  This updates the
//           current bit pattern (in the indicators) as well as calculating
//           the likelihood of the transition.  It keeps total_bits_changed
//           up to date and does not modify local_bits_changed.
//    4. Reject with reject_transition() (optional).  The reject function
//           unwinds the previously applied apply_transition().  Nothing
//           need done if the previous state is accepted.
//    5. End the step with end_step().


#include "numerics/cephes.h"
#include "containers/bit_remap.h"
#include "mcmc/mcmc_params.h"
#include "mcmc/starting_state.h"

namespace SAGE
{

namespace MCMC
{

#define  MAX_TRANS               99   //max tunnel through transitions
#define  SM_UTP_weight           10   //weight to choose untyped pivot person

class mcmc_simulator
{
  public:

    mcmc_simulator(const McmcMeiosisMap&   ped, const pedigree_region& pr,
                   mcmc_parameters*        par, cerrorstream&          err);
    virtual ~mcmc_simulator();

    virtual bool start(ostream& o);

    MersenneTwister& get_random()                     const;
    double           get_random_real()                const;
    int              get_random_uniform_int(size_t t) const;

  protected: 

    // transition rules:
    // T01            rule T0 left  (mother) node of the pivot person
    // T02            rule T0 right (father) node of the pivot person
    // T03            rule T0 both.  Both parental nodes are switched.
    // T1             rule T1 
    // T2a            rule T2a
    // T2b            rule T2b
    // ILLEGAL_SWITCH_RULE     
    //
    enum transition_rule { T01=0, T02, T03, T1, T2a, T2b, ILLEGAL_SWITCH_RULE };

    // The transition defines a single transition of the primary
    // types (T0, T1 or T2).
    //
    struct transition
    {
      transition() 
        : marker(0),
          rule(ILLEGAL_SWITCH_RULE),
          person(0),
          extent(0)
      { }

      int              marker;
      transition_rule  rule; 
      size_t           person;
      int              extent;
    };

    //Storage for typed and untyped individuals and offspring for each marker.
    //These are used to choose the pivot person in each transition.
    //
    typedef vector<size_t>   ind_vector;

    struct typed_individual_type
    {
      ind_vector typed_individual;
      ind_vector untyped_individual;

      double typed_probability;
    };

    struct local_individual_type
    {
      ind_vector local;
    };

    //member functions
    //

    // Initialization functions

    virtual
    bool   initialize_analysis(ostream&);
    int    get_valid_markers();

    void   init_storage();
    void   update_storage();

    void   init_local_individuals();

    bool   create_starting_state();

    bool   dememorize(dot_formatter& out, int demem_steps = 0);

    // Step funtions

    void   start_step();

    void   generate_transition(bool do_tunnel = true);
    double apply_transition();
    void   reject_transition();

    void   end_step();

    // Step helper functions 

    //   Choose element

    void   select_a_locus            (transition& t, bool& n) const;
    void   select_local_pivot_person (transition& t)          const;
    void   select_pivot_person       (transition& t)          const;
    void   select_trans_rule         (transition& t)          const;
    void   select_extent             (transition& t)          const;

    //   Run a transition

    void   T0_switch(const transition& t);
    void   T1_switch(const transition& t);
    void   T2_switch(const transition& t);

    void   change_bit(size_t i, bool mother);

    void   apply_bit_pattern_at_marker(int m);

    //   Likelihood generating

    double   log_recomb_prob_ratio() const;

    // Debugging functions

    void   dump_nuc_fams(ostream& o) const;
    void   dump_marker_inds(ostream& o) const;
    void   dump_storage() const;   

    void   print_descent_pattern(size_t m) const;

    //double calculate_inheritance_around_marker(int);
    //double calculate_inheritance_at_interval(int);

    //member data
    //
    const McmcMeiosisMap&                 my_pedigree;
    const pedigree_region&                my_ped_region;
    mcmc_parameters*                      my_parameters;

    mcmc_data_accessor*                   my_data;
    recombination_calculator*             my_recomb_like;
    marker_likelihood_calculator*         my_marker_like;
    starting_state_generator*             my_starting_state;

    mutable MersenneTwister               my_random_generator;

    bool                                  my_loop;
    bool                                  my_multipoint;

    size_t                                my_total_loci;
    size_t                                my_total_useful_loci;
    size_t                                my_current_marker;     //used for single point 

    int                                   my_transition_steps;   //step of tunnel through of last walk

    //book keeping data  
    int                                   my_last_switch_marker;
    int                                   my_last_pivot_person;

    int                                   my_old_last_switch;
    int                                   my_old_last_pivot;

    //The following parameters can be obtained from my_param. 
    //For the sake of efficiency, make the extra copies here since 
    //they are accessed millions times during the simulation.  
    //
    int                                   my_max_tunnel; 
    double                                my_single_marker_ratio;

    vector<double>                        my_log_theta_ratio;    // Ratios of theta probabilities on log scale

    vector<typed_individual_type>         my_marker_individuals; // For T0 rule - typed & untyped ind vectors per marker

    vector<vector<size_t> >               my_nuclear_families;   // For T2 - Nuclear families by first sib.
    vector<local_individual_type>         my_local_individuals;

    // For building a transition
    typedef bit_remap<bit_field, unsigned char> bit_remap_type;

    vector<int>                           my_markers_used;
    vector<bit_field>                     my_bits_changed;
    vector<bit_remap_type>                my_bit_remappers;

    bit_field                             my_transition_pattern;

    cerrorstream                          errors;
};

#include "mcmc/mcmc_simulator.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
