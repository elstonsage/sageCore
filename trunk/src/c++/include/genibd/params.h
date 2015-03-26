#ifndef GENIBD_PARAMS_H
#define GENIBD_PARAMS_H

//========================================================================== 
//  File:     GENIBD parameters
//
//  Author:   Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:  0.01 gcw Initial Interface design                 Jun 15, 1998
//            0.1  gcw Release version                          Jul  8, 1998
//                 yjs Added old_ibd_format option              Aug 30, 2002
//            2.0  yjs Updated to new libraries.                Nov.    2004
//
//  Copyright (c) 1998 R. C. Elston
//  All Rights Reserved
//========================================================================== 

#include "genibd/definitions.h"

#define GENIBD_CUTOFF        18
#define GENIBD_ANALYSIS_NAME "IBD_ANALYSIS"

namespace SAGE
{

namespace GENIBD
{

class genibd_parameters
{
  public:

    genibd_parameters();

    genibd_parameters(const genibd_parameters&);
    
    genibd_parameters& operator=(const genibd_parameters&);

    ~genibd_parameters();

    const string&              title()                     const;
    const string&              output()                    const;

    size_type                  max_exact_size()            const;
    bool                       is_multipoint()             const; 
    bool                       scan_interval()             const;
    double                     interval_distance()         const;
    bool                       allow_loops()               const;
    bool                       output_ibd_state()          const;
    choice_type                allow_simulation()          const;
    choice_type                allow_family_splitting()    const;
    pair_category_type         pair_category()             const;

    genibd_region_iterator     region_begin()              const;
    genibd_region_iterator     region_end()                const;

    size_type                  region_count()              const;

          mcmc_parameters& get_sim_parameters();
    const mcmc_parameters& get_sim_parameters()        const;

    // Modification functions

    void set_title(const string&);
    void set_output(const string&);

    void set_max_exact_size(long);
    void set_multipoint(bool);
    void set_scan_interval(bool);
    void set_interval_distance(double);
    void set_allow_loops(bool);
    void set_output_ibd_state(bool);
    void set_allow_simulation(choice_type);
    void set_allow_family_splitting(choice_type);
    void set_pair_category(pair_category_type);

    void build_analysis_region(genome_description* gd, const genibd_region_list& gl,
                               cerrorstream&  err = sage_cerr);

    void add_region(const string&, const string& = string());
    void remove_region(const string&);
                                       
    // Printout functions
    
    void dump_parameters(ostream& out = cout) const;

  protected:

    void dump_multipoint(SAGE::cerrorstream& out) const;
    void dump_singlepoint(SAGE::cerrorstream& out) const;

    void dump_title(SAGE::cerrorstream& out) const;
    void dump_regions(SAGE::cerrorstream& out) const;

    string               my_title;
    string               my_output;

    size_t               my_exact_size;
    bool                 my_multipoint;
    bool                 my_scan_interval;
    double               my_interval_distance;
    bool                 my_loops;
    bool                 my_output_ibd_state;
    choice_type          my_simulation;
    choice_type          my_family_split;
    pair_category_type   my_pair_type;

    genibd_region_list   my_regions;

    mutable mcmc_parameters my_sim_parameters;
};

#include "genibd/params.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
