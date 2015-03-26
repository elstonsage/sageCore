#ifndef  GENIBD_ANALYSIS_H
#define  GENIBD_ANALYSIS_H

//==========================================================================
//  File:       analysis.h
//
//  Author:     Geoff Wedig & Yeunjoo Song
//
//  History:    Version 1.0   Initial implementation
//                      2.0   Updated to new libraries         yjs  Nov. 03
//
//  Notes:      This class is the main driver performing analysis.
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "genibd/params.h"
#include "genibd/pair_ibd_analysis.h"
#include "genibd/sim_ibd_analysis.h"

namespace SAGE
{

namespace GENIBD
{

class genibd_analysis
{
  public:

    genibd_analysis(const RefMultiPedigree&   mp,
                    const genibd_parameters&  params,
                    genome_description*       genome,
                    ostream&                  info,
                    cerrorstream&             err = sage_cerr);

    ~genibd_analysis();

    bool run_analysis();

    const genibd_parameters*     get_parameters()    const;
    const RefMultiPedigree*      get_multipedigree() const;

  private:

    // ANALYSIS
    //
    bool  build();         //preparations before starting analysis.

    void  init_likelihood_data();

    bool  split_pedigree(const meiosis_map& mm);
    bool  use_single_point(const meiosis_map& mm);
    bool  use_exact(const meiosis_map& mm);
    bool  use_simulation(const meiosis_map& mm);

    bool  check_pedigree(const meiosis_map& mm, const region_type& r, analysis_data& data);

    size_t dump_header(size_t r_title);
    void   dump_bad_pedigrees(size_t r_title);

    bool  do_analysis();

    void  process_pedigree(const string&      title,
                           const string&      output,
                           const RefPedigree& rped,
                           const region_type& r);

    void  process_subpedigree(const string&         title,
                              const string&         output,
                              const string&         name,
                              const meiosis_map&    mm,
                              const region_type&    r,
                              const analysis_data&  adata);

    bool  do_exact_genibd_analysis     (const string&         title,
                                        const meiosis_map&    mm,
                                        const region_type&    r);

    bool  do_single_genibd_analysis    (const string&         title,
                                        const meiosis_map&    mm,
                                        const region_type&    r);

    bool  do_simulation_genibd_analysis(const string&         title,
                                        const meiosis_map&    mm,
                                        const region_type&    r);

    void  add_pairs_exact      (const subped_type& fsubped);
    void  add_pairs_single     (const subped_type& fsubped);
    void  add_pairs_simulation (const subped_type& fsubped);
    void  generate_pairs       (const subped_type& fsubped);

    IBD*  expand_intervals(IBD*, const meiosis_map& mmap);

    //
    // MEMBER DATA
    //
    const RefMultiPedigree*           my_multipedigree;
    const genibd_parameters*          my_parameters;
    genome_description*               my_genome;

    pair_ibd_analysis*                my_pair_ibd_analysis;
    exact_ibd_analysis*               my_exact_ibd_analysis;
    sim_ibd_analysis*                 my_sim_ibd_analysis;

    pedigree_region                   my_ped_region;
    
    relpair_set_type                  my_relpair_set;

    unsigned int                      pair_types;

    RefIBDWriteFile*                  my_ibd_prob_file;
    RefIBDWriteFile*                  my_ibd_state_file;

    vector<string>                    my_bad_fam;
    vector<string>                    my_bad_ped;

    size_t                            my_bad_fam_name_size;
    size_t                            my_bad_ped_name_size;

    ostream&                          info;
    cerrorstream                      errors;
};

#include "genibd/analysis.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
