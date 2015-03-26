#ifndef PAIR_IBD_ANALYSIS_H
#define PAIR_IBD_ANALYSIS_H

//====================================================================
//  File:     pair_ibd_analysis.h
//
//  Author:   Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:  0.1  gcw Relpal's Generator              until July 1998
//            0.2  gcw Adapted for GENIBD              Aug/Sept   1998
//            1.0  yjs Updated to new libraries.             Jan. 2004
//
//  Notes:    Pairwise IBD Analysis - Single point Relative Pair IBD
//             generator.
//
//  Copyright (c) 1998  R.C. Elston
//====================================================================

#include "genibd/peeler.h"

namespace SAGE
{

namespace GENIBD
{

class pair_ibd_analysis
{
  public:

    pair_ibd_analysis(cerrorstream& = sage_cerr);
    
    ~pair_ibd_analysis();
    
    cerrorstream get_errors()  const;
    cerrorstream set_errors(cerrorstream& s);

    bool         built()       const;
    bool         valid()       const;

    IBD*         ibd_adaptor() const;

    bool         build();

    bool         set_pedigree(const meiosis_map& mmap, const pedigree_region&);

    bool         build_ibds();

    bool         add_pair(fmember_const_pointer m1, fmember_const_pointer m2,
                          fmember_const_pointer c1, fmember_const_pointer c2,
                          pair_generator::pair_type t);

    bool         compute(const string& title = string());

  protected:

    void compute_subpedigree_likelihood(const inheritance_model& model);

    void generate_pi_hats(const inheritance_model& model);

    void generate_sib    (size_t current_pair, const inheritance_model& model);
    void generate_hsib   (size_t current_pair, const inheritance_model& model);
    void generate_avunc  (size_t current_pair, const inheritance_model& model);
    void generate_gpar   (size_t current_pair, const inheritance_model& model);
    void generate_cous   (size_t current_pair, const inheritance_model& model);

    double cond_post (fmember_const_pointer    mp,
                      fmember_const_pointer    spouse,
                      fmember_const_pointer    child,
                      const MLOCUS::penetrance_model::phased_penetrance_iterator&      mp_iter,
                      const MLOCUS::penetrance_model::phased_penetrance_iterator&      child_iter,
                      const allele&            al,
                      const inheritance_model& model) const;

    bool   ch_geno1  (const phased_genotype& pg1,
                      const phased_genotype& pg2,
                      bool                   first_is_mother,
                      const phased_genotype& cg,
                      const allele&          al,
                      conditional_genotype&  con_geno) const;

    bool   condition (const phased_genotype& g1,
                      const phased_genotype& g2,
                      allele*                al) const;

    // transmission probabilty conditional upon shared alleles. 
    double cond_tran (const phased_genotype& g, const conditional_genotype& cg) const;

    // Number of alleles of type a in genotype g.
    size_t has_allele(const phased_genotype& g, const allele& al) const;

    // Member
    //
    relative_pairs              my_pairs;

    size_t                      my_current_marker;

    double                      my_subped_likelihood;

    meiosis_map                 my_meiosis_map;

    region_type                 my_region;
    pedigree_region             my_ped_region;

    peeler*                     my_peeler;
    IBD*                        my_ibds;

    text_dot_formatter*         my_dots;

    bool                        my_built;
    bool                        my_valid;
    bool                        my_verbose;

    cerrorstream                errors;
};

#include "genibd/pair_ibd_analysis.ipp"

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
