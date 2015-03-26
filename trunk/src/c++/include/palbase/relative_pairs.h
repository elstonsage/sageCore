#ifndef RELATIVE_PAIRS_H
#define RELATIVE_PAIRS_H

//****************************************************************************
//* File:      relative_pairs.h                                              *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Initial implementation                            kbj         *
//*            Pair information added.                           yjs Jan. 02 *
//*            X-linkage added.                                  yjs May. 02 *
//*                                                                          *
//* Notes:     This header file defines relative_pairs class to store ibd.   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "palbase/pair_info.h"

namespace SAGE    {
namespace PALBASE {

class relative_pairs
{
  public:

    relative_pairs();
    ~relative_pairs();

    // Validity tests
    bool      valid() const;
    void      invalidate();
    void      validate();

    // Set the pedigree
    void set_fped(FPED::Multipedigree& fp);

    // Pedigree data accessors
    const FPED::Multipedigree& get_fped() const;
          FPED::Multipedigree& get_fped();

    // Information stored for all pedigrees relating to traits and covariates
    const FPED::FilteredMultipedigreeInfo& fped_info() const;
          FPED::FilteredMultipedigreeInfo& fped_info();

    const FPED::FilteredPedigreeInfo&      ped_info(size_t p)  const;
          FPED::FilteredPedigreeInfo&      ped_info(size_t p);

    // Pedigree Data
    size_t    pedigree_count()                     const;
    string    pedigree_name(size_t p)              const;
    size_t    pedigree_number(size_t i)            const;

    size_t    member_count(size_t p)               const;
    string    member_name(size_t p, size_t i)      const;

    size_t    marker_count()                       const;
    size_t    marker_find(const std::string &name) const;

    string    marker_name(size_t m)                const;
    double    marker_distance(size_t m)            const;
    gmodel_type marker_genotype_model(size_t m)    const;

    double    get_average_marker_distance()        const;
    bool      valid_distance_exist()               const;

    size_t    trait_count()                        const;
    string    trait_name(size_t t)                 const;
    size_t    trait_find(const std::string &name)  const;

    // Total number of relative pairs
    size_t    pair_count()      const;
    size_t    fsib_pair_count() const;
    size_t    hsib_pair_count() const;
    size_t    mm_pair_count()   const;
    size_t    mf_pair_count()   const;
    size_t    ff_pair_count()   const;

    // Data access operations - single individual
    //    By pedigree, individual and trait index

    double    trait(size_t p, size_t i, size_t t)         const;
    double    trait_missing_code(size_t t)                const;
    bool      trait_missing(size_t p, size_t i, size_t t) const;

    double    trait(const MPED::member_base& m, size_t t)         const;
    bool      trait_missing(const MPED::member_base& m, size_t t) const;

    const rel_pair_data &rels(size_t i) const;

    void      set_f0(size_t i, size_t m, double f0);
    void      set_f2(size_t i, size_t m, double f2);
    void      set_f1mp(size_t i, size_t m, double f1);

    void      prebuild();
    void      build(bool build_sib_cluster = false, bool build_x_type = false);
    void      build_from_pedfile(const RPED::RefMultiPedigree& mped, bool build_sib_cluster = false, bool build_x_type = false);
    bool      built() const;

    bool      is_x_linked(size_t m)    const;
    bool      x_linked_marker_exist()  const;
    bool      autosomal_marker_exist() const;

    size_t    add_marker(const string& name, double dist, gmodel_type mt);

    size_t    add_pair(mem_pointer i1, mem_pointer i2, pair_type pt);

    size_t    find_pair(const mem_pointer& i1, const mem_pointer& i2) const;

    // Return the probability that pair at index i share n alleles IBD at
    // marker m.  Also return the average allele sharing.
    //
    // NOTE: Returns actual IBD sharing if available (or qNaN otherwise)

    double    prob_share(size_t i, size_t m, size_t n) const;
    double    avg_share(size_t i, size_t m)            const;
    double    weighted_share(size_t i, size_t m, double w0 = 0.0, double w1 = 0.5, double w2 = 1.0) const;

    double    prior_prob_share(size_t p, size_t i1, size_t i2,  size_t n) const;
    double    prior_prob_share(size_t i,  size_t n)                       const;
    double    prior_avg_share(size_t p, size_t i1, size_t i2)             const;
    double    prior_avg_share(size_t i)                                   const;
    double    prior_weighted_share(size_t i, double w0 = 0.0, double w1 = 0.5, double w2 = 1.0) const;

    sibship_cluster_iterator  sibship_cluster_begin();
    sibship_cluster_iterator  sibship_cluster_end();

    sibship_cluster_const_iterator  sibship_cluster_begin() const;
    sibship_cluster_const_iterator  sibship_cluster_end()   const;

    // For pair_info_file use
    //
    // Validity tests
    bool      valid_pair_info() const;
    void      invalidate_pair_info();
    void      validate_pair_info();

    size_t    pair_covariate_count() const;
    size_t    pair_weight_count()    const;

    size_t    pair_covariate_info_count() const;
    size_t    pair_weight_info_count()    const;

    int       set_pair_covariate(size_t pn, size_t c, const string& v, const pair_pheno_info& info);
    int       set_pair_weight   (size_t pn, size_t w, const string& v, const pair_pheno_info& info);

    double    get_pair_covariate(size_t pn, size_t c) const;
    double    get_pair_weight   (size_t pn, size_t w) const;

          pair_pheno_info& pair_covariate_info(size_t c);
    const pair_pheno_info& pair_covariate_info(size_t c) const;

          pair_pheno_info& pair_weight_info(size_t w);
    const pair_pheno_info& pair_weight_info(size_t w) const;
    
    string    get_pair_covariate_name(size_t c) const;
    string    get_pair_weight_name(size_t w)    const;

    size_t    pair_covariate_find(const string& name) const;
    size_t    pair_weight_find(const string& name)    const;

    size_t    add_pair_covariate(const string& name, pair_pheno_info::info_use usage = pair_pheno_info::mean,
                                 double value = std::numeric_limits<double>::quiet_NaN());

    size_t    add_pair_weight(const string& name, pair_pheno_info::info_use usage = pair_pheno_info::unknown,
                              double value = std::numeric_limits<double>::quiet_NaN());

    void      resize_pair_covariates(size_t c);
    void      resize_pair_weights(size_t w);

    const sped_pointer get_subped(size_t i) const;
    sped_pointer       get_subped(size_t i);

    bool      set_ibd_state(const sped_pointer p, const ibd_state_info& i_info);
    bool      get_ibd_state(const sped_pointer p, ibd_state_info& i_info) const;

    void      dump_pairs(ostream &out)             const;
    void      dump_sibship_cluster(ostream &out)   const;
    void      dump_pairs_by_pedigree(ostream &out) const;

  protected:

    pair_type     get_pair_type    (const mem_pointer i1, const mem_pointer i2) const;
    pair_sex_type get_pair_sex_type(const mem_pointer i1, const mem_pointer i2) const;
    pair_x_type   get_pair_x_type  (const mem_pointer i1, const mem_pointer i2) const;

    FPED::Multipedigree*                 my_fped;

    size_t                               my_fsib_pair_count;
    size_t                               my_hsib_pair_count;

    size_t                               my_mm_pair_count;
    size_t                               my_mf_pair_count;
    size_t                               my_ff_pair_count;

    bool                                 my_valid;
    bool                                 my_valid_distance_exist;
    bool                                 my_valid_pair_info;

    rel_map                              my_relmap;

    vector<RefPriorIBD>                  my_prior_ibd;

    vector<ibd_marker_info>              my_markers;
    vector<ibd_pair_info>                my_pairs;
    vector<ibd_probability_info>         my_ibd_probs;

    subped_ibd_state_map                 my_ibd_states;

    vector<pair_pheno_info>              my_pair_covariate_info;
    vector<pair_pheno_info>              my_pair_weight_info;
    
    vector< vector <double> >            my_pair_covariates;
    vector< vector <double> >            my_pair_weights;
  
    size_t                               my_pair_covariate_count;
    size_t                               my_pair_weight_count;

    sibship_cluster                      my_sibship_cluster;
};

#include "relative_pairs.ipp"

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
