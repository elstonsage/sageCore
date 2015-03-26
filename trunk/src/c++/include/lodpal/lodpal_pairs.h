#ifndef LODPAL_PAIRS_H
#define LODPAL_PAIRS_H

//****************************************************************************
//* File:      lodpal_pairs.h                                                *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   1. Initial implementation                         yjs         *
//*            2. bad_sib_pair storage & func added.             yjs Mar. 01 *
//*            3. rel_pair_map storage & func added.             yjs Mar. 01 * 
//*            4. seperated from ARPTest.                        yjs Jul. 01 *
//*            5. weight added.                                  yjs Jan. 02 *
//*            6. x-linkage added.                               yjs Apr. 02 *
//*            7. contrast option added.                         yjs Feb. 06 *
//*                                                                          *
//* Notes:     This header file defines lodpal_pairs class.                  *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_params.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     lodpal_pairs                                                 ~
// ~                                                                         ~
// ~ Purpose:   Save pair-specific covariate values for constraints checking ~
// ~             of covariate parameters if covariate_count > 0.             ~
// ~            Re-designed to store individual pair's likelihood.           ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class lodpal_pairs
{
  public:

    lodpal_pairs(RelativePairs& p, lodpal_parameters& par, cerrorstream& err);
    lodpal_pairs(lodpal_pairs& p);

    ~lodpal_pairs();

    struct covariates_info
    {
      covariates_info() { pair_value       = numeric_limits<double>::quiet_NaN();
                          ad_pair_value    = numeric_limits<double>::quiet_NaN();
                          re_ad_pair_value = numeric_limits<double>::quiet_NaN(); }

      double       pair_value;
      double       ad_pair_value;
      double       re_ad_pair_value;
    };
      
    struct prior_x_ibd_type
    {
      prior_x_ibd_type() { subtype  = "";
                           pf0 = numeric_limits<double>::quiet_NaN();
                           pf1 = numeric_limits<double>::quiet_NaN();
                           pf2 = numeric_limits<double>::quiet_NaN(); }

      prior_x_ibd_type(string s, double f0, double f1, double f2)
      { subtype = s; pf0 = f0; pf1 = f1; pf2 = f2; }

      string        subtype;
      double        pf0;
      double        pf1;
      double        pf2;
    };

    struct lodpal_pair_info
    {
      lodpal_pair_info() { }
      lodpal_pair_info(const rel_pair p, const vector<covariates_info> c,
                       double w  = 1.0,
                       bool   af = true,
                       size_t type = 0,
                       double l  = numeric_limits<double>::quiet_NaN(),
                       double f  = numeric_limits<double>::quiet_NaN(),
                       double f1 = numeric_limits<double>::quiet_NaN(),
                       double ff = numeric_limits<double>::quiet_NaN(),
                       bool   r  = false,
                       bool   rx = false,
                       size_t x  = (size_t)-1)

      {lodpal_pair=p; lodpal_cov=c; lodpal_weight=w; concordantly_affected=af; pair_type = type;
       likelihood=l; f0=f; f1mp=f1; f2=ff; removed=r; removed_x = rx; prior_x_ibd_index = x;}
      
      rel_pair                    lodpal_pair;
      vector<covariates_info>     lodpal_cov;

      double                      lodpal_weight;
      
      double                      f0;
      double                      f1mp;
      double                      f2;
      double                      likelihood;

      bool                        removed;   // Removed (temp at a point) in the process of maximization.
      bool                        removed_x; // Not included at all for the x-linkage analysis.

      bool                        concordantly_affected;

      size_t                      pair_type;

      size_t                      prior_x_ibd_index;
    };

    typedef  vector<lodpal_pair_info>            pair_info_type;
    typedef  vector<prior_x_ibd_type>            prior_x_ibd_vector;
    typedef  vector<pair<double, double> >       drp_ibd_info;

    // Storage for each valid pair likelihood from previous point.
    typedef  std::map<rel_pair, double, PALBASE::rel_pair_less>  pair_map_type;

    void     invalidate_build_pairs_info();

    void     build_pairs_info();
    void     build_pairs_info_x();
    void     build_pairs_map();
    void     re_build_pairs_info();
    
    size_t   remove_biggest_pair(bool x = false);
    void     reset_removed_pairs();

    // Generate filter for rel pair iteration.
    void           make_filter();
    pair_filter&   filter();

    // Get the trait value for an individual.                                                    
    double member_trait(const MPED::member_base* m, size_t t) const;                      

    // Compute the value of a covariate/weight for a pair.
    double covariate_value(const rel_pair& pair, const covariate_type& param)     const; 
    double weight_value   (const rel_pair& pair, const lodpal_weight_type& param) const; 

    bool     built_pairs_info()                   const;
    bool     built_pairs_info_x()                 const;
    bool     re_built_pairs_info()                const;
    bool     removed_biggest_pair()               const;
 
    bool     is_mm_pair(size_t prior_index)       const;
    bool     is_mf_pair(size_t prior_index)       const;
    bool     is_ff_pair(size_t prior_index)       const;
    bool     is_ff_sib_pair(size_t prior_index)   const;

    bool     parent_of_origin_allowed()           const;

    size_t   re_built_covariate()                 const;

    size_t   pair_count()                         const;
    size_t   fsib_pair_count()                    const;
    size_t   hsib_pair_count()                    const;
    size_t   other_pair_count()                   const;

    size_t   mm_pair_count()                      const;
    size_t   mm_fsib_pair_count()                 const;
    size_t   mfm_hsib_pair_count()                const;
    size_t   mm_other_pair_count()                const;
    size_t   mm_invalid_pair_count()              const;

    size_t   mf_pair_count()                      const;
    size_t   mf_fsib_pair_count()                 const;
    size_t   mff_hsib_pair_count()                const;
    size_t   mf_other_pair_count()                const;
    size_t   mf_invalid_pair_count()              const;

    size_t   ff_pair_count()                      const;
    size_t   ff_fsib_pair_count()                 const;
    size_t   fff_hsib_pair_count()                const;
    size_t   ff_other_pair_count()                const;
    size_t   ff_invalid_pair_count()              const;

    size_t   max_ped_name_size()                  const;
    size_t   max_ind_name_size()                  const;

    double   get_drp_prob_share(size_t pt, size_t f) const;
    double   get_drp_x_prob_share(size_t pt, size_t f) const;
    
          pair_info_type&      pairs_info();
    const pair_info_type&      pairs_info()       const;
          pair_map_type&       pairs_map();
    const pair_map_type&       pairs_map()        const;

          RelativePairs&       relative_pairs();
    const RelativePairs&       relative_pairs()   const;
          lodpal_parameters&   parameters();
    const lodpal_parameters&   parameters()       const;

          vector<size_t>&      removed_fsib_pairs();
    const vector<size_t>&      removed_fsib_pairs() const;
    
          vector<size_t>&      removed_hsib_pairs();
    const vector<size_t>&      removed_hsib_pairs() const;

          vector<size_t>&      removed_other_pairs();
    const vector<size_t>&      removed_other_pairs() const;

          prior_x_ibd_vector&  prior_x_ibds();
    const prior_x_ibd_vector&  prior_x_ibds() const;

    cerrorstream&              errs(){return errors;}

  protected:

    size_t   get_pair_type(const rel_pair& pair) const;
    size_t   get_x_pair_type(const rel_pair& pair) const;

    RelativePairs&            pairs;
    lodpal_parameters&        params;
    pair_filter               my_filter;
    
    size_t                    my_pair_count;
    size_t                    my_fsib_pair_count;
    size_t                    my_hsib_pair_count;

    size_t                    my_mm_pair_count;
    size_t                    my_mm_fsib_pair_count;
    size_t                    my_mfm_hsib_pair_count;
    size_t                    my_mm_invalid_pair_count;

    size_t                    my_mf_pair_count;
    size_t                    my_mf_fsib_pair_count;
    size_t                    my_mff_hsib_pair_count;
    size_t                    my_mf_invalid_pair_count;

    size_t                    my_ff_pair_count;
    size_t                    my_ff_fsib_pair_count;
    size_t                    my_fff_hsib_pair_count;
    size_t                    my_ff_invalid_pair_count;

    size_t                    my_max_ped_name;
    size_t                    my_max_ind_name;

    pair_info_type            my_pairs_info;
    pair_map_type             my_pairs_map;
    
    bool                      my_built_pairs_info;
    bool                      my_built_pairs_info_x;
    bool                      my_re_built_pairs_info;

    bool                      my_parent_of_origin_allowed;

    size_t                    my_re_built_covariate;

    vector<size_t>            my_removed_fsib_pairs;
    vector<size_t>            my_removed_hsib_pairs;
    vector<size_t>            my_removed_other_pairs;

    prior_x_ibd_vector        my_prior_x_ibd_vector;

    drp_ibd_info              my_drp_ibd_info;
    drp_ibd_info              my_drp_x_ibd_info;
/*
    pair_info_type            my_lodpal_pairs_concordant;
    pair_info_type            my_lodpal_pairs_discordant;
    
    bool                      my_seperation;
    bool                      my_concordant_pair_maximized;
    bool                      my_discordant_pair_maximized;
*/
    cerrorstream&             errors;
};

#include "lodpal/lodpal_pairs.ipp"

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
