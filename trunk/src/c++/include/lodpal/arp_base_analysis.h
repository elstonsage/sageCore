#ifndef LODPAL_ARP_BASE_H
#define LODPAL_ARP_BASE_H

//****************************************************************************
//* File:      ARP_base_analysis.h                                             *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                kbj         *
//*                    1.0 one-parameter model added.            yjs Nov. 00 *
//*                    1.1 covariate added.                      yjs Nov. 00 *
//*                    1.2 bad_sib_pair storage & func added.    yjs Mar. 01 *
//*                    1.3 rel_pair_map storage & func added.    yjs Mar. 01 * 
//*                    1.3 multipoint, singlepoint seperated.    yjs Apr. 01 *
//*                    1.4 dsp, re-parameterization added.       yjs Apr. 01 *
//*                    1.5 evaluate & update_bound for dsp added.yjs May. 01 *
//*                    1.6 diagnostic option added.              yjs Jun. 01 *
//*                    1.7 seperation of pair & analysis.        yjs Jul. 01 *
//*                    1.8 one-parameter model as default.                   *
//*                        two-parameter with covariate removed. yjs Aug. 01 *
//*                    2.0 parent-of-origin added.               yjs Feb. 03 *
//*                                                                          *
//* Notes:     This header file defines base analysis class for lodpal.      *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_result.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     ARP_base_analysis                                            ~
// ~                                                                         ~
// ~ Purpose:   Defines base analysis class for lodpal test                  ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ARP_base_analysis : public MaxFunction
{
  public:

    ARP_base_analysis(lodpal_pairs& l_pairs, cerrorstream& err = sage_cerr);

    void  set_parameters(const lodpal_parameters& p);
    void  set_previous_parameters(const lodpal_parameters& p);

    void  set_good_bound(bool b = true);
    bool  is_good_bound() const;

          lodpal_parameters& parameters();
    const lodpal_parameters& parameters()              const;
          lodpal_parameters& previous_parameters();
    const lodpal_parameters& previous_parameters()     const;

          lodpal_pairs&      pairs_info();
    const lodpal_pairs&      pairs_info()              const;

    const RelativePairs&     relative_pairs()          const;

    size_t                   valid_parameter_count()   const;
    bool                     built()                   const;

    double                   get_lod_mm()              const;
    double                   get_lod_mf()              const;
    double                   get_lod_ff()              const;

          lodpal_result&     get_lodpal_result();
    const lodpal_result&     get_lodpal_result()       const;

    void                invalidate();
    void                build();
    void                re_estimate_parameters();

    const marker_type&  current_marker() const;
    void                set_current_marker(const marker_type& m);

    double              re_evaluate(const vector<double>& theta);

    void                update_result(const vector<double>& param_est, size_t removed_pos);

    virtual void        encode_params(size_t in, maxfun_param_mgr& pm) = 0;

    virtual void        decode_theta(const vector<double>& params_estimate,
                                     const vector<double>& params_first_derivative) = 0;

  protected:

    virtual double      evaluate(vector<double>& theta);
    virtual double      compute_lod_no_covariate  (vector<double>& theta, bool) = 0;
    virtual double      compute_lod_with_covariate(vector<double>& theta, bool) = 0;

    virtual double      compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval = false) = 0;

    virtual int         update_bounds(vector<double>& theta);
    virtual int         update_bounds_no_covariate  (vector<double>& theta) = 0;
    virtual int         update_bounds_with_covariate(vector<double>& theta) = 0;

    virtual void        compute_lambda(vector<double>& theta) = 0;

    // members
    //
    lodpal_parameters         my_parameters;
    lodpal_parameters         my_previous_parameters;

    lodpal_pairs              my_pairs;
    
    marker_type               my_current_marker;

    size_t                    my_valid_parameter_count;

    bool                      my_built;
    bool                      my_good_bound;

    // Results from each Maxfun routine
    //  since Maxfun only returns total_max_lod.
    double                    my_ld_mm;
    double                    my_ld_mf;
    double                    my_ld_ff;

    // Final result
    lodpal_result             my_result;

    cerrorstream              errors;
};

#include "lodpal/arp_base_analysis.ipp"

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
