#ifndef LODPAL_ANALYSIS_H
#define LODPAL_ANALYSIS_H

//****************************************************************************
//* File:      lodpal_analysis.h                                             *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Initial implementation                                        *
//*                                                                          *
//* Notes:     This header file defines analysis driver class for lodpal.    *
//*                                                                          *
//* Copyright (c) 2006 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_parser.h"
#include "lodpal/lodpal_out.h"
#include "lodpal/arp_one_analysis.h"
#include "lodpal/arp_two_analysis.h"
#include "lodpal/arp_po_one_analysis.h"
#include "lodpal/arp_po_two_analysis.h"
#include "lodpal/arp_x_one_analysis.h"
#include "lodpal/arp_x_two_analysis.h"
#include "lodpal/dsp_one_analysis.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     lodpal_analysis                                              ~
// ~                                                                         ~
// ~ Purpose:   Defines analysis driver class for lodpal test                ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class lodpal_analysis
{
  public:

    lodpal_analysis(cerrorstream& err = sage_cerr);

    void     do_lodpal_test(RelativePairs &pairs, lodpal_parser& parser,
                            ostream &out, ostream &dig, ostream &xln);


    bool     try_run();
    void     run();

  protected:

    void     invalidate();

    // Temporary storage of intermediate results
    struct inter_results
    {
      inter_results() {}
      
      maxfun_results         maxfun_info;

      vector<double>         param_estimates;

      char                   method;

      double                 lod_score_with_cap;
      double                 lod_score_without_cap;

      double                 lod_score_mm;
      double                 lod_score_mf;
      double                 lod_score_ff;
    };

    void                run_model();
    inter_results       run_maxfun(maxfun_seq_cfg& sc, maxfun_param_mgr& i_theta);
    inter_results       run_method_A(maxfun_param_mgr& i_theta);
    inter_results       run_method_B(maxfun_param_mgr& i_theta);
    inter_results       run_method_D(maxfun_param_mgr& i_theta);

    void                encode_maxfun_params(vector<maxfun_param_mgr>& init_theta);
    
    void                decode_theta(const maxfun_results& re, vector<double>& params_estimate);
    void                decode_final_result(const inter_results& final_result);

    // members
    //
    typedef boost::shared_ptr<ARP_base_analysis> analysis_ptr;
    analysis_ptr              arp;

    bool                      my_built;
    size_t                    my_method_D_called;


    // Error stream
    cerrorstream              errors;
};

//#include "lodpal/lodpal_analysis.ipp"

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
