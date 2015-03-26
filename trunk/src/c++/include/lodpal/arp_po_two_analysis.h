#ifndef LODPAL_ARP_PO_TWO_H
#define LODPAL_ARP_PO_TWO_H

//****************************************************************************
//* File:      arp_po_two_analysis.h                                         *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Feb. 03 *
//*                                                                          *
//* Notes:     This header file defines ARP_po_two_analysis class.           *
//*                                                                          *
//* Copyright (c) 2003 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_base_analysis.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     ARP_po_two_analysis                                          ~
// ~                                                                         ~
// ~ Purpose:   Defines two-parameter model                                  ~
// ~            ARP analysis for parent-of-origin test.                      ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ARP_po_two_analysis : public ARP_base_analysis
{
  public:

    ARP_po_two_analysis(lodpal_pairs& pairs, cerrorstream& err = sage_cerr);

  protected:

    virtual void   encode_params(size_t in, maxfun_param_mgr& pm);

    virtual void   decode_theta(const vector<double>& params_estimate,
                                const vector<double>& params_first_derivative);

    virtual double compute_lod_no_covariate  (vector<double>& theta, bool);
    virtual double compute_lod_with_covariate(vector<double>& theta, bool);

    virtual double compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval = false);

    virtual int    update_bounds_no_covariate  (vector<double>& theta);
    virtual int    update_bounds_with_covariate(vector<double>& theta);

    virtual void   compute_lambda(vector<double>& theta);
};

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
