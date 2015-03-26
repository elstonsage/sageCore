#ifndef LODPAL_DSP_ONE_H
#define LODPAL_DSP_ONE_H

//****************************************************************************
//* File:      dsp_one_analysis.h                                            *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                yjs Apr. 01 *
//*                                                                          *
//* Notes:     This header file defines DSP_one_analysis class.              *
//*                                                                          *
//* Copyright (c) 2001 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/arp_one_analysis.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     DSP_one_analysis                                             ~
// ~                                                                         ~
// ~ Purpose:   Defines one-parameter model DSP(discordant sib pair) test.   ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class DSP_one_analysis : public ARP_one_analysis
{
  public:

    DSP_one_analysis(lodpal_pairs& pairs, cerrorstream& err = sage_cerr);

  protected:

    virtual void   encode_params(size_t in, maxfun_param_mgr& pm);

    virtual double compute_lod_with_covariate(const vector<double>& theta, bool);

    virtual double compute_a_log_lr(const vector<double>& theta, size_t pos, bool re_eval = false);

    virtual int    update_bounds_with_covariate(const vector<double>& theta);
};

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
