#ifndef LODPAL_RESULT_H
#define LODPAL_RESULT_H

//****************************************************************************
//* File:      lodpal_result.h                                               *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                        yjs Nov. 06 *
//*                                                                          *
//* Notes:     This header file defines analysis result class for lodpal.    *
//*                                                                          *
//* Copyright (c) 2006 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "lodpal/lodpal_pairs.h"

namespace SAGE   {
namespace LODPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     lodpal_result                                                ~
// ~                                                                         ~
// ~ Purpose:   Defines result storage class for lodpal test                 ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class lodpal_result
{
  public:

    lodpal_result() {}

    void  clear();

    void  set_maxfun_result(const maxfun_results& re);
    void  set_method(char c);                           

    void  set_lod_score(const pair<double, double>& ld);

    void  set_lod_score_mm(double mm);
    void  set_lod_score_mf(double mm);
    void  set_lod_score_ff(double mm);

    double                   lod_score()               const;
    int                      function_evaluations()    const;
    int                      last_error()              const;
    int                      ivfl()                    const;
    int                      igage()                   const;
    char                     method()                  const;
    const Matrix2D<double>&  var_cov_matrix()          const;

    double                   lod_score_mm()            const;
    double                   lod_score_mf()            const;
    double                   lod_score_ff()            const;

  protected:

    maxfun_results            my_maxfun_result;
    char                      my_method;
    Matrix2D<double>          my_vc_matrix;

    double                    my_lod_score_with_cap;   
    double                    my_lod_score_without_cap;

    double                    my_lod_score_mm;
    double                    my_lod_score_mf;
    double                    my_lod_score_ff;
};

#include "lodpal/lodpal_result.ipp"

} // end of namespace LODPAL
} // end of namespace SAGE

#endif
