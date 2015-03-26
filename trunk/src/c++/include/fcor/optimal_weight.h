#ifndef OPTIMAL_WEIGHT_H
#define OPTIMAL_WEIGHT_H

//****************************************************************************
//* File:      optimal_weight.h                                              *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Mar 04 *
//*                                                                          *
//* Notes:     This header file defines class for finding the optimal weight.*
//*                                                                          *
//* Copyright (c) 2004 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/definitions.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     optimal_weight_finder                                        ~
// ~                                                                         ~
// ~ Purpose:   Finds the optimal weights.                                   ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class optimal_weight_finder
{
  public:
    
    optimal_weight_finder();

    void estimate_variance_by_quadradic(pairset_result& result) const;    

    void estimate_weighted_covariance(size_t t,
                                      const pairset_result&             result1,
                                      const pairset_result&             result2,
                                      const Matrix2D< vector<double> >& ac_vec_matrix,
                                            Matrix2D< double >&         ac_matrix) const;

    // Compute optimal weight w1 & w2.
    //
    void find_optimal_weights(const pairset_vector&        pairsets,
                              const pairset_result_vector& results);

    const weight_matrix_vector&  get_weight_matrix() const;

  protected:

    void find_optimal_weight(const pairset_by_pedigree_type&  pairset,
                             const pairset_result&            result,
                                   weight_matrix_by_pedigree& weight);

    void debug_view(const weight_matrix& m) const;

    weight_matrix_vector      my_weight_matrix;
};

#include "fcor/optimal_weight.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
