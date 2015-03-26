#ifndef HTEST_H
#define HTEST_H

//****************************************************************************
//* File:      htest.h                                                       *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0                                                   *
//*                                                                          *
//* Notes:     This header file defines functions to calculate chi-square &  *
//*            p-value for equivalence hypothesis test of correlations.      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/accal.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     Htest                                                        ~
// ~                                                                         ~
// ~ Purpose:   Calculates chi-square & p-values of correlations equivalence ~
// ~            test.                                                        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class Htest
{
  public:

    Htest(const PairSetData& g, const CorrelationCal& c);
    ~Htest();

    void   compute_homogeneity_test(const pairset_result_vector& sub_results,
                                          Htest_result_vector&   H_results);

  protected:

    bool   compute_asymptotic_covariance(main_to_sub_type_const_iterator       rel_group,
                                         const pairset_result_vector&          sub_results,
                                         vector< vector< Matrix2D<double> > >& ac_matrix)  const;

    void   compute_linear_correlation   (main_to_sub_type_const_iterator     rel_group,
                                         const pairset_result_vector&        sub_results,
                                         Matrix2D<internal_real_type>&       y,
                                         size_t                              t)          const;

    void   compute_omega      (const vector< vector< Matrix2D<double> > >&   ac_matrix,
                                     Matrix2D<internal_real_type>&           omega,
                               size_t                                        trait_cnt)  const;

    void   compute_omega_intra(const Matrix2D<internal_real_type>&           old_omega,
                                     Matrix2D<internal_real_type>&           omega,
                               size_t                                        df,
                               size_t                                        trait_cnt)  const;

    double compute_chi_square (const Matrix2D<internal_real_type>&           y,
                               const Matrix2D<internal_real_type>&           omega)      const;

    //void   debug_view(const Matrix2D<internal_real_type>& m) const;
    //void   debug_view(const Matrix2D<double>&             m) const;

    const PairSetData*            my_pairsetdata;
    const CorrelationCal*         my_correlation;

    bool                          my_conservative;
    bool                          my_individual_homog;
};

//#include "fcor/htest.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
