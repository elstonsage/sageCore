#ifndef ACCAL_H
#define ACCAL_H

//****************************************************************************
//* File:      accal.h                                                       *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Oct 99 *
//*            1.0 Maintype computation added                     yjs May 01 *
//*            1.1 Conservative & robust modified                 yjs Sep 01 *
//*                                                                          *
//* Notes:     This header file defines class for calculating the asymptotic *
//*            variance/covariances of familiar correlations.                *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/subcalM.h"
#include "fcor/optimal_weight.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     ACCal                                                        ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic var/covariances of correlations.        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class ACCal
{
  public:
    
    ACCal(const CorrelationCal* c);

    void  compute_covariance(const pairset_by_pedigree_type&  rel_type1,
                             const pairset_by_pedigree_type&  rel_type2,
                             const vector<CorrelationInfo>&   cor1,
                             const vector<CorrelationInfo>&   cor2,
                             const weight_matrix_by_pedigree& optimal_w1,
                             const weight_matrix_by_pedigree& optimal_w2,
                             Matrix2D< vector<double> >&      ac_matrix,
                             Matrix2D< vector<bool> >&        replaced_matrix,
                             Matrix2D< vector<size_t> >&      least_pair_matrix) const;

    void  compute_variance  (const pairset_by_pedigree_type&  rel_type1,
                             const vector<CorrelationInfo>&   cor1,
                             const weight_matrix_by_pedigree& optimal_w1,
                             Matrix2D< vector<double> >&      ac_matrix,
                             Matrix2D< vector<bool> >&        replaced_matrix,
                             Matrix2D< vector<size_t> >&      least_pair_matrix) const;

  protected:

    const CorrelationCal*    my_correlation;
};


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     StdErrCal                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic standard errors of correlations.        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class StdErrCal
{
  public:
    
    StdErrCal(const CorrelationCal& c);

    void  compute_standard_errors(pairset_result_vector& results);
    void  compute_standard_errors(const weight_matrix_vector&  weights,
                                        pairset_result_vector& results);

  protected:

    void  set_standard_error(pairset_result& pr, const Matrix2D< vector<double> >  ac_matrix,
                                                 const Matrix2D< vector<bool> >    re_matrix,
                                                 const Matrix2D< vector<size_t> >  lp_matrix) const;

    void  compute_pooled_cross_correlation(pairset_result&           pr,
                                           const Matrix2D< double >& ac_matrix) const;

    const CorrelationCal*    my_correlation;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     VarCovCal                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate variance-covariance matrix of correlations.        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class VarCovCal
{
  public:

    VarCovCal(const CorrelationCal& sc, const CorrelationCal& mc);
    
    void   compute_variance_covariances(const pairset_result_vector& sub_results,
                                        const pairset_result_vector& main_results,
                                              var_cov_result_vector& vc_results);

  protected:

    void   compute_variance_covariance(size_t i, size_t p1, size_t p2,
                                       const pairset_result_vector& results,
                                       const CorrelationCal*        corcal);

    size_t find_reltype(const pairset_info_vector& pinfos, string reltype) const;

    const CorrelationCal*       my_sub_correlation;
    const CorrelationCal*       my_main_correlation;

    var_cov_result_vector       my_var_covs;
};

#include "fcor/accal.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
