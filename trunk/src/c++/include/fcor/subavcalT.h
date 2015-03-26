#ifndef SUBAVCALT_H
#define SUBAVCALT_H

//****************************************************************************
//* File:      subavcalT.h                                                   *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Oct 99 *
//*                                                                          *
//* Notes:     This header file calculates the asymptotic variance of        *
//*            a multipedigree correlation(XY & XY).                         *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/corcal.h"

namespace SAGE {
namespace FCOR {
    
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubAVCalT                                                     ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic variance of a correlation(XY & XY).     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SubAVCalT
{
  public:
    
    // New pair types created from given two XY & XY pair types. 
    enum pair_type { XX=0, XY, YX, YY };

    SubAVCalT(const SubPairSetData* p, const SubCorrelationCal* c);

    void  compute_asymptotic_variance(const pairset_by_pedigree_type& pairset,
                                      const CorrelationInfo&          xy,
                                      TriangleMatrix<double>&         av,
                                      TriangleMatrix<size_t>&         replaced);

  private:
    
    void    init();
    void    compute(TriangleMatrix<double>& av);

    void    compute_we();

    void    build_weight_triangle(TriangleMatrix<double>& weight, size_t i, size_t p);

    void    compute_sub(const TriangleMatrix<internal_real_type>& sum_p,
                      vector< TriangleMatrix<internal_real_type> >& ikj,
                              TriangleMatrix<internal_real_type>&   kl,
                        size_t p);

    void    module1_8  (const TriangleMatrix<internal_real_type>& s1,
                        const TriangleMatrix<internal_real_type>& s2,
                              TriangleMatrix<internal_real_type>& result,
                        size_t i, size_t j, size_t p)                             const;

    void    module9    (const TriangleMatrix<internal_real_type>& s_xx,
                        const TriangleMatrix<internal_real_type>& s_yy,
                        const TriangleMatrix<internal_real_type>& s_xy,
                        const TriangleMatrix<internal_real_type>& s_yx,
                              TriangleMatrix<internal_real_type>& result,
                        size_t i, size_t j, size_t p)                             const;
    
    void    sub_module (      TriangleMatrix<internal_real_type>& sub,
                        size_t i, size_t j, size_t p, pair_type p_t)              const;
    
    void    pedigree_weight(TriangleMatrix<internal_real_type>& result, size_t p) const;

    size_t  pedigree_count()                                                      const;
    size_t  trait_count()                                                         const;
    size_t  xy_pair_count(size_t p)                                               const;

    double  get_corinfo(size_t i, size_t j, size_t p, size_t t1, size_t t2)       const;

    const SubPairSetData*                                       my_pairsetdata;
    const SubCorrelationCal*                                    my_correlation;
    const CorrelationInfo*                                      my_corinfo_xy;

    // Pointer to a pairset.
    const pairset_by_pedigree_type*                             my_pairset_xy;
    
    // Conservative or robust.
    bool                                                        my_conservative;
    
    // Weight triangle for each member pair.
    vector< vector< TriangleMatrix< double > > >              	my_weight_w;

    // Common weight factor used in modules of calculation formular of asymptotic covariance. 
    TriangleMatrix< internal_real_type >	                my_we;
      
    // Common sum of correlations used throughout the calculating process.
    // i, j, k, l indicate loop index for each summation.
    vector< TriangleMatrix< internal_real_type > >              my_sum_kl;
    vector< vector< TriangleMatrix< internal_real_type > > >    my_sum_ikj;
};

#include "fcor/subavcalT.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
