#ifndef SUBCALM_H
#define SUBCALM_H

//****************************************************************************
//* File:      subcalM.h                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 0.0 Initial implementation                 yjs Oct 99 *
//*            1.0 Reconstructing the inheritance                 yjs 2000   *
//*            1.1 Conservative & robust resolved                 yjs Sep 01 *
//*                                                                          *
//* Notes:     This header file defines the classes to calculates the        *
//*            asymptotic variance & covariance of a multipedigree           *
//*            correlation.                                                  *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/corcal.h"

namespace SAGE {
namespace FCOR {
    
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubCalBaseM                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines the base class for SubACCalM & SubAVCalM.            ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SubCalBaseM
{
  public:

    // New pair types created from given two XY & XY pair types and two UV & XY pair types.
    enum pair_type { XX=0, XY, YX, YY, UX, UY, VX, VY };

    SubCalBaseM();
    SubCalBaseM(const CorrelationCal* cc);
    ~SubCalBaseM();

  protected:
    
    void    compute_te_we(Matrix2D< vector<internal_real_type> >& wt,
                          const vector<CorrelationInfo>*          type);

    void    build_weight_matrix(Matrix2D< vector<double> >&     weight,  const size_t& i,
                                const pairset_by_pedigree_type* pairset, const size_t& ped,
                                const weight_matrix&            opt_weight );

    void    pedigree_weight(const vector< Matrix2D< vector<double> > >& weight,
                            Matrix2D< vector<internal_real_type> >&     sum_p)   const;

    size_t   pedigree_count()                                                    const;
    size_t   trait_count()                                                       const;
    size_t   pair_count(const pairset_by_pedigree_type* pairset, size_t p)       const;

    vector< pair<double, size_t> > get_corinfo(const pairset_by_pedigree_type* pairset1,
                                               const pairset_by_pedigree_type* pairset2,
                                               size_t i, size_t j, size_t p, pair_type p_t,
                                               size_t t1, size_t t2)             const;

    const CorrelationCal*                                       my_correlation;

    const vector<CorrelationInfo>*                              my_corinfo_xy;

    // Pointer to pairset.
    const pairset_by_pedigree_type*                             my_pairset_xy;

    // Conservative or robust.                                                  
    bool                                                        my_conservative;

    // Weight triangle for each member pair.
    vector< vector< Matrix2D< vector<double> > > >              my_weight_w;

    // Common weight factor used in modules of calculation formular of asymptotic covariance. 
    Matrix2D< vector<internal_real_type> > 	    	        my_we;

    Matrix2D< vector<size_t> >                                  my_least_pair_count;
    Matrix2D< vector<size_t> >                                  my_extended_least_pair_count;

    size_t                                                      my_pedigree_count;
    size_t                                                      my_trait_count;
    size_t                                                      my_weight_count;
};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubAVCalM                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic variance of a correlation(XY & XY).     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SubAVCalM : SubCalBaseM
{
  public:
    
    SubAVCalM(const CorrelationCal* cc);
    ~SubAVCalM();

    void  compute_asymptotic_variance(const pairset_by_pedigree_type&  pairset,
                                      const vector<CorrelationInfo>&   xy,
                                      const weight_matrix_by_pedigree& optimal_w,
                                      Matrix2D< vector<double> >&      ac,
                                      Matrix2D< vector<bool> >&        replaced,
                                      Matrix2D< vector<size_t> >&      least_pair);
                                   
  private:
    
    // Common sum of correlations used throughout the calculating process.
    // i, j, k, l indicate loop index for each summation.
    // xx, xy, yx, yy indicate new pairsets from existing xy & xy pairsets.  
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_xx;
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_xy;
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_yy;
    
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ikj_x;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ik_xy;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_kj_xy;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ikj_y;

    void    init(const weight_matrix_by_pedigree& optimal_w);
    void    compute(Matrix2D< vector<double> >& ac);

    void    compute_sub(const Matrix2D< vector<internal_real_type> >&   sum_p,
                      vector< Matrix2D< vector<internal_real_type> > >& ikj,
                              Matrix2D< vector<internal_real_type> >&    kl,
                        size_t p, pair_type p_t);

    void    module1_8  (const Matrix2D< vector<internal_real_type> >& s1,
                        const Matrix2D< vector<internal_real_type> >& s2,
                              Matrix2D< vector<internal_real_type> >& result,
                        size_t i, size_t j, size_t p) const;

    void    module9    (const Matrix2D< vector<internal_real_type> >& s_xx,
                        const Matrix2D< vector<internal_real_type> >& s_yy,
                        const Matrix2D< vector<internal_real_type> >& s_xy,
                        const Matrix2D< vector<internal_real_type> >& s_yx,
                              Matrix2D< vector<internal_real_type> >& result,
                        size_t i, size_t j, size_t p) const;

    void    sub_module (      Matrix2D< vector<internal_real_type> >& sub,
                        size_t i, size_t j, size_t p, pair_type p_t);

    void    sub_xx_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_xy_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_yx_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_yy_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

};

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     SubACCalM                                                    ~
// ~                                                                         ~
// ~ Purpose:   Calculate asymptotic covariance of two correlations(XY & UV).~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class SubACCalM : SubCalBaseM
{
  public:
    
    SubACCalM(const CorrelationCal* cc);
    ~SubACCalM();

    void  compute_asymptotic_covariance(const pairset_by_pedigree_type& p1,
                                        const pairset_by_pedigree_type& p2,
                                        const vector<CorrelationInfo>&  t1,
                                        const vector<CorrelationInfo>&  t2,
                                        const weight_matrix_by_pedigree& optimal_w1,
                                        const weight_matrix_by_pedigree& optimal_w2,
                                        Matrix2D< vector<double> >&     ac,
                                        Matrix2D< vector<bool> >&       replaced,
                                        Matrix2D< vector<size_t> >&     least_pair);
                                   
  private:
    
    const vector<CorrelationInfo>*                              my_corinfo_uv;

    // Pointer to two pairsets.
    const pairset_by_pedigree_type*                             my_pairset_uv;

    // Weight triangle for each member pair.
    vector< vector< Matrix2D< vector<double> > > >              my_weight_t;

    // Common weight factor used in modules of calculation formular of asymptotic covariance. 
    Matrix2D< vector<internal_real_type> >	                my_te;
      
    // Common sum of correlations used throughout the calculating process.
    // i, j, k, l indicate loop index for each summation.
    // ux, uy, vx, vy indicate new pairsets from existing xy & uv pairsets.  
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_ux;
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_uy;
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_vx;
    vector< Matrix2D< vector<internal_real_type> > >                    my_sum_kl_vy;
    
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ik_ux;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ik_uy;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ik_vx;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_ik_vy;

    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_kj_ux;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_kj_uy;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_kj_vx;
    vector< vector< Matrix2D< vector<internal_real_type> > > >          my_sum_kj_vy;
    
    void    init(const weight_matrix_by_pedigree& optimal_w, const weight_matrix_by_pedigree& w2);
    void    compute(Matrix2D< vector<double> >& ac);


    void    compute_ik (const Matrix2D< vector<internal_real_type> >& sum_pw,
                              Matrix2D< vector<internal_real_type> >& ik,
                        size_t i, size_t p, pair_type p_t);

    void    compute_kj (const Matrix2D< vector<internal_real_type> >& sum_pt,
                              Matrix2D< vector<internal_real_type> >& kj,
                        size_t j, size_t p, pair_type p_t);
                                                   
    void    compute_kl (const Matrix2D< vector<internal_real_type> >& sum_pt,
                        const Matrix2D< vector<internal_real_type> >& sum_pw,
                              Matrix2D< vector<internal_real_type> >& kl,
                        size_t p, pair_type p_t);

    void    module1_8  (const Matrix2D< vector<internal_real_type> >& s1,
                        const Matrix2D< vector<internal_real_type> >& s2,
                              Matrix2D< vector<internal_real_type> >& result,
                        size_t i, size_t j, size_t p) const;

    void    module9    (const Matrix2D< vector<internal_real_type> >& s_ux,
                        const Matrix2D< vector<internal_real_type> >& s_vy,
                        const Matrix2D< vector<internal_real_type> >& s_uy,
                        const Matrix2D< vector<internal_real_type> >& s_vx,
                              Matrix2D< vector<internal_real_type> >& result,
                        size_t i, size_t j, size_t p) const;

    void    sub_module (      Matrix2D< vector<internal_real_type> >& sub,
                        size_t i, size_t j, size_t p, pair_type p_t);

    void    sub_ux_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_uy_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_vx_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);

    void    sub_vy_type(      Matrix2D< vector<internal_real_type> >& sub, 
                        const Matrix2D< vector<double> >&             cor, 
                        const Matrix2D< vector<internal_real_type> >* ik,  
                        const Matrix2D< vector<internal_real_type> >* kj,  
                        const Matrix2D< vector<internal_real_type> >* kl);
};

#include "fcor/subcalM.ipp"

} // end of namespace FCOR
} // end of namespace SAGE

#endif
