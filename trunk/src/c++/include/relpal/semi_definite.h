#ifndef  RELPAL_CONSTRAINTS_H
#define  RELPAL_CONSTRAINTS_H

//==========================================================================
//  File:       semi_definite.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial Implementation.                              Sep. 07
//
//  Notes:      This file contains two methods to constrain the matrix to be
//              positive semi-definite.
//
//  Copyright (c) 2007 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "maxfun/maxfun.h"
#include "relpal/definitions.h"

namespace SAGE
{
namespace RELPAL
{

class semi_definite : public MaxFunction
{
  public :

    semi_definite(const matrix& U, const matrix& si, size_t t);
    
    void           compute_eigenvalue();
    void           initialize();

    const matrix&  get_estimates() const;
    bool           is_invalid() const;

  protected :

    virtual double evaluate     (vector<double>& theta);
    virtual int    update_bounds(vector<double>& theta);

    virtual double evaluate_derived     (vector<double>& theta) = 0;
    virtual int    update_bounds_derived(vector<double>& theta) = 0;


    const matrix&  my_vector_U;
    const matrix&  my_sigma_i;
    size_t         my_trait_count;

    matrix         my_current_estimate;  // w for cholesky_decom, b for cutting_plane, b = ww'

    matrix         my_initial_b;
    vector<double> my_current_e_value;
    matrix         my_current_e_vector;
    int            my_current_info;

    bool           my_invalid;
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Class:     cholesky_decom
//
//  Notes:     This class defines the Cholesky decomposition method to
//             constrain the matrix to be positive semi-definite.
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class cholesky_decom : public semi_definite
{
  public :

    cholesky_decom(const matrix& U, const matrix& si, size_t t);
    
    void           initialize();

  protected :

    virtual double evaluate_derived     (vector<double>& theta);
    virtual int    update_bounds_derived(vector<double>& theta);

};

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//  Class:     cutting_plane
//
//  Notes:     This class defines the Cutting Plane Algorithm to find
//             feasible estimates (positive semidefinite) in covariance
//             component analysis.
//             Shaw & Geyer, Biometrika(1997), 84, 95-102
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class cutting_plane : public semi_definite
{
  public :

    cutting_plane(const matrix& U, const matrix& si, size_t t);
    
    bool           is_positive_semidefinite() const;

    void           reset_b();
    void           add_constraints();

  protected :

    virtual double evaluate_derived     (vector<double>& theta);
    virtual int    update_bounds_derived(vector<double>& theta);

    vector<matrix> my_v;
};

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
