#ifndef RELPAL_CALULATOR_H
#define RELPAL_CALULATOR_H

//=============================================================================
// File:    two_level_calculator.h
//
// Author:  Yeunjoo Song
//
// History: Version 0.0 Initial implementation
//
// Notes:   This file contains definition for following data structure.
//            class calculator_base
//            class relpal_least_square
//            class relpal_score
//
// Copyright (c) 2007   R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "relpal/parser.h"

namespace SAGE   {
namespace RELPAL {

class calculator_base
{
  public:

    calculator_base(size_t m = 0);

    virtual ~calculator_base();

    virtual void reset() = 0;
    virtual void reset(size_t m) = 0;

    virtual void add_block     (const matrix& y, const trimatrix& W, const matrix& A, trimatrix& Wi) = 0;
    virtual void add_block_kron(const matrix& y, const trimatrix& b, const matrix& A) = 0;

    virtual bool compute() = 0;

    void     set_debug_out(ostream* o, bool R_format=false);

    matrix&  triXZ (const trimatrix& X, const matrix& Z, matrix& o) const;
    matrix&  XtriZ (const matrix& X, const trimatrix& Z, matrix& o) const;
    matrix&  XTtriZ(const matrix& X, const trimatrix& Z, matrix& o) const;

    void     get_tri_kron(const trimatrix& b, trimatrix& wi) const;

    size_t   parameters;
    size_t   observation_count;
    size_t   cluster_count;
    size_t   svd_return_code;

    // Decompositions
    SAGE::SVD  svd;
    SAGE::LUP  lup;
    
  protected:

    void touch_up_matrix(matrix& m) const;

    int  get_tri_inverse(const trimatrix& w, trimatrix& wi) const;

    ostream*  my_debug_out;
};

class relpal_least_square : public calculator_base
{
  public:

    relpal_least_square(size_t m = 0);

    ~relpal_least_square();

    virtual void reset();
    virtual void reset(size_t m);

    virtual void add_block     (const matrix& y, const trimatrix& W, const matrix& A, trimatrix& Wi);
    virtual void add_block_kron(const matrix& y, const trimatrix& b, const matrix& A);

    virtual bool compute();

    void build_residuals(const matrix& y, const matrix& A, matrix& r) const;

    matrix   beta;
    matrix   information;
    matrix   Variance;
    matrix   residual_variance;
    matrix   total_variance;
    
  protected:

    matrix   BA;                // partial B A
    matrix   By;                // partial B y
    matrix   yWy;               // weighted sum of squares
    matrix   AWA;               // sum(A' W A)
    matrix   AWy;               // sum(A' W y)
    matrix   temp1,temp2,temp3;
};


class relpal_score : public calculator_base
{
  public:

    relpal_score(size_t m = 0);

    ~relpal_score();

    virtual void  reset();
    virtual void  reset(size_t m);

    virtual void  add_block     (const matrix& y, const trimatrix& W, const matrix& A, trimatrix& Wi);
    virtual void  add_block_kron(const matrix& y, const trimatrix& b, const matrix& A);

    virtual bool  compute();

    void     add_block_IBD(const matrix& B, const trimatrix& varIBD);

    bool     is_naive_var()       const;
    bool     is_sandwich_var()    const;
    bool     is_alternative_var() const;
    bool     is_IBD_var()         const;

    void     set_naive_var      (bool v = true);  
    void     set_sandwich_var   (bool v = true);  
    void     set_alternative_var(bool v = true);  
    void     set_IBD_var        (bool v = true);  

    matrix   U_star;            // (sum(D' W^ D))^U

    matrix   sigma_na;          // (sum(D' W^ D))^
    matrix   sigma_na_i;        // sigma_na*^-1
    matrix   sigma_na_sqrt;     // sigma_na*^0.5

    matrix   sigma_sw;          // sum((D' W^ D)^(sigma)(D' W^ D)^)
    matrix   sigma_sw_i;        // sigma_sw*^-1
    matrix   sigma_sw_sqrt;     // sigma_sw*^0.5

    matrix   sigma_al;          // sum(D' W^ D)^sum((U-EU)(U-EU)')sum(D' W^ D)^
    matrix   sigma_al_i;        // sigma_al*^-1
    matrix   sigma_al_sqrt;     // sigma_al*^0.5

    matrix   sigma_ib;          // sum(D' W^ D)^sum(B Var(phi) B')sum(D' W^ D)^
    matrix   sigma_ib_i;        // sigma_ib*^-1
    matrix   sigma_ib_sqrt;     // sigma_ib*^0.5

    matrix   T_na;              // U*'sigma_na^U*
    matrix   T_sw;              // U*'sigma_sw^U*
    matrix   T_al;              // U*'sigma_al^U*
    matrix   T_ib;              // U*'sigma_ib^U*

    size_t   df_na;
    size_t   df_sw;
    size_t   df_al;
    size_t   df_ib;

  protected:

    void     compute_U_sigma(const matrix& D, const trimatrix& Wi, const matrix& S);

    matrix   U;                 // sum(D' W^ S)
    matrix   sigma;             // sum((D' W^ S)(S' W^ D))

    matrix   DWD;               // sum(D' W^ D)
    matrix   BPB;               // sum(B var(phi) B')

    bool     my_naive_var;
    bool     my_sandwich_var;
    bool     my_alternative_var;
    bool     my_IBD_var;

    vector<matrix> Us;          // Place to save U = DWS foe each pedigree for sigma_alt
    vector<matrix> DWDs;        // Place to save DWD for each pedigree for sigma_alt
};

#include "relpal/two_level_calculator.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
