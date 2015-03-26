#ifndef COV_CALULATOR_H
#define COV_CALULATOR_H

//=============================================================================
// File:    covariance_calculator.h
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

#include "relpal/definitions.h"

namespace SAGE   {
namespace RELPAL {

class covariance_calculator
{
  public:

    covariance_calculator();

    ~covariance_calculator();

    bool             is_valid();
    void             compute(FPED::SubpedigreeConstPointer sp);

    const trimatrix& get_covariance_matrix() const;
    void             get_covariance_matrix(const vector<mem_pointer>& mem,
                                           trimatrix&                 c) const;
    void             get_covariance_matrix(const vector<mem_pointer>& mem,
                                           const ibd_state_info&      is_info,
                                           size_t                     mi,
                                           trimatrix&                 c) const;


  protected:

    void             compute_kinship_2p();
    void             compute_kinship_3p();
    void             compute_kinship_4p();
    void             compute_covariance();

    double           get_3p_value(int a, int b, int c);
    double           get_4p_value(int a, int b, int c, int d);

    bool             is_constant_pair(mem_pointer mid1, mem_pointer mid2) const;
    bool             is_ancestor(mem_pointer mid1, mem_pointer p) const;

    const subpedigree*            my_sped;

    trimatrix                     my_kinships_2p;
    vector<trimatrix>             my_kinships_3p;
    vector< vector<trimatrix> >   my_kinships_4p;

    trimatrix                     my_cov_matrix;

    bool                          my_valid;
};

#include "relpal/covariance_calculator.ipp"

} // end of namespace RELPAL
} // end of namespace SAGE

#endif
