#ifndef CORCAL_H
#define CORCAL_H

//****************************************************************************
//* File:      corcal.h                                                      *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   0.0 Initial implementation                         yjs Oct 99 *
//*            1.0 Revised to implement MainCorCal                yjs May 01 *
//*                                                                          *
//* Notes:     This header file computes  correlations for each type of      *
//*            pairset in PairSetData.                                       *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "fcor/pairsetdata.h"

namespace SAGE {
namespace FCOR {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     CorrelationCal                                               ~
// ~                                                                         ~
// ~ Purpose:   Calculate correlations of PairSetData.                       ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class CorrelationCal
{
  public:

    CorrelationCal();
    CorrelationCal(const FcorParser& p);

    void  set_parser(const FcorParser& p);

    void  compute_correlations(const pairset_vector&        pairset,
                               const pairset_info_vector&   pinfo,
                                     pairset_result_vector& result);

    void  compute_correlations(const pairset_vector&        pairset,
                               const pairset_info_vector&   pinfo,
                               const weight_matrix_vector&  weights,
                                     pairset_result_vector& result);

    bool                           is_valid()             const;

    const FcorParser*              get_parser()           const;

    const pairset_vector*          get_pairset()          const;
    const pairset_info_vector*     get_pairset_info()     const;

    const corinfo_vector&          get_corinfo()          const;
    const pair_to_corinfo_map&     get_corinfo_map()      const;

    void  view_corinfo(ostream& out, bool see_class=true) const;

  private:

    void  compute_correlations(pairset_result_vector& result);
    void  compute_correlations(const weight_matrix_vector&  weights,
                                     pairset_result_vector& result);

    void  set_correlation(pairset_result& pr, const vector<CorrelationInfo>& cor);

    void  build_pair_to_corinfo_map();
    
    void  compute_pairset_correlation        (const pairset_by_pedigree_type& pair,
                                                    corinfo_by_weight_type&   cor) const;


    //void  compute_pairset_optimal_correlation(const pairset_by_pedigree_type&  pair,
    //                                          vector<CorrelationInfo>&         cor,
    //                                          const weight_matrix_by_pedigree& w)        const;

    const FcorParser*          my_parser;
    const pairset_vector*      my_pairset;
    const pairset_info_vector* my_pairset_info;

    corinfo_vector             my_corinfo;
    pair_to_corinfo_map        my_corinfo_map;

    bool                       my_valid;
};

#include "fcor/corcal.ipp"    

} // end of namespace FCOR
} // end of namespace SAGE

#endif




