#ifndef SIB_CLUSTER_H
#define SIB_CLUSTER_H

//=============================================================================
// File:    sib_cluster.h
//
// Author:  Yeunjoo Song
//
// History: Initial implementation                                 May. 06
//
// Notes:   This file contains definition for object to hold sib_cluster.
//
// Copyright (c) 2006 R.C. Elston
//   All Rights Reserved
//=============================================================================

#include "sibpal/definitions.h"

namespace SAGE   {
namespace SIBPAL {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     sib_cluster                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines the object to store sib_cluster.                     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class sib_cluster
{
  public:

    sib_cluster();

    sib_cluster(const vector<size_t>& fs, const vector<size_t>& hs,
                relative_pairs& rp, bool fs, bool hs,
                const pair_filter& p, bool x_linked = false, bool mm = true, bool mf = true, bool ff = true);

    size_t  valid_sib_count() const;

    // Number of sib pairs
    size_t valid_pair_count()      const;
    size_t valid_fsib_pair_count() const;
    size_t valid_hsib_pair_count() const;

    size_t valid_fsib_mm_pair_count() const;
    size_t valid_fsib_mf_pair_count() const;
    size_t valid_fsib_ff_pair_count() const;

    size_t valid_hsib_mm_pair_count() const;
    size_t valid_hsib_mf_pair_count() const;
    size_t valid_hsib_ff_pair_count() const;

    sib_pair operator[](size_t i) const;

    bool operator==(const sib_cluster &rhs) const;
    bool operator!=(const sib_cluster &rhs) const;

    const sib_set& get_valid_sibs() const;

    void dump() const;

  protected:

    const pair_filter* my_filter;

    relative_pairs*    my_data;

    vector<size_t>     my_valid_fsib_pairs;
    vector<size_t>     my_valid_hsib_pairs;

    vector<size_t>     my_valid_fsib_mm_pairs;
    vector<size_t>     my_valid_fsib_mf_pairs;
    vector<size_t>     my_valid_fsib_ff_pairs;

    vector<size_t>     my_valid_hsib_mm_pairs;
    vector<size_t>     my_valid_hsib_mf_pairs;
    vector<size_t>     my_valid_hsib_ff_pairs;

    sib_set            my_sibs;
};

#include "sib_cluster.ipp"

} //end of namespace SIBPAL
} //end of namespace SAGE

#endif
