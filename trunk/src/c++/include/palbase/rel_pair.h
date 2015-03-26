#ifndef A_REL_PAIR_H
#define A_REL_PAIR_H

//****************************************************************************
//* File:      rel_pair.h                                                    *
//*                                                                          *
//* Author:    Kevin Jacobs                                                  *
//*                                                                          *
//* History:   Initial implementation                            kbj         *
//*                                                                          *
//* Notes:     This header file defines access class for each rel pair.      *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "palbase/relative_pairs.h"

namespace SAGE    {
namespace PALBASE {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     rel_pair                                                     ~
// ~                                                                         ~
// ~ Purpose:   Public access class to each rel pair.                        ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class rel_pair
{
  public:

    friend class relative_pairs;
    friend class sib_cluster;

    rel_pair() : pair_num(0), data(NULL) {}

    rel_pair(size_t n, relative_pairs *pd) : pair_num(n), data(pd) {}

    void set_pair(size_t n);

    bool is_fsib_pair()    const;
    bool is_hsib_pair()    const;

    bool is_mm_pair()      const;
    bool is_mf_pair()      const;
    bool is_ff_pair()      const;

    size_t pedigree_number() const;
    size_t pair_number()     const;

    const rel_pair_data &rels() const;

    FPED::PedigreeConstPointer pedigree() const;
    relative_pairs*            pair_data() const;

    size_t marker_count()                  const;
    size_t marker_find(const string& name) const;
    string marker_name(size_t i)           const;

    double avg_share(size_t m)                       const;
    double prob_share(size_t m, size_t n)            const;
    double weighted_share(size_t m, double w0 = 0.0,
                                    double w1 = 0.5,
                                    double w2 = 1.0) const;

    double prior_avg_share()                     const;
    double prior_prob_share(size_t n)            const;
    double prior_weighted_share(double w0 = 0.0,
                                double w1 = 0.5,
                                double w2 = 1.0) const;

    bool operator==(const rel_pair &rhs) const;

  protected:

    size_t                   ped;
    size_t                   pair_num;
    relative_pairs*          data;
};

//
// rel_pair_less
//
struct rel_pair_less : public std::binary_function<rel_pair, rel_pair, bool>
{
  rel_pair_less() {}
  bool operator()(rel_pair p1, rel_pair p2) const;
};

#include "rel_pair.ipp"

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
