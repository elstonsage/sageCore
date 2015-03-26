#ifndef PAIR_FILTER_H
#define PAIR_FILTER_H

//****************************************************************************
//* File:      pair_filter.h                                                 *
//*                                                                          *
//* Author:    Kevin Jacobs & Yeunjoo Song                                   *
//*                                                                          *
//* History:   Initial implementation                            kbj         *
//*            Re-organization                                   yjs May. 06 *
//*                                                                          *
//* Notes:     This header file defines pair_filter class.                   *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "palbase/rel_pair.h"

namespace SAGE    {
namespace PALBASE {

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     pair_filter                                                  ~
// ~                                                                         ~
// ~ Purpose:   Defines the object to filter out the invalid relative pair.  ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class pair_filter
{
  public:

    pair_filter();
    ~pair_filter();

    bool filter_any() const;

    bool valid          (const rel_pair& sp, bool sex_info = false, double cp = 0.0) const;

    bool valid_marker   (const rel_pair& sp)                const;
    bool valid_trait    (const rel_pair& sp, double cp=0.0) const;
    bool valid_subset   (const rel_pair& sp)                const;
    bool valid_covariate(const rel_pair& sp)                const;
    bool valid_sex_info (const rel_pair& sp)                const;

    bool is_concordant_aff_pair  (const rel_pair& sp, double cp)  const;
    bool is_concordant_unaff_pair(const rel_pair& sp, double cp)  const;
    bool is_discordant_pair      (const rel_pair& sp, double cp)  const;

    void add_marker(size_t m, double tolerance = std::numeric_limits<double>::infinity() );

    void add_trait(size_t t,
                   double min = -std::numeric_limits<double>::infinity(),
                   double max =  std::numeric_limits<double>::infinity());
    void add_trait(size_t t, size_t a,
                   double min = -std::numeric_limits<double>::infinity(),
                   double max =  std::numeric_limits<double>::infinity());
    void add_trait(size_t t, const bool* a,
                   double min = -std::numeric_limits<double>::infinity(),
                   double max =  std::numeric_limits<double>::infinity());

    void add_subset(size_t t,
                    double min = -std::numeric_limits<double>::infinity(),
                    double max =  std::numeric_limits<double>::infinity());
    void add_subset(size_t t, size_t a,
                    double min = -std::numeric_limits<double>::infinity(),
                    double max =  std::numeric_limits<double>::infinity());
    void add_subset(size_t t, const bool* a,
                    double min = -std::numeric_limits<double>::infinity(),
                    double max =  std::numeric_limits<double>::infinity());

    void add_covariate(size_t t,
                       double min = -std::numeric_limits<double>::infinity(),
                       double max =  std::numeric_limits<double>::infinity());

    void add_pair_covariate(size_t t,
                            double min = -std::numeric_limits<double>::infinity(),
                            double max =  std::numeric_limits<double>::infinity());

  private:

    class filter_trait
    {
      public:

        filter_trait();

        filter_trait(size_t t,
                     double min = -std::numeric_limits<double>::infinity(),
                     double max =  std::numeric_limits<double>::infinity());

        filter_trait(size_t t, size_t a,
                     double min = -std::numeric_limits<double>::infinity(),
                     double max =  std::numeric_limits<double>::infinity());

        filter_trait(size_t t, const bool* a,
                     double min = -std::numeric_limits<double>::infinity(),
                     double max =  std::numeric_limits<double>::infinity());

        void clear();

        bool operator==(const filter_trait& m) const;
        bool operator<(const filter_trait& m) const;

        size_t         trait;
        mutable double min;
        mutable double max;
        mutable bool   affection[3];
    };

    class filter_marker
    {
      public:

        filter_marker();

        filter_marker(size_t m, double tolerance = std::numeric_limits<double>::infinity());

        void clear();

        bool operator==(const filter_marker& m) const;
        bool operator<(const filter_marker& m) const;

        size_t         marker;
        mutable double tolerance;
    };

    typedef set<filter_marker>                                marker_set;
    typedef marker_set::iterator                              marker_iterator;
    typedef pair<marker_iterator, bool>                       marker_ibtype;

    typedef set<filter_trait>                                 trait_set;
    typedef trait_set::iterator                               trait_iterator;
    typedef pair<trait_iterator, bool>                        trait_ibtype;

    marker_set             my_markers;

    trait_set              my_traits;
    trait_set              my_subsets;
    trait_set              my_covariates;
    trait_set              my_pair_covariates;
};

#include "pair_filter.ipp"

} // end of namespace PALBASE
} // end of namespace SAGE

#endif
