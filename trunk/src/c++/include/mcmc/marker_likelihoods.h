#ifndef MARKER_LIKELIHOODS_H
#define MARKER_LIKELIHOODS_H

//==========================================================================
//  File:    marker_likelihoods.h
//
//  Author:  Geoff Wedig
//
//  History: Version 0.90
//           1.0 Updated to new libraries                        yjs May. 04
//
//  Notes:
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include <math.h>
#include <stdio.h>
#include "containers/cache_map.h"
#include "mcmc/hash.h"
#include "mcmc/mcmc_data_accessor.h"
#include "mcmc/founder_allele_graph.h"

namespace SAGE
{

namespace MCMC
{

class marker_likelihood_calculator
{
  public:

    marker_likelihood_calculator(const pedigree_region&  pr,
                                 const McmcMeiosisMap& ped,
                                 mcmc_data_accessor&     dat);
    
    ~marker_likelihood_calculator();
    
    bool   valid_likelihood()                         const;  // Are all markers valid?
    bool   valid_likelihood(size_t m)                 const;  // Is marker m valid?

    double log_likelihood()                           const;  // Returns total likelihood (qnan if invalid)
    double log_likelihood(size_t m)                   const;  // Returns likelihood of a given marker (or qnan)

    double log_likelihood(size_t m, const bit_field&) const;  // Returns likelihood of bits 
                                                              // at given marker (or qnan)
    void   dump_graphs(ostream& o) const;

  protected:

    typedef vector<FounderAlleleGraph>                    graph_vector;
    typedef SAGE::cache_map<bit_field, log_double, hash1> cache_type;
    typedef cache_type::const_iterator                    hash_iterator;

    struct CachingData
    {
      public:
      
        CachingData();
        
        cache_type cache;
        size_t     hit_count;
    };

    typedef vector<CachingData> CacheDataVector;
    
    mutable CacheDataVector my_cache_data;
    mutable graph_vector    my_graphs;
    
    /// Basic Data
    //@{
    
    const pedigree_region&  my_region;
    mcmc_data_accessor*     my_data;

    //@}
};


#include "mcmc/marker_likelihoods.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
