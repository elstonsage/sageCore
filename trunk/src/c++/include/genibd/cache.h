#ifndef GENIBD_CACHE_H
#define GENIBD_CACHE_H

//============================================================================
// File:      cache.h
//
// Author:    Dan Baechle
//
// History:   10/14/2 - created.                                 djb
//            Modified for GENIBD                                yjs Apr 2004
//
// Notes:     defines a peeling cache class for calculating subpedigree 
//            likelihood.
//
// Copyright (c) 2002 R.C. Elston
// All Rights Reserved
//============================================================================

#include "peeling/cache3.h"
#include "numerics/log_double.h"
#include "genibd/definitions.h"

using std::vector;
using std::pair;

namespace SAGE
{

namespace peeling
{

typedef pair<size_t, vector<log_double> >  posterior_vector;

template<class T1, class T2> 
struct pv_equal : public std::binary_function<T1, T2, bool>
{
  bool operator()(T1 first, T2 second) const;
};

//----------------------------------------------------------------------------
//  Class:    individual_cache<ind_pen_iter, log_double>
//
//  Purpose:  specialization of a class to store information about an
//            individual's anterior and posterior values for use when calc-
//            ulating pedigree likelihoods by the method of Fernando, Stricker
//            and Elston.  1993.
//
//----------------------------------------------------------------------------
//
template<>
class individual_cache<MLOCUS::penetrance_model::phased_penetrance_iterator,
                       log_double>
{
  public:

    typedef log_double                                           result_type;
    typedef MLOCUS::penetrance_model::phased_penetrance_iterator data_type;
      
    typedef FPED::Member member_type;
    
    individual_cache();
    void build(const member_type& ind, size_t mph_id, const MLOCUS::inheritance_model& mm);
  
    bool  anterior_cached(const data_type& upi) const;
    bool  posterior_cached(const data_type& upi) const;
    bool  posterior_with_mate_cached(const member_type& mate_index, const data_type& upi) const;
    bool  posterior_except_mate_cached(const member_type& mate_index, const data_type& upi) const;
    
    log_double&  anterior(const data_type& upi);
    log_double&  posterior(const data_type& upi);
    log_double&  posterior_with_mate(const member_type& mate_index, const data_type& upi);
    log_double&  posterior_except_mate(const member_type& mate_index, const data_type& upi);
    
    const log_double&  anterior(const data_type& upi) const;
    const log_double&  posterior(const data_type& upi) const;
    const log_double&  posterior_with_mate(const member_type& mate_index, const data_type& upi) const;
    const log_double&  posterior_except_mate(const member_type& mate_index, const data_type& upi) const;

  private:

    size_t  genotype_index(const data_type& upi) const;

    vector<log_double>&  posteriors_with_mate(size_t mate) const;
    vector<log_double>&  posteriors_except_mate(size_t mate) const;
    
    // Data members.
    //
    std::map<size_t, size_t>               m_indices;  // <marker genotype id, matrix col>
    
    mutable vector<log_double>        my_anteriors;
    mutable vector<log_double>        my_posteriors;
    mutable vector<posterior_vector>  my_posteriors_with_mate;
    mutable vector<posterior_vector>  my_posteriors_except_mate;
};

#include "genibd/cache.ipp"

} // end of GENIBD namespace
} // end of SAGE namespace

#endif
