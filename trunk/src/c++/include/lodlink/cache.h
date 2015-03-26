#ifndef LODLINK_CACHE_H
#define LODLINK_CACHE_H
//============================================================================
// File:      cache.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   10/14/2 - created.                                   djb
//                                                                          
// Notes:     defines a peeling cache class for calculating subpedigree 
//            likelihood.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>
#include <map>
#include <utility>
#include <functional>
#include <algorithm>
#include "fped/fped.h"
#include "peeling/cache3.h"
#include "peeling/peeler3.h"
#include "numerics/log_double.h"
#include "mlocus/penmodel.h"
#include "lodlink/definitions.h"
#include "lodlink/trans_calculator.h"

using std::vector;
using std::pair;
using SAGE::LODLINK::joint_pen_iter;

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
//  Class:    individual_cache<joint_pen_iter, log_double>
//                                                                          
//  Purpose:  specialization of a class to store information about an
//            individual's anterior and posterior values for use when calc-
//            ulating pedigree likelihoods by the method of Fernando, Stricker
//            and Elston.  1993.
//                                                                          
//----------------------------------------------------------------------------
//
template<>
class individual_cache<joint_pen_iter, log_double>
{
  public:
    typedef log_double  result_type;
    typedef joint_pen_iter  data_type;
      
    typedef FPED::FilteredMultipedigree::member_type  member_type;
    
    individual_cache();
    void build(const member_type& ind, size_t tph_id, size_t mph_id, 
               const MLOCUS::penetrance_model& tm, const MLOCUS::penetrance_model& mm);
  
    bool  anterior_cached(const joint_pen_iter& jpi) const;
    bool  posterior_cached(const joint_pen_iter& jpi) const;
    bool  posterior_with_mate_cached(const member_type& mate_index, const joint_pen_iter& jpi) const;
    bool  posterior_except_mate_cached(const member_type& mate_index, const joint_pen_iter& jpi) const;
    
    log_double&  anterior(const joint_pen_iter& jpi);
    log_double&  posterior(const joint_pen_iter& jpi);
    log_double&  posterior_with_mate(const member_type& mate_index, const joint_pen_iter& jpi);
    log_double&  posterior_except_mate(const member_type& mate_index, const joint_pen_iter& jpi);
    
    const log_double&  anterior(const joint_pen_iter& jpi) const;
    const log_double&  posterior(const joint_pen_iter& jpi) const;
    const log_double&  posterior_with_mate(const member_type& mate_index, const joint_pen_iter& jpi) const;
    const log_double&  posterior_except_mate(const member_type& mate_index, const joint_pen_iter& jpi) const;

  private:
    size_t  phenoset_index(const joint_pen_iter& jpi) const;
    vector<log_double>&  posteriors_with_mate(size_t mate) const;
    vector<log_double>&  posteriors_except_mate(size_t mate) const;
    
    // Data members.
    std::map<size_t, size_t>  t_indices;         // <trait genotype id, matrix row x col count>
    std::map<size_t, size_t>  m_indices;         // <marker genotype id, matrix col>
    
    mutable vector<log_double>  my_anteriors;
    mutable vector<log_double>  my_posteriors;
    mutable vector<posterior_vector>  my_posteriors_with_mate;
    mutable vector<posterior_vector>  my_posteriors_except_mate;
};

#include "lodlink/cache.ipp"
}
}

#endif

