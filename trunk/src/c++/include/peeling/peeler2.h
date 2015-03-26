#ifndef PEELER2_H
#define PEELER2_H
//
//  Peeling infrastructure for pedigree algorithms
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   (APP)
//             0.1 gcw Initial Implementation                Aug 28 96
//             1.0 gcw Updated and debugged for beta release Dec 30 96
//             2.0 gcw Redesign to use the new Genotype      Sep 18 97
//                     Elimination Routines
//             (peeler)
//             0.1 gcw Generalized to handle a variety of    Apr    01
//                     algorithms
//             gcw/djb created peeler2.h for subpedigrees    Oct 24 02
//
//  Copyright (c) 2001  R.C. Elston
//


#include "rped/rped.h"
#include "peeling/cache2.h"

namespace SAGE
{

namespace peeling
{

/// the main object in the pedigree peeling algorithm.
/** It is an infrastructure for creating algorithms based upon the Fernado,
 *  Stricker, and Elston paper on efficient calculations of Elston
 *  Stewart-based algorithms using anteriors and posteriors.
 * 
 *  As it is infrastructure, it is not, in itself, a complete algorithm, but
 *  includes a variety of functions which must be filled in by the algorithm
 *  developer.  It is also templatized to allow multiple different algorithms
 *  to be computed.  Different algorithm require different data to be passed
 *  up and down anterior/posterior chains.  The peeler, and its helper
 *  classes are templatized upon a Data class, the data which is needed to
 *  compute or reference a particular set of anterior/posterior values, and a
 *  Result class, the data which is returned from an anterior/posterior
 *  calculation.
 * 
 *  The public interface of the peeler consists of inline functions for the
 *  caching (using the peeling_cache from cache.h) and lookup of previously
 *  computed results.  All actual computation is done using the private functions
 *  which must be filled in by the user.  This is accomplished by specializing
 *  the class based upon the Data and Result type and adding those functions. 
 */
template <class Data, class Result, class IndCache = individual_cache<Data, Result> >
class peeler
{
public:

  typedef Data   data_type;
  typedef Result result_type;

  typedef RPED::RefMultiPedigree::subpedigree_type  subped_type;
  typedef RPED::RefMultiPedigree::member_type       member_type;

  peeler(const subped_type&);

  virtual ~peeler();

  /// @name External functions.
  /** These call the internal functions and cache results.  There should be
   *  no reason to change these when specializing the peeler for a particular
   *  algorithm. */
  //@{
  const result_type& anterior    (const member_type& ind, const data_type& g);
  const result_type& posterior   (const member_type& ind, const data_type& g);

  const result_type& posterior_with_mate   (const member_type& ind,
                                            const member_type& mate,
                                            const data_type& g);
  const result_type& posterior_except_mate (const member_type& ind, 
                                            const member_type& mate,
                                            const data_type& g);
  //@}

  const subped_type& get_subpedigree () { return my_subpedigree; }

protected:

  typedef IndCache                             individual_cache_type;
  typedef peeling_cache<individual_cache_type> peeling_cache_type;

  const subped_type&  my_subpedigree;
  
  peeling_cache_type  my_cache;

  /// Functions needing to be implemented by the algorithm developer
  //@{
  
  virtual const result_type& internal_anterior    (const member_type& ind, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior   (const member_type& ind, const data_type& g, result_type&) =0;

  virtual const result_type& internal_posterior_with_mate
                                 (const member_type& ind, const member_type& mate, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior_except_mate
                                 (const member_type& ind, const member_type& mate, const data_type& g, result_type&) =0;

  // These also need implemented.  They can safely be split out from the
  // algorithm development.  It can be assumed that all anteriors and
  // posteriors actually contain a valid subpedigree.

  virtual const result_type& internal_anterior_terminal(const member_type& ind, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior_terminal(const member_type& ind, const data_type& g, result_type&) = 0;
  
  //@}
};

/// makes it possible to create a peeler with no cache without overhead of tests
class no_cache {};

/// specialization of peeler for no_cache
template <class Data, class Result>
class peeler<Data, Result, no_cache>
{
public:

  typedef Data   data_type;
  typedef Result result_type;

  typedef RPED::RefMultiPedigree::subpedigree_type  subped_type;
  typedef RPED::RefMultiPedigree::member_type    member_type;

  peeler(const subped_type&);

  virtual ~peeler();

  /// @name External functions.
  /// These call the internal functions and cache results.  There should be
  /// no reason to change these when specializing the peeler for a particular
  /// algorithm.
  //@{
  result_type anterior    (const member_type& ind, const data_type& g);
  result_type posterior   (const member_type& ind, const data_type& g);

  result_type posterior_with_mate   (const member_type& ind, const member_type& mate, const data_type& g);
  result_type posterior_except_mate (const member_type& ind, const member_type& mate, const data_type& g);
  //@}

protected:

  const subped_type& my_subpedigree;
  
  /// Functions needing to be implemented by the algorithm developer
  //@{
  
  virtual const result_type& internal_anterior    (const member_type& ind, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior   (const member_type& ind, const data_type& g, result_type&) =0;

  virtual const result_type& internal_posterior_with_mate
                                 (const member_type& ind, const member_type& mate, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior_except_mate
                                 (const member_type& ind, const member_type& mate, const data_type& g, result_type&) =0;

  // These also need implemented.  They can safely be split out from the
  // algorithm development.  It can be assumed that all anteriors and
  // posteriors actually contain a valid subpedigree.

  virtual const result_type& internal_anterior_terminal(const member_type& ind, const data_type& g, result_type&) =0;
  virtual const result_type& internal_posterior_terminal(const member_type& ind, const data_type& g, result_type&) = 0;
  
  //@}
};

#include "peeling/peeler2.ipp"

} // end namespace peeling

} // end namespace SAGE

#endif

