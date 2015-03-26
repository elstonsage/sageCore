#ifndef PEELING_CACHE2_H
#define PEELING_CACHE2_H
//
//  Caching intermediate peeling results
//
//  Author: Geoff Wedig (wedig@darwin.cwru.edu)
//
//  History:   (APPelements)
//             0.1 gcw Initial Implementation                Aug 28 96
//             1.0 gcw Updated and debugged for beta release Dec 30 96
//             (caching)
//             0.1 gcw Converted to a general framework      Apr 18 01
//             gcw/djb created cache2.h for subpedigrees     Oct 24 02
//
//  Copyright (c) 2001  R.C. Elston
//


#include "mped/mp.h"

namespace SAGE
{

// The caching object store intermediate results of anterior and posterior
// values based upon the Fernado, Stricker, and Elston (1993) paper.  The
// caches are heavily templatized, and require specialization.  Like the
// peeling classes themselves, they provide infrastructure rather than a
// complete implementation.  Subclassing of these classes, filling in the
// functions required, will be necessary to implement a specific
// computation.

namespace peeling
{

/// stores the results of the anterior and posterior values of an individual

/** The individual cache is primarily a base class.  We may create our own variants
 *  by either specialization or independent instance since it is a template
 *  parameter of the peeling_cache and the peeler.  Any class that supports
 *  the twelve public functions can be used as a cache.
 *
 *  Because it is a base class, and because the data type is undefined,
 *  storage of the results is up to the class.  Results are passed as const
 *  references to allow for large result data structures in this base class.
 * 
 *  The individual cache is created and stored by the peeling cache.  It is
 *  mainly used to reference any individual specific data which may be
 *  required. Note that it will generally not be necessary for algorithms to
 *  call the functions directly.  That's handled by the peeling object
 *  itself.  However individual algorithms must implement the twelve access
 *  functions. Primarily, these just store the intermediate results, and
 *  should be relatively easy to implement.
 * 
 *  The access functions can be grouped into three categories:
 * 
 *  -#  Test for existance - boolean test to see if there is a cached value
 *      for the data or not.
 *  -#  Cached result lookup - returns the cached result.
 *  -#  Cached result setting - Set the cached result.
 * 
 *  There are 4 functions in each category:
 * 
 *  -# Anterior           - The cached result for the part of the pedigree
 *                          anterior to the individual
 *  -# Posterior          - The cached result for the part of the pedigree
 *                          posterior to the individual, including all mates.
 *  -# Posterior w. mate  - The cached result for the part of the pedigree
 *                          posterior to the individual, but only with a
 *                          particular mate.
 *  -# Posterior w/o mate - The cached result for the part of the pedigree
 *                          posterior to the individual, but with all mates except
 *                          a particular mate.
 */
template <class Data, class Result>
class individual_cache
{
public:

  typedef Data   data_type;
  typedef Result result_type;

  typedef RPED::RefMultiPedigree::member_type member_type;

  individual_cache() { }
  
  ~individual_cache() { }

  // Functions to be implemented by Algorithm:

    // Information retrieval functions

  bool anterior_cached              (const data_type&) const;
  bool posterior_cached             (const data_type&) const;
  bool posterior_with_mate_cached   (const member_type& mate, const data_type&) const;
  bool posterior_except_mate_cached (const member_type& mate, const data_type&) const;

  result_type& anterior              (const data_type&);
  result_type& posterior             (const data_type&);
  result_type& posterior_with_mate   (const member_type& mate, const data_type&);
  result_type& posterior_except_mate (const member_type& mate, const data_type&);

  const result_type& anterior              (const data_type&)   const;
  const result_type& posterior             (const data_type&)   const;
  const result_type& posterior_with_mate   (const member_type& mate, const data_type&) const;
  const result_type& posterior_except_mate (const member_type& mate, const data_type&) const;

protected:

  // Actual data to store results are filled in by individual specializations

};

/// stores the cached individual values.
/** As such, it can work with any individual cache without specific
 *  specialization.  It stored a vector equal to the number of individuals
 *  in the pedigree section which it is given.  Individual caches are
 *  accessed by the individual index in the RPED::pedigree_section.
 */
template <class IndCache>
class peeling_cache
{
public:

  typedef IndCache                       individual_cache_type;
  typedef typename IndCache::data_type   data_type;
  typedef typename IndCache::result_type result_type;

  typedef RPED::RefMultiPedigree::subpedigree_type subped_type;
  typedef RPED::RefMultiPedigree::member_type   member_type;

  peeling_cache(const subped_type&);
  
  ~peeling_cache();

  individual_cache_type& get_individual_cache(const member_type& index);

protected:

  typedef vector<individual_cache_type> ind_data;

  const subped_type& my_subpedigree;

  ind_data my_data;
};

template <class IndCache>
inline
peeling_cache<IndCache>::peeling_cache(const subped_type& p)
  : my_subpedigree(p)
{
  my_data.resize(my_subpedigree.member_count());
}

template <class IndCache>
inline
peeling_cache<IndCache>::~peeling_cache()
{ }

template <class IndCache>
inline typename peeling_cache<IndCache>::individual_cache_type&
peeling_cache<IndCache>::get_individual_cache(const member_type& index)
{
  return my_data[index.subindex()];
}

}

}

#endif
