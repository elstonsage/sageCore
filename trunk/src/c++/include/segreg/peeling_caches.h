#ifndef PEELING_CACHES_H
#define PEELING_CACHES_H
//===================================================================
//
//  File:	peeling_caches.h
//
//  Author:	Stephen Gross
//
//  History:	sag Initial implementation		Jul 30 2001
//
//  Copyright (c) 2001 R. C. Elston
//  All rights reserved
//
//===================================================================


#include "error/internal_error.h"
#include "peeling/cache3.h"
#include "numerics/log_double.h"
#include "segreg/segreg_datatypes.h"
#include "segreg/sub_model_base.h"
#include "segreg/types/TypeDescription.h"

namespace SAGE
{
namespace peeling
{

template<>
class individual_cache<SEGREG::TypeDescription::State,log_double>
{
  public:
    typedef SEGREG::TypeDescription::State data_type;
    typedef log_double                     result_type;

    typedef FPED::Member member_type;

    individual_cache();
    
    ~individual_cache()
    { }

    void set_member(const member_type& m);

    bool anterior_cached(const data_type&) const;
    bool anterior_with_mate_cached(const member_type& mate, const data_type&, const data_type&) const;
    bool posterior_cached(const data_type&) const;
    bool posterior_with_mate_cached(const member_type& mate, const data_type&) const;
    bool posterior_except_mate_cached(const member_type& mate, const data_type&) const;
    bool ppl_cached(const member_type& mate, const data_type&) const;
    bool ppl_with_mate_cached(const member_type& mate, const data_type&, const data_type&) const;

    result_type & anterior(const data_type &);
    result_type & anterior_with_mate( const member_type& mate, const data_type &, const data_type&);
    result_type & posterior(const data_type &);
    result_type & posterior_with_mate (const member_type& mate, const data_type&);
    result_type & posterior_except_mate (const member_type& mate, const data_type&);
    result_type & ppl(const member_type& mate, const data_type&);
    result_type & ppl_with_mate(const member_type& mate, const data_type&, const data_type&);

    const result_type & anterior( const data_type &) const;
    const result_type & anterior_with_mate(const member_type& mate, const data_type &, const data_type&) const;
    const result_type & posterior(const data_type &) const;
    const result_type & posterior_with_mate(const member_type& mate, const data_type&) const;
    const result_type & posterior_except_mate(const member_type& mate, const data_type&) const;
    const result_type & ppl(const member_type& mate, const data_type&) const;
    const result_type & ppl_with_mate(const member_type& mate, const data_type&, const data_type&) const;

  protected:

    struct mate_data
    {
      mate_data(const member_type* mate_param = NULL)
      {
        mate = mate_param;
        
        double qNaN = numeric_limits<double>::quiet_NaN();

        for(int i=0; i<3; ++i)
        {
          partial_parental_likelihood[i] = qNaN;
          posterior_with_mate        [i] = qNaN;
          posterior_except_mate      [i] = qNaN;

          for(int j = 0; j < 3; ++j)
          {
            anterior_with_mate[i][j] = qNaN;
            ppl_with_mate     [i][j] = qNaN;
          }
        }
      }

      mate_data(const mate_data& md)
      {
        mate = md.mate;
        
        for(int i=0; i<3; ++i)
        {
          partial_parental_likelihood[i] = md.partial_parental_likelihood[i] ;
          posterior_with_mate        [i] = md.posterior_with_mate        [i] ;
          posterior_except_mate      [i] = md.posterior_except_mate      [i] ;

          for(int j = 0; j < 3; ++j)
          {
            anterior_with_mate[i][j] = md.anterior_with_mate[i][j];
            ppl_with_mate     [i][j] = md.ppl_with_mate     [i][j];
          }
        }
      }

      mate_data& operator=(const mate_data& md)
      {
        if(this == &md) return *this;

        mate = md.mate;
        
        for(int i=0; i<3; ++i)
        {
          partial_parental_likelihood[i] = md.partial_parental_likelihood[i] ;
          posterior_with_mate        [i] = md.posterior_with_mate        [i] ;
          posterior_except_mate      [i] = md.posterior_except_mate      [i] ;

          for(int j = 0; j < 3; ++j)
          {
            anterior_with_mate[i][j] = md.anterior_with_mate[i][j];
            ppl_with_mate     [i][j] = md.ppl_with_mate     [i][j];
          }
        }

        return *this;
      }

      mutable const member_type* mate;
    
      /// partial parental likelihood(g) = anterior(g) * posterior_except_mate(g)
      mutable log_double partial_parental_likelihood[3];

      /// partial parental likelihood(g,h) = anterior_with_mate(g,h) * posterior_except_mate(g)
      mutable log_double ppl_with_mate[3][3];

      mutable log_double anterior_with_mate[3][3];

      mutable log_double posterior_with_mate[3];
      mutable log_double posterior_except_mate[3];
    };

    typedef vector<mate_data>::const_iterator   mate_data_const_iterator;
    typedef vector<mate_data>::iterator         mate_data_iterator;

    mate_data_const_iterator find(const member_type& mate) const;
    mate_data_iterator       find(const member_type& mate);

    log_double                      my_anterior[3];
    log_double                      my_posterior[3];
 
    vector<mate_data> my_mate_data;

};

typedef SEGREG::genotype_index genotype_index;
typedef SEGREG::genetic_info   genetic_info;

template<>
class individual_cache<genetic_info,log_double>
{
  public:
    typedef genetic_info data_type;
    typedef log_double   result_type;

    typedef FPED::Member member_type;

    individual_cache();
    individual_cache(const individual_cache&);

    bool anterior_cached              (                   const data_type&) const;
    bool posterior_cached             (                   const data_type&) const;
    bool posterior_with_mate_cached   (const member_type& mate, const data_type&) const;
    bool posterior_except_mate_cached (const member_type& mate, const data_type&) const;

    result_type & anterior              (                   const data_type&);
    result_type & posterior             (                   const data_type&);
    result_type & posterior_with_mate   (const member_type& mate, const data_type&);
    result_type & posterior_except_mate (const member_type& mate, const data_type&);

    const result_type & anterior              (                   const data_type&) const;
    const result_type & posterior             (                   const data_type&) const;
    const result_type & posterior_with_mate   (const member_type& mate, const data_type&) const;
    const result_type & posterior_except_mate (const member_type& mate, const data_type&) const;

  protected:
    log_double my_anterior [3][MAX_POLYGENOTYPE];
    log_double my_posterior[3][MAX_POLYGENOTYPE];

    struct polygenotype_mate_info
    {
      polygenotype_mate_info(const member_type* mate_param = NULL);
      polygenotype_mate_info(const polygenotype_mate_info&);
      polygenotype_mate_info& operator=(const polygenotype_mate_info&);
      const member_type* mate;
      log_double with_mate  [3][MAX_POLYGENOTYPE];
      log_double except_mate[3][MAX_POLYGENOTYPE];
    };

    vector<polygenotype_mate_info> posterior_on_mate;

    typedef vector<polygenotype_mate_info>::const_iterator gen_const_iterator;
    typedef vector<polygenotype_mate_info>::iterator       gen_iterator;
};

} // End namespace
}

#include "segreg/peeling_caches.ipp"

#endif

