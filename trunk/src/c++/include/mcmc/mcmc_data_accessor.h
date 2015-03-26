#ifndef  MCMC_DATA_ACCESSOR_H
#define  MCMC_DATA_ACCESSOR_H

//==========================================================================
//  File:    mcmc_data_accessor.h
//
//  Author:  Geoff Wedig
//           Yeunjoo Song
//
//  History: Version 0.90
//           1.0 Changed the class name MCMCDataStrategy
//               to mcmc_data_accessor &
//               updated to new libraries                        yjs May. 04
//
//  Notes:   This header defines accesses to individual's raw data 
//           for MCMC IBD simulatior.
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "mcmc/mcmc_meiosis_map.h"

namespace SAGE
{

namespace MCMC
{

class mcmc_data_accessor
{
  public:

    typedef vector<bit_field> indicator_type;

    mcmc_data_accessor(const McmcMeiosisMap& mmap, size_t locus_count); 
    ~mcmc_data_accessor() { }

    //input, output

    bool       is_valid_locus(size_t loc)               const;
    bool       is_geno_miss  (size_t id, size_t loc)    const;

    int        mother_bit(size_t ind, size_t locus)     const;
    int        father_bit(size_t ind, size_t locus)     const;

    size_t     individual_count()                       const;
    size_t     locus_count()                            const;

    void       set_valid_locus(size_t locus, bool = true);

    indicator_type&         get_indicator();
    const indicator_type&   get_indicator() const;

    bit_field&              get_indicator(size_t m);
    const bit_field&        get_indicator(size_t m)     const;
    
    const McmcMeiosisMap& get_mcmc_meiosis_map()      const;

    void                    dump_accessor(ostream& o)   const;

  private:

    const McmcMeiosisMap& my_mcmc_meiosis_map;

    indicator_type          my_indicator;

    size_t                  my_locus_count;

    vector<bool>            my_valid_loci;
};

#include "mcmc/mcmc_data_accessor.ipp"

} // end of namespace MCMC

} // end of namespace SAGE

#endif
