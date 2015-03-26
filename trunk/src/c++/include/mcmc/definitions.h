#ifndef MCMC_DEFS_H
#define MCMC_DEFS_H

//==========================================================================
//  File:      definitions.h
//
//  Author:    Yeunjoo Song
//
//  History:   Initial implementation.                              May. 04
//
//  Notes:     This header file contains various type definition statements
//             used through out the mcmc files.
//
//  Copyright (c) 2004 R.C. Elston
//    All Rights Reserved
//==========================================================================

#include "gelim/pedigree_region.h"
#include "fped/fped.h"
#include "pairs/relpair.h"
#include "rped/genome_description.h"
#include "containers/muxmap.h"
#include "containers/bitfield.h"
#include "numerics/print_util.h"
#include "util/dots.h"
#include "globals/SAGEConstants.h"

namespace SAGE
{

namespace MCMC
{

// Default mcmc_parameter options.
//
#define DEFAULT_TT_MAX_TRANS      40
#define DEFAULT_T0_WEIGHT         0.55
#define DEFAULT_T1_WEIGHT         0.4
#define DEFAULT_T2_WEIGHT         0.05
#define DEFAULT_LOCAL_WEIGHT      0.75 //weight to choose local marker/person

#define DEFAULT_MCMC_STEPS        200000
#define DEFAULT_DEMEMO_STEPS      50000
#define DEFAULT_BATCH_COUNT       100

enum mode_type { SINGLEPOINT, MULTIPOINT };

inline void print_bit_field(ostream& o, const bit_field& b)
{
  for( size_t i = 0; i != b.size(); ++i )
    if( b[i] ) o << 'X';
    else       o << '.';

  return;
}

} // end of namespace MCMC

} // end of namespace SAGE

#endif
