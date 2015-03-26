#ifndef SEGREG_DATATYPES_H
#define SEGREG_DATATYPES_H
//===================================================================
//
//  File:	segreg_datatypes.h
//
//  Author:	Stephen Gross
//
//  History:    sag Initial implementation		Jul 26 01
//
//  Copyright (c) 2001, R. C. Elston
//  All rights reserved
//===================================================================

#include "rped/rped.h"
#include "segreg/sub_model_base.h"
#include "segreg/ascertainment_sub_model.h"

namespace SAGE {
namespace SEGREG {

enum relationship_type  { undef = -1, no_rel = 0, sib = 1, parent_off = 2, mate = 3 };

class segreg_errors
{
  public:

    enum error_code { EVAL_OK = 0,
                      MCC_FAILED_TRANSFORM,
                      MCC_VARIANCE_INVALID,
                      BAD_INIT_PARAM,
                      ZERO_LIKELIHOOD,
                      BAD_THRESHOLDS,
                      BAD_LIKELIHOOD,
                      BAD_LEX_VALUE,
                      MLM_CORR_BOUNDARY
                    };
};


#define MAX_POLYGENOTYPE 11
#define NO_POLYGENIC_DATA 999

struct genetic_info
{
  genetic_info(genotype_index genotype_param     = index_AA,
               size_t         polygenotype_param = NO_POLYGENIC_DATA)
  { genotype = genotype_param; polygenotype = polygenotype_param; }
  genotype_index genotype;
  size_t         polygenotype;
};

// End namespace
}
}

#endif
