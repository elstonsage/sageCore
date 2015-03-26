#ifndef ASSOC_DATATYPES_H
#define ASSOC_DATATYPES_H
//======================================================
//
//  File:	Datatypes.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================


#include <stdexcept>
#include "mped/mp.h"
#include "mped/sp.h"
#include "rped/rped.h"
#include "sampling/sampling.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846

#endif

namespace SAGE  {
namespace ASSOC {

const size_t  effect_count = 5;
const string  effect_names[] = { "Random", "Polygenic", "Family", "Sibling", "Marital" };

inline bool
fixed_effect(const string& effect_name)
{
  for(size_t i = 0; i < effect_count; ++i)
  {
    if(effect_names[i] == effect_name)
    {
      return  true;
    }
  }
  
  return  false;
}

typedef SAMPLING::PartitionedMemberDataSample  Sampledata;

// The pre-compiler define "ASSOC_FORMULAS_DEBUG", when defined, enables
// structural/algorithmic reporting:
//
// #define ASSOC_FORMULAS_DEBUG  

// ASSOC_FORMULAS_DEBUG_l2 enables reporting of calculation portion of equations:
//
// #define ASSOC_FORMULAS_DEBUG_L2  

class BadLikelihood : public runtime_error
{ 
  public:
    BadLikelihood(const string& model_name)
          : runtime_error(model_name)
    {}
};

class InsufficientData : public runtime_error
{ 
  public:
    InsufficientData(const string& model_name)
          : runtime_error(model_name)
    {}
};

class NonConvergenceWithTransformation : public runtime_error
{ 
  public:
    NonConvergenceWithTransformation(const string& model_name)
          : runtime_error(model_name)
    {}
};

enum DependentTraitTypeEnum
{
  QUANTITATIVE = 0,
  BINARY       = 1
};

const double sqrt_2   = sqrt(2.0);
const double sqrt_2PI = sqrt(2.0 * M_PI);
const double PI_2 = 6.28318530717958647692;

}
} 

#endif
