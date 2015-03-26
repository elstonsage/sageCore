#ifndef AO_DATATYPES_H
#define AO_DATATYPES_H
//=====================================================
//
//  File:	Datatypes.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//=====================================================


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "numerics/log_double.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "rped/rped.h"

namespace SAGE {
namespace AO   {

//=====================================================
//  Enums:
//=====================================================

enum trait_type    { type_suscept                   = 0, 
                     type_mean                      = 1, 
                     type_var                       = 2 };

//=====================================================
// Analysis type definitions and functions:
//=====================================================

const int SUSCEPTIBILITIES_MASK   = 1; // Binary: First bit
const int SUSCEPTIBILITIES_EQUAL  = 0; // Binary: x0
const int SUSCEPTIBILITIES_FREE   = 1; // Binary: x1
const int TRUNCATION_MASK         = 2; // Binary: Second bit
const int NO_TRUNCATION           = 0; // Binary: 0x
const int USE_TRUNCATION          = 2; // Binary: 1x
const int num_of_analysis_types   = 4;

inline bool SusceptibilitiesEqual  (int t) { return (t & SUSCEPTIBILITIES_MASK) == SUSCEPTIBILITIES_EQUAL; }
inline bool SusceptibilitiesFree   (int t) { return (t & SUSCEPTIBILITIES_MASK) == SUSCEPTIBILITIES_FREE;  }
inline bool UseTruncation          (int t) { return (t & TRUNCATION_MASK)       == USE_TRUNCATION;         }
inline bool NoTruncation           (int t) { return (t & TRUNCATION_MASK)       == NO_TRUNCATION;          }

string ConvertAnalysisType(int analysis_type, bool pooled);

//=====================================================
// Constants definitions:
//=====================================================

const bool INCLUDE_DETAILED = true;
const bool EXCLUDE_DETAILED = false;

#define AO_CLASS_TYPE                   "_CLASS_"
#define AO_DEFAULT_ANALYSIS_TITLE       "AGEON Analysis "
#define AO_DEFAULT_PEDIGREE_FILE        "_random_AO_pedigree.ped"
#define AO_DEFAULT_USE_TRUNCATION       false
#define AO_DEFAULT_USE_ADJUSTMENT       false

#define AO_DEFAULT_EPSILON              0.00001

#define AO_DEFAULT_USE_SUSCEPT_BASE     true
#define AO_DEFAULT_USE_MEAN_BASE        true
#define AO_DEFAULT_USE_VAR_BASE         true

#define AO_DEFAULT_LAMBDA1_TYPE         MAXFUN::Parameter::INDEPENDENT
#define AO_DEFAULT_LAMBDA2_TYPE         MAXFUN::Parameter::FIXED

#define AO_DEFAULT_NUM_OF_CLASSES       6

#define AO_LOWER_BOUND_GEN_SUSCEPT     -std::numeric_limits<double>::infinity()
#define AO_DEFAULT_GEN_SUSCEPT          0.0000
#define AO_UPPER_BOUND_GEN_SUSCEPT      std::numeric_limits<double>::infinity()

#define AO_LOWER_BOUND_MEAN            -std::numeric_limits<double>::infinity()
#define AO_DEFAULT_MEAN	                0.0
#define AO_UPPER_BOUND_MEAN             std::numeric_limits<double>::infinity()

#define AO_LOWER_BOUND_VARIANCE         0.0000
#define AO_DEFAULT_VARIANCE             0.0005
#define AO_UPPER_BOUND_VARIANCE         std::numeric_limits<double>::infinity()

#define AO_LOWER_BOUND_LAMBDA1         -std::numeric_limits<double>::infinity()
#define AO_DEFAULT_LAMBDA1              1.00000
#define AO_UPPER_BOUND_LAMBDA1          std::numeric_limits<double>::infinity()

#define AO_LOWER_BOUND_LAMBDA2         -std::numeric_limits<double>::infinity()
#define AO_DEFAULT_LAMBDA2              0.05000
#define AO_UPPER_BOUND_LAMBDA2          std::numeric_limits<double>::infinity()

#define AO_LOWER_BOUND_COVARIATE       -std::numeric_limits<double>::infinity()
#define AO_DEFAULT_COVARIATE            0.00000
#define AO_UPPER_BOUND_COVARIATE        std::numeric_limits<double>::infinity()

const double ONE_OVER_TWO_PI = 1.0 / sqrt(2.0 * M_PI);

//==================================================================
// Typedefs:
//==================================================================

typedef RPED::RefMultiPedigree::member_const_pointer       member_const_pointer;
typedef RPED::RefMultiPedigree::member_const_iterator      member_const_iterator;
typedef RPED::RefMultiPedigree::pedigree_const_iterator    pedigree_const_iterator;
typedef RPED::RefMultiPedigree::subpedigree_const_pointer  subpedigree_const_pointer;
typedef RPED::RefMultiPedigree::subpedigree_const_iterator subpedigree_const_iterator;

//=====================================================
//  forward class declarations:
//=====================================================

class Model;

//=====================================================
//  class ao_errors
//=====================================================
class ao_errors
{
  public:

    enum error_code 
    {
      BOO = 0
    };
};

}} // End namespace

#endif
