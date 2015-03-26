#ifndef SEGREG_SUB_MODEL_BASE_H
#define SEGREG_SUB_MODEL_BASE_H
//============================================================================
// File:      sub_model_base.h
//                                                                          
// Author:    Geoff Wedig (wedig@darwin.cwru.edu)
//                                                                          
// History:   0.1 gcw Initial Implementation                    Apr 2001   
//                djb reformatted, added  misc. global 
//                    declarations/definitions                  Jun 2001
//                                                                          
// Notes:     The purpose of this file is to create an easy method of inclusion of
//            different sub-models into maxfun.  While it is currently in the SEGREG
//            hierarchy, it is actually quite general, and after some more development,
//            we may wish to include it in the maxfun library as a useful model generator.
//
//            In essence, a 'sub-model' is a class that contains a set of parameters for
//            maxfun to use.  There is only a minimal public interface to the sub-model; most of 
//            the interface is supplied in the derived class.  The sub-model provides the initial 
//            estimate, bounds, and status type for each parameter to a sub-model sequencer.  The
//            sub-model sequencer then inserts this into the list of parameters on its maxfun
//            object.
//
//            In addition to the construction phase, the sub-model sequencer wrappers the
//            function calls maxfun makes to MaxFunction.  This wrapper copies out the
//            data from the maxfun parameter list back into the sub-models, which should do
//            any dependency checking.  The evaluate() function can then use the sub-model
//            instead of the raw vector supplied to it from maxfun.  This allows a
//            nearly transparent use of the sub-models, without the need for manual
//            sequencing, which is prone to errors.
//
//            In addition, the modular framework of the sub-model makes it quite easy to
//            maintain. With proper use (each sub-model dealing with only a few
//            parameters), it becomes nearly trivial to take out one sub-model and replace
//            it with another, without disturbing code that does not rely on that
//            sub-model. With the raw vector, this is not the case.  Knowing which parameter is
//            which index is of prime importance.  The sub-model interface takes care of that
//            detail so that the programmer doesn't have to.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



#include <math.h>
#include <assert.h>
#include <cmath>
#include <list>
#include <vector>
#include <functional>
#include <algorithm>
#include <limits>
#include <ostream>
#include <sstream>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "LSF/parse_ops.h"
#include "app/SAGEapp.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/maxfunapi.h"
#include "rped/rped.h"
#include "segreg/segreg_options.h"
#include "util/StringUtils.h"

using std::string;
using std::ostream;

namespace SAGE
{
namespace SEGREG
{

/// \brief SEGREG Specific Maxfun Submodel
///
/// This class extends the MAXFUN::Submodel to include things needed for SEGREG.
/// At present (2004-09-21), this means allowing ostreams access to certain
/// internals for debugging and a few functions from the original, maxfun
/// sub_model which were not included in the MAXFUNAPI::Submodel
class SegregSubmodel : public MAXFUN::Submodel
{
    friend ostream& operator<< (ostream& out, const SegregSubmodel& sm);

  public:
  
    /// Default Constructor
    ///
    SegregSubmodel (cerrorstream& errors=sage_cerr);

    /// Copy Constructor
    ///
    SegregSubmodel (const SegregSubmodel &other);
    
    /// Copy Operator
    ///
    SegregSubmodel& operator= (const SegregSubmodel &other);    

    /// Destructor.  Virtual due to virtual interface of MAXFUN::Submodel
    ///
    virtual ~SegregSubmodel();

    virtual string name()               const = 0;
    virtual string option_description() const = 0;
};

ostream& operator<< (ostream& out, const SegregSubmodel& sm);
  
const int  DUMP_PRECISION = 12;

const int NUM_OF_TYPES = 3;

// - misc. enumerations.
//
enum gi_type { index_AA = 0, index_AB, index_BB, index_INVALID };

class genotype_index
{
  public:

    inline genotype_index(unsigned int);
    inline genotype_index(gi_type = index_AA);
    inline genotype_index(const genotype_index&);

    inline genotype_index& operator=(unsigned int);
    inline genotype_index& operator=(gi_type);
    inline genotype_index& operator=(const genotype_index&);
  
    inline operator unsigned int () const;

    inline genotype_index& operator ++();
    inline genotype_index  operator ++(int);

    inline bool operator==(const genotype_index&) const;
    inline bool operator!=(const genotype_index&) const;

    inline bool operator==(gi_type) const;
    inline bool operator!=(gi_type) const;

  protected:

    gi_type my_type;

  //lint --e{1739} <- We don't want == and != operators as non-members
};

enum genotype_info { no_geno = 0, 
                     AA = 1, AB = 2, BB = 4,
                     AA_AB = 3, AA_BB = 5, AB_BB = 6,
                     all = 7 };

// ??? What's this for? GCW 2003-05-05
enum freq_index { index_freq_A = 3, index_corr };                 

genotype_info  index_2_info(genotype_index index);
genotype_info  total_info(const model_input& mi_AA, const model_input& mi_AB, 
                          const model_input& mi_BB);
genotype_index  third_index(genotype_index index_one, genotype_index index_two);

//lint -e{1761}
ostream&  operator<<(ostream& out, const model_input& mi);
string  value_phrase(const maxfun_parameter& mi);

}
}

#include "segreg/sub_model_base.ipp"

#endif


