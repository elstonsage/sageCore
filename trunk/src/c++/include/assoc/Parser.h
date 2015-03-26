#ifndef ASSOC_PARSER_H
#define ASSOC_PARSER_H
//=============================================================================
// File:      Parser.h
//                                                                          
// Author:    Stephen Gross
//
// History    11/14/7  Added code to parse summary_display subblock.    -djb
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//=============================================================================


#include <vector>
#include <string>
#include <limits>
#include <algorithm>
#include <iostream>
#include <sstream>
#include "globals/SAGEConstants.h"
#include "rped/rped.h"
#include "app/aparser.h"
#include "LSF/LSF.h"
#include "LSF/Attr.h"
#include "LSF/LSFfile.h"
#include "maxfunapi/maxfunapi.h"
#include "assoc/Configuration.h"
#include "app/ParsingFunctions.h"
#include "fped/fped.h"

namespace SAGE  {
namespace ASSOC {

class Parser
{
  public:

    // NOTE: Assumes that the given LSFBase* is non-null and has a non-null List().
    //
    static Configuration parseAssocBlock(size_t analysis_num, const LSFBase* param, 
                                         const RPED::MultiPedigree& mp, 
                                         std::ostream& info, cerrorstream& errors);
            
  private:
  
    static void parseTitle(Configuration& config, const LSFBase* param, cerrorstream& errors);
    static void parsePrimaryTrait(Configuration& config, const RPED::MultiPedigree& mp, const LSFBase* param, cerrorstream& errors);
    static void parseBatch(Configuration& config, const RPED::MultiPedigree& mp, cerrorstream& errors);
      static bool modelNameInUse(const string& name, Configuration& config, cerrorstream& errors);
    static void parseCovariate(Configuration& config, const RPED::MultiPedigree& mp, 
                               const LSFBase* param, std::ostream& info, cerrorstream& errors);
    static void parseIntercept(Configuration& config, const LSFBase* param, std::ostream& info, cerrorstream& errors);
    static void parseFixedEffect(const LSFBase* lsf_param, Parameter& param,
                                 cerrorstream& errors, bool& use_effect);
    static void parseUserEffect(const LSFBase* lsf_param, Configuration& config, const RPED::RefMPedInfo& rmp, cerrorstream& errors); 
      static bool  isCategoricalTrait(const string& trait_name, const RPED::RefMPedInfo& rmp);    
    static void parseTransSubModel(Configuration& config, 
                                    const LSFBase* param, ostream& info, cerrorstream& errors);
      static void parseLambda1(MFSUBMODELS::Transformation::Configuration& trans_config, 
                               const LSFBase* param, ostream& info, cerrorstream& errors);
      static void parseLambda2(Configuration& config, MFSUBMODELS::Transformation::Configuration& trans_config, 
                               const LSFBase* param, ostream& info, cerrorstream& errors);
                           
    static void parseResiduals(Configuration& config, const LSFBase* param, cerrorstream& errors);
      static void parseModel(Configuration& config, const LSFBase* param, cerrorstream& errors);
    static void parseAllowAveraging(Configuration& config, const LSFBase* param, ostream& info, cerrorstream& errors);
    static void  parseSummaryDisplay(Configuration& config, const LSFBase* param, cerrorstream& errors);
      static void  parseOrder(Configuration& config, const LSFBase* param, cerrorstream& errors);
      static void  parseFilters(Configuration& config, const LSFBase* param, cerrorstream& errors);
        static void  parseAll(Configuration& config, const LSFBase* param, cerrorstream& errors);
        static void  parseWald(Configuration& config, const LSFBase* param, cerrorstream& errors);
        static void  parseLRT(Configuration& config, const LSFBase* param, cerrorstream& errors);
        static void  parseLimitNumber(Configuration& config, const LSFBase* param, cerrorstream& errors);
    
    // - Added 6-6-7. djb
    //
    static void  prohibitBinaryTransformation(Configuration& config, cerrorstream& errors);
};

} 
} 

#endif
