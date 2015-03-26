#ifndef ASSOC_OUTPUTFORMATTER_H
#define ASSOC_OUTPUTFORMATTER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include "numerics/cephes.h"
#include "numerics/functions.h"
#include "LSF/parse_ops.h"
#include "output/Output.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/CovarianceMatrix.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/Results.h"
#include "maxfunapi/Parameter.h"

namespace SAGE  {
namespace ASSOC {

/** \internal
  * \brief Formats Maxfun output.
  *
  * Given a Results object (and in some cases, two Results objects), this class formats
  * the requested portion of the output, and returns said output as a string. 
  *
  */
class OutputFormatter
{
public:
    static OUTPUT::Table convertEstimates(const MAXFUN::Results& results, double residuals_scaling_factor, bool detailed = true);
    static OUTPUT::Table convertMatrix(const MAXFUN::Results& results, double residuals_scaling_factor);

private:

    static double  covariance_adjustment(const string& group_name1, const string& param_name1,
                                         const string& group_name2, const string& param_name2,
                                         double residuals_scaling_factor);
    static string getWasSkippedMessage();
};

} // End namespace ASSOC
} // End namespace SAGE

#endif
