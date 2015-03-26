#ifndef OUTPUTFORMATTER_H
#define OUTPUTFORMATTER_H

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

namespace SAGE   {
namespace MAXFUN {

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

  /// @name Primary output sections
  //@{

    static OUTPUT::Table convertEstimates(const Results& results, bool detailed = true, bool p_val = true);

    static OUTPUT::Table convertMatrix(const Results& results);

  //@}

private:

    static string getWasSkippedMessage();
};

class JointTest
{
  public:
  
    JointTest(const Results & results1, const Results & results2);
  
    JointTest(const JointTest & other);
    
    JointTest & operator= (const JointTest & other);
  
    OUTPUT::Table summarizeAsTable(string title = "") const;

    double getPValue() const { return my_p_value; } 

  private:

    std::string my_H0_name;
    std::string my_H1_name;

    bool my_H0_valid;
    bool my_H1_valid;
    
    double my_H0_val;
    double my_H1_val;

    size_t my_degrees_of_freedom;
    double my_comp_val;
    double my_p_value;

};


} // End namespace MAXFUN
} // End namespace SAGE

#endif
