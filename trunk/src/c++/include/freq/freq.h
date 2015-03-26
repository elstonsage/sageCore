#ifndef FREQ_H
#define FREQ_H

#include <cassert>   
#include <functional>
#include <fstream>   
#include <iomanip>   
#include <iostream>
#include <limits>
#include <string>

#include "LSF/LSF.h"
#include "LSF/LSFfile.h"
#include "LSF/LSFsymbol.h"
#include "LSF/parse_ops.h"
#include "app/SAGEapp.h"  
#include "rped/rpfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h" 
#include "data/SAGEdata.h"
#include "freq/Output.h"
#include "freq/Parser.h"

namespace SAGE {
namespace FREQ {

/// \brief Defines freq application derived from SAGEapp.
///  
class freq : public APP::SAGEapp
{
public:

  /// @name Constructors
  //@{
  
    ///
    /// Constructor
    freq(int argc = 0, char **argv = NULL);
    
  //@}

  // Run the application
  virtual int main();
};

/// \brief Freq-specialized version of APP::SAGE_Data
///
class freq_data : public APP::SAGE_Data
{
public:

  freq_data(const std::string & program_name, bool debug);

  void  input(int argv, char** argv);
  virtual bool read_analysis() { return true; }
};

}// End namespace FREQ
} // End namespace SAGE
#endif
