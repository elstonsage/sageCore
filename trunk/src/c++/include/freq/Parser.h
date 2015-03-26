#ifndef FREQ_PARSER_H
#define FREQ_PARSER_H

#include <vector>
#include <string>
#include <set>
#include <sstream>
#include "LSF/LSF.h"
#include "app/ParsingFunctions.h"
#include "error/errorstream.h"
#include "error/errormanip.h" 
#include "LSF/LSFsymbol.h"
#include "freq/Estimator.h"
#include "freq/Configuration.h"

namespace SAGE {
namespace FREQ {

/// \brief Extends the APP::BasicParser for parsing input files
///
class Parser
{
public:

  ///
  /// Parses the options in the given LSFBase* and sets them accordingly in the given Estimator.
  ///
  /// NOTE: Assumes that the given LSFBase* is non-null and has a non-null List().
  static Configuration parseFreqBlock(size_t analysis_num, const LSFBase * param, const FPED::FilteredMultipedigree & mp, cerrorstream & errors = SAGE::sage_cerr);

  ///
  /// Creates a default configuration object.
  static Configuration createDefaultConfig(const FPED::FilteredMultipedigree & mp);
};

}// End namespace FREQ
} // End namespace SAGE

#endif
