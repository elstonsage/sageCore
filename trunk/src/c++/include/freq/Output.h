#ifndef FREQ_OUTPUT
#define FREQ_OUTPUT

#include <iostream>
#include "output/Output.h"
#include "freq/Results.h"

namespace SAGE {
namespace FREQ {

class Output
{
public:

  static void generateOutput(const Results & results, const FPED::FilteredMultipedigree & fp);

private:

  static void generateAnalysisOutput(const Results & results, const FPED::FilteredMultipedigree & fp, bool include_detailed);

  static void generateLocusDescription(const Results & results, const FPED::FilteredMultipedigree & fp);
};

} // End namespace FREQ
} // End namespace SAGE

#endif
