#ifndef EVALFUNC_H
#define EVALFUNC_H

#include "util/RegexUtils.h"

#include "func/Function.h"
#include "func/FunctionParser.h"
#include "func/AdjustedTraitCreator.h"
#include "func/tai.h"
#include "func/twp.h"

namespace SAGE {
namespace FUNC {
 
class FunctionEvaluator
{
public:

  ///
  /// Evalutes a function block and populates a multipedigree with the calculated data.
  ///
  /// Function blocks generally take the form:
  /// \verbatim
  /// function
  /// {
  ///  trait = x, expression = "..."
  /// }
  /// \endverbatim
  /// 
  /// (See SAGE user manual for complete information on function block syntax).
  ///
  /// \param p The multipedigree
  /// \param errors The errorstream for messages
  /// \param params The pointer to the function block
  static void evaluateFunction(RPED::RefMultiPedigree& p, cerrorstream errors, const LSFBase* params);

private:

  ///
  /// Evaluates the specific function. This function is called by evaluateFunction because evaluateFunction
  /// actually (potentially) batches together a bunch of function blocks.
  static void evaluateSingleFunction(
    RPED::MultiPedigree & mp,
    cerrorstream errors,
    const FunctionParser::TraitData & trait_data);
};

} // End namespace FUNC
} // End namespace SAGE

#endif

