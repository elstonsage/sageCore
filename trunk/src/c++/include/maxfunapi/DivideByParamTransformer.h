#ifndef DIVIDE_BY_PARAM_TRANSFORMER_H
#define DIVIDE_BY_PARAM_TRANSFORMER_H

#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

/** \class DivideByParamTransformer
 *
 * This transformer will automatically divide one parameter's current estimate by another parameter's current
 * estimate.
 */
class DivideByParamTransformer : public Transformer
{
public:
  ///
  /// Constructor.
  /// \param param_idx The index number of the parameter whose current estimate will
  /// serve as the denominator.
  explicit DivideByParamTransformer(int param_idx) { my_param_idx = param_idx; }

  DivideByParamTransformer(const DivideByParamTransformer & other) { my_param_idx = other.my_param_idx; }

  DivideByParamTransformer& operator=(const DivideByParamTransformer & other) { my_param_idx = other.my_param_idx; return *this;}

  virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const;
  virtual double transformToReportedScale     (const ParameterMgr * mgr, double val) const;

private:
  int my_param_idx;
};

} // End namespace MAXFUN
} // End namespace SAGE

#endif
