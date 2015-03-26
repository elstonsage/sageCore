#ifndef COVARIATETRANSFORMER_H
#define COVARIATETRANSFORMER_H

#include "maxfunapi/ParameterMgr.h"

namespace SAGE   {
namespace MAXFUN {

/** \class CovariateTransformer
 */
class CovariateTransformer : public Transformer
{
public:
  explicit CovariateTransformer(double val) { my_stdev = val; }
  CovariateTransformer(const CovariateTransformer & other) { my_stdev = other.my_stdev; }
  CovariateTransformer& operator=(const CovariateTransformer & other) { my_stdev = other.my_stdev; return *this;}

  virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const { return val * my_stdev; }
  virtual double transformToReportedScale     (const ParameterMgr * mgr, double val) const { return val / my_stdev; }

private:
  double my_stdev;
};

} // End namespace MAXFUN
} // End namespace SAGE

#endif
