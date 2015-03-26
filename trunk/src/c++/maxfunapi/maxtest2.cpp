#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "numerics/normal_pdf.h"
#include "maxfunapi/maxfunapi.h"
#include "boost/function.hpp"
#include "boost/lambda/bind.hpp"

namespace SAGE {

int calculateVar(SAGE::MAXFUN::ParameterMgr & mgr)
{
  double stdev = mgr.getParameter("global", "stdev").getCurrentEstimate();

  mgr.getParameter("global", "var").setCurrentEstimate(stdev * stdev);

  return 0;    
}

double getPdf(SAGE::MAXFUN::ParameterMgr & mgr)
{
  double mean = mgr.getParameter("global", "mean") .getCurrentEstimate(),
         var  = mgr.getParameter("global", "var")  .getCurrentEstimate(),
         x    = mgr.getParameter("global", "x")    .getCurrentEstimate();
         
  return SAGE::NUMERICS::normal_pdf(x, mean, var);
}

void go()
{
  SAGE::MAXFUN::ParameterMgr m;
  SAGE::MAXFUN::Function     f(m);

  m.addParameter("global", "mean",  SAGE::MAXFUN::Parameter::FIXED,        10.0);
  m.addParameter("global", "stdev", SAGE::MAXFUN::Parameter::FIXED,        5.0);
  m.addParameter("global", "var",   SAGE::MAXFUN::Parameter::DEPENDENT,    5.0);
  m.addParameter("global", "x",     SAGE::MAXFUN::Parameter::INDEPENDENT,  0.0);

  m.addParameter("special", "numerator",             SAGE::MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 3.0);
  m.addParameter("special", "denominator",           SAGE::MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 5.0);
  m.addParameter("special", "numerator/denominator", SAGE::MAXFUN::Parameter::DEPENDENT);
  m.addParameter("special", "2.0 / denominator",     SAGE::MAXFUN::Parameter::DEPENDENT);
  m.addParameter("special", "numerator/ 2.0",        SAGE::MAXFUN::Parameter::DEPENDENT);
  m.addParameter("special", "n*d",                   SAGE::MAXFUN::Parameter::DEPENDENT);

  f.addDependentParamStep("special", "numerator/denominator",
    MF_PARAM(m.getParamID("special", "numerator")) / MF_PARAM(m.getParamID("special", "denominator")));
  
  f.addDependentParamStep("special", "2.0 / denominator", 
    MF_CONSTANT(2.0) / MF_PARAM(m.getParamID("special", "denominator")));
  
  f.addDependentParamStep("special", "numerator/ 2.0", 
    MF_PARAM(m.getParamID("special", "numerator")) / MF_CONSTANT(2.0));
  
  f.addDependentParamStep("special", "n*d",
    MF_PARAM(m.getParamID("special", "numerator")) * MF_PARAM(m.getParamID("special", "denominator")));

  f.addStep(boost::bind(calculateVar, _1), "calculate variance");

  f.setEvaluator(boost::bind(getPdf, _1));
  
  SAGE::MAXFUN::Results results =
    SAGE::MAXFUN::Maximizer::maximize(f, SAGE::MAXFUN::SequenceCfg(), SAGE::MAXFUN::DebugCfg(SAGE::MAXFUN::DebugCfg::COMPLETE));

  std::cout << SAGE::MAXFUN::OutputFormatter::convertEstimates (results)
            << SAGE::MAXFUN::OutputFormatter::convertMatrix    (results);
}

} // End namespace SAGE

int main(int argc, char* argv[])
{
  SAGE::go();

  return 0;
}
