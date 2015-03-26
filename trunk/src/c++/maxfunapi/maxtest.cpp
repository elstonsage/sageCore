#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "maxfunapi/maxfunapi.h"

using namespace std;
using namespace SAGE;

// SUBMODEL

class MySubmodel : public MAXFUN::Submodel
{
public:

  double y_intercept;

  virtual int update()
  {
    y_intercept = getParameterMgr()->getParameter("submodel", "y").getCurrentEstimate();
    return 0;
  }

  virtual int finalizeConfiguration()
  {
    my_parameters.push_back(MAXFUN::ParameterInput("submodel", "y", MAXFUN::Parameter::FIXED, 2.0, -1.0, 3.0));

    return 0;
  }
};

// MAXFUNCTION

class SpecialTransformer : public MAXFUN::Transformer
{
public:
  virtual double transformToMaximizationScale (const MAXFUN::ParameterMgr * mgr, double val) const { return val / 2; }
  virtual double transformToReportedScale     (const MAXFUN::ParameterMgr * mgr, double val) const { return val * 2; }
};

class MyFunction : public MaxFunction
{
public:
  MyFunction();

  virtual double evaluate      (parameter_vector& theta);
  virtual int    update_bounds (parameter_vector& theta);

  MAXFUN::ParameterMgr m;
  int m1;
  int m2;
  MySubmodel sm;
};

MyFunction::MyFunction()
{
  m1 = m.addParameter("global", "X1", MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, -1.0);

//  m.getParameter("global", "X1").setIncludeInOutput(false);

//  m.addParameter("tr. global", "X1", MAXFUN::Parameter::DEPENDENT, -1.0, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);

  m2 = m.addParameter("global", "X2", MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 1.0);

//  m.addTransformer(MAXFUN::TransformerShCstPtr(new SpecialTransformer()), "global", "X1", "tr. global", "X1");

  int m3 = m.addParameter("global", "X3", MAXFUN::Parameter::DEPENDENT);

  vector<int> ids;
  ids.push_back(m1);
  ids.push_back(m2);

  m.addParamCalculator(m3, MAXFUN::ParamCalculatorShPtr(new MAXFUN::AdditiveParamCalculator(ids)));
}

double MyFunction::evaluate(parameter_vector& tr)
{
  double v = exp(-((m(m1)-1) * (m(m1)-1))) +
             exp(-(m(m2) * m(m2))) +
             sm.y_intercept;

  return v;
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
  return 0;
}

int main(int argc, char* argv[])
{
  MyFunction my_func; 

  my_func.m.createSubmodel(MAXFUN::SMType<MAXFUN::FooSubmodel> ());

  cout << "Evaluate once: " << MAXFUN::evaluateOnce(my_func.m, my_func).first << endl;

  MAXFUN::SequenceCfg cfg(MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, "Maxtest", "neat function");
  MAXFUN::DebugCfg    dbg(MAXFUN::DebugCfg::COMPLETE);

  MAXFUN::Results my_results = MAXFUN::maximize(my_func.m, my_func, cfg, dbg);

  std::cout << MAXFUN::OutputFormatter::convertEstimates (my_results)
            << MAXFUN::OutputFormatter::convertMatrix    (my_results)
            << MAXFUN::JointTest(my_results, my_results).summarizeAsTable();

  return 0;
}
