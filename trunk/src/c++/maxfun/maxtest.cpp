#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "maxfun/MaxfunAPI.h"

using namespace std;
using namespace SAGE;

// TRANSLATOR

class mytransformer : public MAXFUN::ParameterTransformer
{
  virtual ParameterTransformer* clone() { ParameterTransformer * tmp = new mytransformer(*this); return tmp; }

  virtual MaxfunParameter transformToInternal (const MaxfunParameter & param) 
  { 
    MaxfunParameter newparam = param;

    newparam.setInitialEstimate(newparam.InitialEstimate() * 0.5);

    return newparam; 
  }

  virtual MaxfunParameter transformToExternal (const MaxfunParameter & param)
  { 
    MaxfunParameter newparam = param;

    newparam.FinalEstimate() *= 2.0;

    return newparam; 
  }
};

// SUBMODEL

class MySubmodel : public MAXFUN::MaxfunSubmodel
{
public:

  double y_intercept;

  virtual int update()
  {
    y_intercept = getMaxfunInfo()->GetParameter("submodel", "y").CurrentEstimate();
    return 0;
  }

  virtual int finalizeConfiguration()
  {
    my_parameters.push_back(MAXFUN::ParameterInput("submodel", "y", MAXFUN::type_fixed, 2.0, -1.0, 3.0));

    return 0;
  }

  virtual void setAdvancedParameterOptions()
  {
    mytransformer tr;
    getMaxfunInfo()->GetParameter("submodel", "y").addTransformer(tr);
  }
};

// MAXFUNCTION

class MyFunction : public MaxFunction
{
public:
  MyFunction();

  virtual double evaluate      (parameter_vector& theta);
  virtual int    update_bounds (parameter_vector& theta);

  MAXFUN::MaxfunInfo m;
  int m1;
  MySubmodel sm;
};

MyFunction::MyFunction()
{
  m1 = m.AddParameter("global", "X1", MAXFUN::type_independent, -1.0, -d_INFINITY, d_INFINITY);
}

double MyFunction::evaluate(parameter_vector& tr)
{
  double v = exp(-(m(m1) * m(m1))) + sm.y_intercept;

  return v;
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
  return 0;
}

int main(int argc, char* argv[])
{
  MyFunction my_func; 

  MaxfunConfig cfg(MAXFUN::max_default, "Maxtest");
  MaxfunDebug dbg(MAXFUN::debug_complete);

  my_func.m.AddSubModel(&my_func.sm);

  MaxfunResults my_results = Maximize(my_func.m, my_func, cfg, dbg);
//  MaxfunResults my_results = MaximizeDefault(my_func.m, my_func);

  my_func.m.removeSubmodels();

  cout << FormatEstimates(my_results) << FormatMatrix(my_results);

//  SAGEXML::Document doc("test");
//  MaxfunOutput::addEstimates(my_results, doc);
//  doc.generateOutput();

  return 0;
}
