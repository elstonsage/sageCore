#include <string>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "numerics/normal_pdf.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "app/SAGEapp.h"
#include "data/SAGEdata.h"
#include "rped/rped.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "mfsubmodels/mfsubmodels.h"
#include "maxfunapi/maxfunapi.h"

namespace SAGE {

//============
// MyFunction
//============

class MyFunction : public MaxFunction
{
public:
  MyFunction();

  virtual double evaluate      (parameter_vector& theta);
  virtual int    update_bounds (parameter_vector& theta);

  MAXFUN::ParameterMgr m;
  MFSUBMODELS::TransformationSubmodel sm;

  vector<double> base_dataset;

  vector<double> working_dataset;
};

MyFunction::MyFunction()
{
  // The following dataset represents a set of values taken from a normal
  // distribution with mean = 0.0, stdev = 0.27386. The values have been
  // subsequently squared, then been added to two.

  // Therefore, this program should estimate lambda2 = -2, and lambda1 = 0.5.

  base_dataset.push_back(2.00);
  base_dataset.push_back(2.00);
  base_dataset.push_back(2.01);
  base_dataset.push_back(2.02);
  base_dataset.push_back(2.04);
  base_dataset.push_back(2.06);
  base_dataset.push_back(2.09);
  base_dataset.push_back(2.12);
  base_dataset.push_back(2.16);
  base_dataset.push_back(2.20);
  base_dataset.push_back(2.00);
  base_dataset.push_back(2.00);
  base_dataset.push_back(2.01);
  base_dataset.push_back(2.02);
  base_dataset.push_back(2.04);
  base_dataset.push_back(2.06);
  base_dataset.push_back(2.09);
  base_dataset.push_back(2.12);
  base_dataset.push_back(2.16);
  base_dataset.push_back(2.20);

  sm.set(MFSUBMODELS::TransformationSubmodel::george_elston, model_input(0.5), model_input(2, true), QNAN, QNAN);
}

double MyFunction::evaluate(parameter_vector& tr)
{
  double v = 1, t = 0;

  for(vector<double>::iterator i = working_dataset.begin(); i != working_dataset.end(); ++i)
  {
    t = NUMERICS::normal_pdf(*i - 2, 0.27386);
    v *= t;
  }

  return v;
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
  int err_code = 0;

  working_dataset = base_dataset;

  if(sm.calculate_geom_mean(working_dataset))
  {
    err_code = sm.transform(working_dataset);
    err_code = 0; // Override error
  }
  else
  {
    err_code = 100;
  }

  return err_code;
}

int testMyFunction()
{
  cout << "Testing transformation maximization stuff..." << endl << endl;

  MyFunction my_func; 

  cout << "Evaluate once: " << MAXFUN::evaluateOnce(my_func.m, my_func).first << endl;

  MAXFUN::SequenceCfg cfg(MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, "Maxtest", "neat function");
  MAXFUN::DebugCfg    dbg(MAXFUN::DebugCfg::COMPLETE);

  my_func.m.addSubModel(&my_func.sm);

  MAXFUN::Results my_results = MAXFUN::maximize(my_func.m, my_func, cfg, dbg);

  my_func.m.removeSubmodels();

  cout << MAXFUN::OutputFormatter::convertEstimates(my_results) << MAXFUN::OutputFormatter::convertMatrix(my_results);

  return 0;
}

void testP1()
{
  // Set up variables:
  std::istringstream input_stream;
  input_stream.str("OPTION=NONE;");
  LSFBase * param1 = new LSFBase ("transformation");
  std::cout << input_stream.str() << std::endl;
  LSF_input(input_stream).input_to(param1, false);
  std::cout << MFSUBMODELS::Transformation::Parser::parse_block(param1, std::cout, sage_cerr).dump();
}

void testP2()
{
  // Set up variables:
  std::istringstream input_stream;
  input_stream.str("OPTION=BOX_COX; lambda1,val=0.3,upper_bound=9.0,lower_bound=-10,fixed=true;");
  LSFBase * param1 = new LSFBase ("transformation");
  std::cout << input_stream.str() << std::endl;
  LSF_input(input_stream).input_to(param1, false);
  std::cout << MFSUBMODELS::Transformation::Parser::parse_block(param1, std::cout, sage_cerr).dump();
}

void testP3()
{
  // Set up variables:
  std::istringstream input_stream;
  input_stream.str("OPTION=GEORGE_ELSTON; lambda2,val=2,fixed=false");
  LSFBase * param1 = new LSFBase ("transformation");
  std::cout << input_stream.str() << std::endl;
  LSF_input(input_stream).input_to(param1, false);
  std::cout << MFSUBMODELS::Transformation::Parser::parse_block(param1, std::cout, sage_cerr).dump();
}

class Max
{
  public:
  
    explicit Max(MFSUBMODELS::Transformation::Facade & _facade) : facade(_facade)
    {
    }

    int go(MAXFUN::ParameterMgr & mgr)
    {
      std::vector<double> v;
      v.push_back(1);
      v.push_back(2);
      v.push_back(3);
      v.push_back(4);
      v.push_back(5);
      
      facade.calculate_geometric_mean(v);
      facade.transform(v);
      
      std::cout << "l1 = " << mgr("Transformation", "Lambda1") << " l2 = " << mgr("Transformation", "Lambda2") << " ";

      for(size_t i = 0; i < v.size(); ++i)
        std::cout << v[i] << " ";

      std::cout << std::endl;
      
      return 0;
    }
    
    MFSUBMODELS::Transformation::Facade & facade;
};

void testM1()
{
  MAXFUN::ParameterMgr m;
  MAXFUN::Function     f(m);

  MFSUBMODELS::Transformation::Configuration transf_config;
  transf_config.set_type(MFSUBMODELS::Transformation::BOX_COX);  

  MFSUBMODELS::Transformation::Facade facade(transf_config, f);
  
  Max max(facade);
  
  f.addStep(boost::bind(&Max::go, &max, _1), "");
  f.setEvaluator(boost::lambda::constant(5.0));
  
  MAXFUN::Maximizer::maximize(f);
}

void test2()
{
  LSFInit();
  
  testP1();
  testP2();
  testP3();

  testM1();
}


} // End namespace SAGE

//==================================================================
//
//                 MAIN IMPLEMENTATION

int main(int argc, char* argv[])
{
  SAGE::testMyFunction();

  SAGE::test2();

  return 0;
}

