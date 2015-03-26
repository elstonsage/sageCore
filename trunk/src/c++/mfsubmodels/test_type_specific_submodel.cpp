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

//=====================================================================
//=====================================================================
//
//                    MAXIMIZATION STUFF
//
//=====================================================================
//=====================================================================

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
private:
  MyFunction(const MyFunction &);
  MyFunction& operator=(const MyFunction &);
};

MyFunction::MyFunction()
{
  MFSUBMODELS::NewTypeSpecificSubmodel & mean_sm = *(static_cast<MFSUBMODELS::NewTypeSpecificSubmodel*>(m.createSubmodel(MFSUBMODELS::new_type_specific_submodel).get()));
  MFSUBMODELS::NewTypeSpecificSubmodel & var_sm  = *(static_cast<MFSUBMODELS::NewTypeSpecificSubmodel*>(m.createSubmodel(MFSUBMODELS::new_type_specific_submodel).get()));

  mean_sm.setCategory       (MFSUBMODELS::NewTypeSpecificSubmodel::MEAN);
  mean_sm.setTypeConstraint (MFSUBMODELS::TypeConstraint::THREE_INC);
  mean_sm.setTypeValues     (0, -3, false, false);
  mean_sm.setTypeValues     (1,  0, false, false);
  mean_sm.setTypeValues     (2,  3, false, false);

  var_sm.setCategory       (MFSUBMODELS::NewTypeSpecificSubmodel::VARIANCE);
  var_sm.setTypeConstraint (MFSUBMODELS::TypeConstraint::ONE);
  var_sm.setTypeValues     (0, 1, true, true);
}

double MyFunction::evaluate(parameter_vector& tr)
{
  const MFSUBMODELS::NewTypeSpecificSubmodel & mean_sm = m.getSubmodel("mean",     MFSUBMODELS::new_type_specific_submodel);
  const MFSUBMODELS::NewTypeSpecificSubmodel & var_sm  = m.getSubmodel("variance", MFSUBMODELS::new_type_specific_submodel);

  double v = NUMERICS::normal_pdf(mean_sm.getValue(0) + 1, var_sm.getValue(0)) +
             NUMERICS::normal_pdf(mean_sm.getValue(1)    , var_sm.getValue(0)) +
             NUMERICS::normal_pdf(mean_sm.getValue(2) - 1, var_sm.getValue(0));

  return v;
} 

int MyFunction::update_bounds(parameter_vector& tr)
{
  return 0;
}

int testMyFunction()
{
  cout << "Testing maximization stuff..." << endl << endl;

  MyFunction my_func; 

  cout << "Evaluate once: " << MAXFUN::evaluateOnce(my_func.m, my_func).first << endl;

  MAXFUN::SequenceCfg cfg(MAXFUN::SequenceCfg::DEFAULT_MAXIMIZATION, "Maxtest", "neat function");
  MAXFUN::DebugCfg    dbg(MAXFUN::DebugCfg::COMPLETE);

  MAXFUN::Results my_results = MAXFUN::maximize(my_func.m, my_func, cfg, dbg);

  cout << MAXFUN::OutputFormatter::convertEstimates(my_results) << MAXFUN::OutputFormatter::convertMatrix(my_results);

  return 0;
}

//===============================================================================================
//===============================================================================================
//===============================================================================================
//
//                                 MEAN / VARIANCE INCOMPATIBILTIY STUFF
//
//===============================================================================================
//===============================================================================================
//===============================================================================================

int testMeanVarIncompatibility()
{
  cout << endl << "Testing mean/variance incompatibility stuff..." << endl << endl;

  MAXFUN::ParameterMgr m;

  MFSUBMODELS::NewTypeSpecificSubmodel & mean_sm = *(static_cast<MFSUBMODELS::NewTypeSpecificSubmodel*>(m.createSubmodel(MFSUBMODELS::new_type_specific_submodel).get()));
  MFSUBMODELS::NewTypeSpecificSubmodel & var_sm  = *(static_cast<MFSUBMODELS::NewTypeSpecificSubmodel*>(m.createSubmodel(MFSUBMODELS::new_type_specific_submodel).get()));

  mean_sm .setCategory       (MFSUBMODELS::NewTypeSpecificSubmodel::MEAN);
  var_sm  .setCategory       (MFSUBMODELS::NewTypeSpecificSubmodel::VARIANCE);
  mean_sm .setTypeConstraint (MFSUBMODELS::TypeConstraint(MFSUBMODELS::TypeConstraint::ONE));
  var_sm  .setTypeConstraint (MFSUBMODELS::TypeConstraint(MFSUBMODELS::TypeConstraint::THREE));

  cerrorstream temp_err(cout);

  MFSUBMODELS::NewTypeSpecificParser::areMeanAndVarianceCompatible(mean_sm, var_sm, temp_err);

  return 0;
}

//===============================================================================================
//===============================================================================================
//===============================================================================================
//
//                                 PARSING STUFF
//
//===============================================================================================
//===============================================================================================
//===============================================================================================

//=========
// TestApp
//=========

class TestApp : public APP::SAGEapp
{
  public:
    TestApp(int argc=0, char **argv=NULL);

    virtual int  main         ();
    virtual void print_help   (std::ostream &);

  private:
    void perform_analyses ();
};

//==========
// TestData
//==========

class TestData : public APP::SAGE_Simple_Data
{
  public:

    TestData(const string & program_name, bool debug);
   
    void process_input(int argc, char** argv);

    virtual bool read_analysis();
};

//========================
// TESTAPP IMPLEMENTATION
//========================

TestApp::TestApp(int argc, char** argv)
     : APP::SAGEapp(APP::APP_AGEON, false, argc, argv)
{
  if(arg_count != 2)
  {
    print_help (cerr);
    exit       (0);
  }

  LSFInit();
}

void TestApp::print_help(ostream & o)
{
  o << "usage: " << argv[0] << " <parameters> <pedigree>" << endl
    <<                                                       endl 
    <<                                                       endl  
    << "Command line parameters:           "              << endl
    << "  parameters   - Parameter File    "              << endl
    << "  pedigree     - Pedigree Data File"              << endl
    <<                                                       endl;
}

int TestApp::main()
{
  perform_analyses();

  return EXIT_SUCCESS;
}

void
TestApp::perform_analyses()
{
  // 1. Create TestData object:
  TestData data(name, debug());
    
  // 2. Print title:
  print_title(data.info());
  data.process_input(argc, argv); 
}

//================================================================
//
//                 TESTDATA IMPLMENTATION

TestData::TestData(const string& program_name, bool debug) :
	APP::SAGE_Simple_Data(program_name, debug)
{}

void
TestData::process_input(int argc, char** argv)
{
  read_parameter_file   (argv[1]);
  read_family_data_file (argv[2], true, false, false, true, true);
  evaluate_functions    ();
  read_analysis         ();
}

bool
TestData::read_analysis()
{
  cerrorstream temp_err(cout);

  for(LSFList::const_iterator i = my_params->List()->begin(); i != my_params->List()->end(); ++i)
  {
    if(!*i) continue;

    string param_name = toUpper((*i)->name());

    if(param_name.substr(0,5) != "TYPE_")
      continue;

    string subblock_name = param_name.substr(5);

    MFSUBMODELS::NewTypeSpecificSubmodel::CategoryEnum c;

         if(subblock_name == "MEAN")    c = MFSUBMODELS::NewTypeSpecificSubmodel::MEAN;
    else if(subblock_name == "VAR")     c = MFSUBMODELS::NewTypeSpecificSubmodel::VARIANCE;
    else if(subblock_name == "SUSCEPT") c = MFSUBMODELS::NewTypeSpecificSubmodel::SUSCEPTIBILITY;
    else continue;
    
    cout <<                                                        endl
         << "=================================================" << endl
         <<                                                        endl 
         << "ENCOUNTERED TYPE_" << subblock_name << " SUBMODEL" << endl
         <<                                                        endl;

    LSFDump(cout).dump(*i);

    MAXFUN::ParameterMgr                   m;
    string                                 sm_name = m.createSubmodel(MFSUBMODELS::new_type_specific_submodel)->getName();
    MFSUBMODELS::NewTypeSpecificSubmodel & sm      = m.getSubmodel(sm_name, MFSUBMODELS::new_type_specific_submodel);

    if(MFSUBMODELS::NewTypeSpecificParser::parse(sm, *i, temp_err))
      sm.dump(cout);
  }    

  return true;
}

int testParsing(int argc, char* argv[])
{
  cout << "Testing Parsing stuff..." << endl << endl;

  TestApp app(argc, argv);

  app.main();

  return 0;
}

} // End namespace SAGE

//==================================================================
//
//                 MAIN IMPLEMENTATION

int main(int argc, char* argv[])
{
  SAGE::testParsing(argc, argv);
  SAGE::testMeanVarIncompatibility();
  SAGE::testMyFunction();

  return 0;
}
