#include "app/SAGEapp.h"
#include "data/SAGEdata.h"

#include <string>
#include <fstream>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "app/SAGEapp.h"
#include "rped/rped.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "util/AutoTrace.h"

namespace SAGE {
namespace APP  {

//=====================================================================
//  class AppData
//=====================================================================
class AppData : public APP::SAGE_Simple_Data
{
  public:

    AppData(const string & program_name, bool debug);
   
    void process_input(int argc, char** argv);

    virtual bool read_analysis();

};

//=======================================================================
//  AppData(...) CONSTRUCTOR
//=======================================================================
AUTOTRACE_CT1(NOARGS,
AppData::AppData(const string& program_name, bool debug) :
	APP::SAGE_Simple_Data(program_name, debug)
{) }

//=======================================================================
//  process_input(...)
//=======================================================================
void
AppData::process_input(int argc, char** argv)
{
  read_parameter_file   (argv[1]);
  read_family_data_file (argv[2], true, true, false, false, true);
  evaluate_functions    ();
  read_analysis         ();
}

//=======================================================================
//  read_analysis()
//=======================================================================
bool
AppData::read_analysis()
{
  return true;
}

//======================================================================
//  class TestApp
//======================================================================
class TestApp : public APP::SAGEapp
{
  public:
    TestApp(int argc=0, char **argv=NULL);

    virtual int  main         ();
    virtual void print_help   (std::ostream &);
};

//======================================================================
//  TestApp(...) CONSTRUCTOR
//======================================================================
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

//======================================================================
//  TestApp::print_help(...)
//======================================================================
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

//======================================================================
//  TestApp::main()
//======================================================================
int TestApp::main()
{
  // 1. Create AppData object:
  AppData data(name, debug());
    
  // 2. Print title:
  print_title(data.info());
  data.process_input(argc, argv); 

  return 0;
}

} // End namespace APP
} // End namespace SAGE

//======================================================================
//
//  MAIN(...)
//
//======================================================================
int main(int argc, char* argv[])
{
  SAGE::APP::TestApp test_app(argc, argv);

  test_app.main();

  return 0;
}
