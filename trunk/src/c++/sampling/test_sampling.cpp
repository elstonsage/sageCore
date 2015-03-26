#include <string>
#include <fstream>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "app/SAGEapp.h"
#include "data/SAGEdata.h"
#include "rped/rped.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "sampling/sampling.h"

namespace SAGE {

//=====================================================================
//
//            CLASS DECLARATIONS


class TestApp : public APP::SAGEapp
{
  public:
    TestApp(int argc=0, char **argv=NULL);

    virtual int  main();
    virtual void print_help(std::ostream&);

  private:
    void perform_analyses();
};

class TestData : public APP::SAGE_Simple_Data
{
  public:

    TestData(const string& program_name, bool debug);
   
    void process_input(int argc, char** argv);

    virtual bool read_analysis() { return true; }
};

//=================================================================================
//
//               TESTAPP IMPLEMENTATION


TestApp::TestApp(int argc, char** argv)
     : APP::SAGEapp(APP::APP_AGEON, false, argc, argv)
{
  if(arg_count != 2)
  {
    print_help(cerr);
    exit(0);
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
  // Create TestData object:
  TestData data(name, debug());
    
  // Print title:
  print_title(data.info());
  data.process_input(argc, argv); 

  FPED::FilteredMultipedigree f(data.pedigrees());
  
  FPED::MPFilterer::add_multipedigree(f, data.pedigrees());
  
  f.construct();
     
  // Create the analyses and execute them:
  SAMPLING::PartitionedMemberDataSample sample(f, data.errors());

  // Import the traits:
  sample.importField("sqrtdbh", "Foo", "sqrtdbh", SAMPLING::Field::ALLOW_AVERAGING);
  sample.importField("ABO",     "Foo", "ABO",     SAMPLING::Field::ALLOW_AVERAGING | 
                                                  SAMPLING::Field::MEAN_ADJUST);
  sample.importField("comt",    "Foo", "comt",    SAMPLING::Field::ALLOW_AVERAGING | 
                                                  SAMPLING::Field::MEAN_ADJUST     | 
                                                  SAMPLING::Field::STDEV_ADJUST);
  sample.importField("KELL",    "Foo", "KELL",    SAMPLING::Field::MEAN_ADJUST);
  sample.importField("P",       "Foo", "P",       SAMPLING::Field::MEAN_ADJUST     | 
                                                  SAMPLING::Field::STDEV_ADJUST);
  sample.importField("FY",      "Foo", "FY");

  sample.finalizeData();

  sample.finalizeUserCreatedData();

  // Dump summary:
  std::cout << sample.getSummaryTable();
  
  /*
  // - Added 6-12-7. djb
  //
  cout << "Valid singletons:  " << sample.getValidSingletonCount() << endl;
  cout << "Invalid singletons:  " << sample.getInvalidSingletonCount() << endl;
  */

  // Dump the sample's contents:
  sample.dumpTraitValues();
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
  read_parameter_file(argv[1]);
  read_family_data_file(argv[2], true, false, false, true, true);
  evaluate_functions();
  read_analysis();
}

} // End namespace AO

//==================================================================
//
//                 MAIN IMPLEMENTATION

int main(int argc, char* argv[])
{
  SAGE::TestApp app(argc, argv);

  app.main();

  return 0;
}
