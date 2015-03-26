#include "LSF/LSFinit.h"
#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "mlod/mlod.h"
#include "mlod/lod_table.h"
#include "mlod/analysis.h"

namespace SAGE
{
namespace MLOD
{

int MLOD::main()
{
  init();

  my_data.input(argc, argv);
  
  // Test for analyses
  if(my_data.get_analyses().size() != 0)
  {
    run_analyses();

    //print_inf_banner(my_data.get_ostreams().messages());
    print_inf_banner(cout);

    return 0;
  }
  else
  {
    my_data.errors() << priority(critical) << "No valid analyses detected! "
                     << "Program aborting." << endl;
    return 1;
  }

  return 0;
}

void MLOD::init()
{
  LSFInit();
}

void MLOD::run_analyses()
{
  cout << "Processing Analyses...." << endl;


  Analyzer a(my_data.pedigrees(), *my_data.genome(), my_data.get_ostreams());

  for(vector<AnalysisParameters>::const_iterator i = my_data.get_analyses().begin();
      i != my_data.get_analyses().end(); ++i)
  {
    AnalysisData d = a.run_analysis(*i);
  }
}

}
}

int main(int argc, char **argv)
{
  SAGE::MLOD::MLOD *mlod = new SAGE::MLOD::MLOD(argc,argv);
  
  assert(mlod != NULL);

  int r = mlod->main();

  delete mlod;

  return r;
}

