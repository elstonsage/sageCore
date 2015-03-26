#include <string>
#include <cassert>
#include <functional>
#include <vector>
#include <iomanip>
#include "LSF/LSFfile.h"
#include "LSF/LSFinit.h"
#include "LSF/LSFsymbol.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/bufferederrorstream.h"
#include "rped/rped.h"
#include "mped/mp.h"
#include "mped/mp_utilities.h"
#include "mlocus/mfile.h"
#include "segreg/test_segreg.h"
#include "segreg/segreg_utilities.h"
#include "segreg/regressive_peeler.h"
#include "segreg/polygenic_penetrance_calculator.h"
#include "segreg/SL_calculator.h"
#include "segreg/member_calculator.h"
#include "segreg/parser.h"
#include "segreg/segreg_calculator.h"

using namespace std;
using namespace SAGE;
using namespace SEGREG;

#define DEBUG(x)

void
  test_segreg_components
    (const model&            mod,
     RPED::RefMultiPedigree& rped,
     bool                    ascer,
     APP::Output_Streams&    output);
      
SEGREGTEST::SEGREGTEST(int argc, char **argv)
  : APP::SAGEapp(APP::APP_SEGREG, true, argc, argv) 
{
  LSFInit();
}

int SEGREGTEST::main()
{
  segreg_data data(name, debug());

  print_title(data.info());

  cerrorstream& errors = data.errors();

  data.input(argc, argv);

  // Run the tests

  if(!data.analyses().size())
  {
    errors << priority(fatal) << "No valid analyses specified.  Program cannot continue." 
           << endl;

    exit(EXIT_FAILURE);
  }
//  else if(!evaluate_pedigree_validity(data))
//  {
//    errors << priority(fatal) << "No valid pedigrees specified. Program cannot continue."
//           << endl;
//
//    exit(EXIT_FAILURE);
//  }
  else
  {
    const vector<model>& analyses = data.analyses();

    for(vector<model>::const_iterator i = analyses.begin(); i != analyses.end(); ++i)
    {
      cout << "RUNNING MODEL" << endl;
      cout << "====================================================" << endl << endl;

      test_segreg_components(*i, data.pedigrees(), false, data.get_ostreams());

      if(i->ascer_sub_model.s_option() != ascertainment_sub_model::none)
        test_segreg_components(*i, data.pedigrees(), true, data.get_ostreams());

      cout << endl << endl
           << "====================================================" << endl << endl;
    }
  }

  cout << endl;


  return 0;
}

using namespace SAGE;

int main(int argc, char **argv)
{
  free(malloc(1));

  SEGREGTEST segreg_test(argc, argv);

  segreg_test.main();

  return 0;
}

void
  test_segreg_components
    (const model&            mod,
     RPED::RefMultiPedigree& rped,
     bool                    ascer,
     APP::Output_Streams&    output)
{
  PedigreeDataSet ped_data(rped, mod, output);
  
  if(!ped_data.is_valid()) return;
  
  cout << "Model has sex effect: " << mod.has_sex_effect() << endl;

  const FPED::Multipedigree& fped = *ped_data.get_raw_data();
  
  model_class  m = mod.get_model_class();

  switch(m)
  {
    case model_A    :
      SEGREG::segreg_utilities::test_continuous_MC(fped,mod, ascer);
      SEGREG::segreg_utilities::test_PENETRANCE(fped,mod,ascer);
      SEGREG::segreg_utilities::test_Regressive(fped,mod);

      break;

    case model_D    :
      SEGREG::segreg_utilities::test_continuous_MC(fped,mod, ascer);
      SEGREG::segreg_utilities::test_PENETRANCE(fped,mod,ascer);
      SEGREG::segreg_utilities::test_Regressive(fped,mod);

      break;

    case model_MLM  :
      SEGREG::segreg_utilities::test_binary_MC       (fped,mod, ascer);
      SEGREG::segreg_utilities::test_binary_fam_resid(fped,mod, ascer);
      SEGREG::segreg_utilities::test_MLM             (fped,mod, ascer);
      SEGREG::segreg_utilities::test_mlm_pairs       (ped_data,mod);
      SEGREG::segreg_utilities::test_mlm_correlations(ped_data,mod);

      SEGREG::segreg_utilities::test_prev_calculation(mod);
      //cout<<"run MLM\n";

      break;

    case model_FPMM :
    {
      primary_type t = mod.get_primary_trait_type();

      switch(t)
      {
        case pt_CONTINUOUS :
          SEGREG::segreg_utilities::test_continuous_MC(fped,mod,ascer);
          break;
        
        case pt_BINARY :
          SEGREG::segreg_utilities::test_binary_MC(fped,mod,ascer);
          SEGREG::segreg_utilities::test_prev_calculation(mod);
          break;
        
        case pt_ONSET :
          SEGREG::segreg_utilities::test_onset_MC(fped,mod,ascer);
          SEGREG::segreg_utilities::test_prev_calculation(mod);
          break;
        
        case pt_NONE:
          throw std::exception();
      }

      SEGREG::segreg_utilities::test_polygenic_penetrance(fped,mod,ascer);
      break;
    }
       
    case model_INVALID :
    
      throw std::exception();

    //  SEGREG::segreg_utilities::test_POLYGENIC_TRANSITION(fped,mod);
    //  SEGREG::segreg_utilities::test_FPMM(fped,mod);

      break;
  }

//  SEGREG::segreg_utilities::test_ASCER_CALC(fped,mod);

//  SEGREG::segreg_utilities::test_SEGREG_CALCULATOR(fped,mod);
}

