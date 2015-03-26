
//======================================================================
//
//
//  File:	ageon.cpp
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================

#include "ageon/ageon.h"

namespace SAGE {
namespace AO   {

//=======================================================================
//
//  ConvertAnalysisType(...)
//
//=======================================================================
string ConvertAnalysisType(int t, bool pooled)
{
  string x = "";

  //x += SusceptibilitiesEqual (t) ? "susceptibilities equal, " : "susceptibilities free, ";

  if( SusceptibilitiesEqual(t) )
    x += "susceptibilities equal, ";
  else if( pooled )
    x += "susceptibilities as pooled, ";
  else
    x += "susceptibilities free, ";

  x += NoTruncation          (t) ? "no truncation"            : "using truncation";

  return x;
}

//======================================================================
//
//  AgeonApp(...) CONSTRUCTOR
//
//======================================================================
AgeonApp::AgeonApp(int argc, char** argv)
  : APP::SAGEapp(APP::APP_AGEON, true, argc, argv)
{
  if((arg_count != 2) && (arg_count != 4))
  {
    exit       (0);
  }

  LSFInit();
}


//======================================================================
//
//  AgeonApp::main()
//
//======================================================================
int AgeonApp::main()
{
  perform_analyses();

  print_inf_banner(cout);

  return EXIT_SUCCESS;
}

//======================================================================
//
//  AgeonApp::perform_analyses(...)
//
//======================================================================
void
AgeonApp::perform_analyses()
{
  // 1. Create AppData object:
  AppData data(name, debug());

  // 2. Print title:
  print_title(data.info());
  data.process_input(argc, argv);

  // 4. Verify that there are valid analyses:

  if(!data.analyses().size())
  {
    data.errors() << priority(fatal)
                      << "No valid analyses specified.  Program cannot continue."
                      << endl;

    exit(EXIT_FAILURE);
  }

  // 6. Create the analyses and execute them:

  for(vector<Model>::const_iterator model_itr = data.analyses().begin(); model_itr != data.analyses().end(); ++model_itr)
  {
    FPED::FilteredMultipedigree f(data.pedigrees());

    FPED::MPFilterer::add_multipedigree(f, data.pedigrees());

    f.construct();

    SAMPLING::PartitionedMemberDataSample sample(f, data.errors());

    (const_cast<Model &>(*model_itr)).setupSample(sample);

    AnalysisOutput output(sample);

    AnalysisWrapper wrapper1(*model_itr, data.pedigrees(), sample, NO_TRUNCATION  | SUSCEPTIBILITIES_EQUAL, output),
                    wrapper2(*model_itr, data.pedigrees(), sample, NO_TRUNCATION  | SUSCEPTIBILITIES_FREE,  output),
                    wrapper3(*model_itr, data.pedigrees(), sample, USE_TRUNCATION | SUSCEPTIBILITIES_EQUAL, output),
                    wrapper4(*model_itr, data.pedigrees(), sample, USE_TRUNCATION | SUSCEPTIBILITIES_FREE,  output);

    string ofilename = output.get_models()[0].get_ofilename();

    ofstream ofile_sum (string(ofilename + ".sum").c_str());
    ofstream ofile_det (string(ofilename + ".det").c_str());

    ofile_sum << getReleaseString() << endl;
    ofile_det << getReleaseString() << endl;

    ofile_sum << "Remember you have agreed to add an appropriate statement (including the" << endl
              << "NIH grant number) under \"acknowledgments\" in any publication of results" << endl
              << "obtained by using this program. Suggested wording is:" << endl << endl
              << "\"(Some of)The results of this paper were obtained by using the software" << endl
              << "package S.A.G.E., which was supported by a U.S. Public Health Service" << endl
              << "Resource Grant (RR03655) from the National Center for Research Resources.\"" << endl
              << endl;

    ofile_det << "Remember you have agreed to add an appropriate statement (including the" << endl
              << "NIH grant number) under \"acknowledgments\" in any publication of results" << endl
              << "obtained by using this program. Suggested wording is:" << endl << endl
              << "\"(Some of)The results of this paper were obtained by using the software" << endl
              << "package S.A.G.E., which was supported by a U.S. Public Health Service" << endl
              << "Resource Grant (RR03655) from the National Center for Research Resources.\"" << endl
              << endl;

    ofile_sum << output.generate_output(EXCLUDE_DETAILED);
    ofile_det << output.generate_output(INCLUDE_DETAILED);

    ExtraOutput::generateFile(sample, output.get_model_traits(), ofilename, data.info());
  }
}

} // End namespace AO
} // End namespace SAGE

//======================================================================
//
//  MAIN(...)
//
//======================================================================
int main(int argc, char* argv[])
{
  SAGE::AO::AgeonApp ageon_app(argc, argv);

  ageon_app.main();

  return 0;
}
