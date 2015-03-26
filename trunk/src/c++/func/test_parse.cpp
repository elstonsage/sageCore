//============================================================================
// File:      test_parse.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 2/7/2001.                                                   
//            7/12/2  Modified to allow user to specify missing, binary,
//                    affected and unaffected attributes.  -djb
//                                                                          
// Notes:     Tests the funcparser class.
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================
#include <string>
#include <functional>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "app/SAGEapp.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "output/Output.h"
#include "func/FunctionParser.h"
#include "func/Function.h"

using namespace std;
using namespace SAGE;

int main(int argc, char* argv[])
{
  if (argc != 2)
  {
    cerr << "usage: " << argv[0] << " <parameters> "
         << endl << endl
         << "Command line parameters:"             << endl
         << "  parameters   - Parameter File"  << endl
         << endl << endl;
    exit(EXIT_FAILURE);
  }

  LSFInit();

  sage_cerr << prefix("%%TESTPARSE-%P: ");

  // - Read parameter file.
  //
  cout << "Loading parameters..............";
  LSFBase *params = loadLSFfile(argv[1], "TESTPARSE Parameter file", sage_cerr, false);
  cout << "done." << endl;

  if(!params)
  {
    sage_cerr << priority(fatal) << "Error reading parameter file.... aborting." << endl;
    exit(EXIT_FAILURE);
  }
  
  ofstream info_file;
  info_file.open("testparse.inf");

  if(!info_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testparse.inf.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  ofstream out_file;
  out_file.open("testparse.out");

  if(!out_file)
  {
    sage_cerr << priority(fatal) 
              << "Cannot open output file: testparse.out.  Exiting..." << endl;
    exit(EXIT_FAILURE);
  }

  //print_title(info_file);
  //print_title(out_file);

  SAGE::cerrormultistream errors;
  cerrorstream error_file(info_file);
  error_file.prefix("%%TESTPARSE-%P:");
  errors.insert(sage_cerr);
  errors.restrict(r_ge, error);
  errors.insert(error_file);
  errors.restrict(r_ge, information);
  

  AttrVal                    a;
  SAGE::FUNC::FunctionParser parser;
  SAGE::OUTPUT::Section      s("Parsing information");

  for(LSFList::const_iterator i = params->List()->begin(); i != params->List()->end(); ++i)
  {
    if(!*i) 
      continue;

    if(toUpper((*i)->name()) == "FUNCTION")
    {
      SAGE::FUNC::FunctionParser::TraitData par_data = SAGE::FUNC::FunctionParser::create_trait_data(*i, errors);
      
      SAGE::OUTPUT::Table t("Trait");
      
      t << (SAGE::OUTPUT::TableRow() << "Trait name" << par_data.trait_name)
        << (SAGE::OUTPUT::TableRow() << "Expression" << par_data.expr)
        << (SAGE::OUTPUT::TableRow() << "Time limit" << par_data.time_limit)
        << (SAGE::OUTPUT::TableRow() << "Missing"    << par_data.missing)
        << (SAGE::OUTPUT::TableRow() << "Binary"     << par_data.binary)
        << (SAGE::OUTPUT::TableRow() << "Affected"   << par_data.affected)
        << (SAGE::OUTPUT::TableRow() << "Unaffected" << par_data.unaffected)
        << (SAGE::OUTPUT::TableRow() << "Skip"       << par_data.skip);
      
      const SAGE::FUNC::FunctionParser::TraitData::constant_vector&  constants = par_data.constants;

      for(SAGE::FUNC::FunctionParser::TraitData::constant_vector::const_iterator iter = constants.begin(); iter != constants.end(); ++iter)
      {
        t << (SAGE::OUTPUT::TableRow() << "Constant" << iter->first << iter->second);
      }
      
      s << t;
    }
  }  
  
  out_file << s;
}


