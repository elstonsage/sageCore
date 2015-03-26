#ifndef APP_DATA_H
#define APP_DATA_H
//======================================================================
///
///  File:	AppData.h
///
///  Author:	Stephen Gross
///
///  Copyright 2001 R. C. Elston
//======================================================================

#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"
#include "assoc/Datatypes.h"
#include "assoc/Model.h"
#include "assoc/Parser.h"

namespace SAGE  {
namespace ASSOC {

class AppData : public APP::SAGE_Simple_Data
{
  public:
    AppData(const string& program_name, bool debug);
   
    vector<Configuration>&  getAnalyses()  { return my_analyses; }
    
    void process_input(int argc, char** argv);
    virtual bool read_analysis();    

  private:
    // Data members
    vector<Configuration>  my_analyses;
};

} 
} 

#endif
