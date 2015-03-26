#ifndef AO_DATA_H
#define AO_DATA_H
//=====================================================================
//
//  File:	AppData.h
//
//  Author:	Stephen Gross
//
//  Copyright 2001 R. C. Elston
//=====================================================================


#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"
#include "ageon/Datatypes.h"
#include "ageon/Model.h"
#include "ageon/Parser.h"

namespace SAGE {
namespace AO   {

//=====================================================================
//
//  class AppData
//
//=====================================================================
class AppData : public APP::SAGE_Simple_Data
{
  public:
    //=================================================================
    // Constructor:
    //=================================================================

    AppData(const string & program_name, bool debug);

    //=================================================================
    // Public utility functions:
    //=================================================================

    void process_input(int argc, char** argv);

    virtual bool read_analysis();

    //=================================================================
    // Public accessors:
    //=================================================================

    const vector<Model> & analyses();

  private:
    //=================================================================
    // Data members:
    //=================================================================

    vector<Model> my_analyses;
};

inline const vector<Model> & AppData::analyses() { return my_analyses; }

}} // End namespace

#endif
