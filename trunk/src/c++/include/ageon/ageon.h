#ifndef AO_AGEONSET_H
#define AO_AGEONSET_H
//======================================================================
//
//  File:	ageon.h
//
//  Author:	Stephen Gross
//
//  Copyright 2002 R. C. Elston
//======================================================================


#include <string>
#include <fstream>
#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "app/SAGEapp.h"
#include "rped/rped.h"
#include "sampling/sampling.h"
#include "mped/mp.h"
#include "mped/sp.h"
#include "ageon/AnalysisOutput.h"
#include "ageon/AnalysisWrapper.h"
#include "ageon/Calculator.h"
#include "ageon/AppData.h"
#include "ageon/Model.h"
#include "ageon/Validator.h"
#include "ageon/ExtraOutput.h"

namespace SAGE {
namespace AO   {

//======================================================================
//
//  class AgeonApp
//
//======================================================================
class AgeonApp : public APP::SAGEapp
{
  public:
    AgeonApp(int argc=0, char **argv=NULL);

    virtual int  main         ();

  private:
    void perform_analyses ();
};

} // End namespace AO
} // End namespace SAGE

#endif
