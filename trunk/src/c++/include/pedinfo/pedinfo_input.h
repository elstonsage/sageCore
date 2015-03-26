#ifndef PEDINFO_INPUT_H
#define PEDINFO_INPUT_H

#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"

namespace SAGE
{

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~ Class:     pedinfo_data                                                 ~
// ~                                                                         ~
// ~ Purpose:   Defines pedinfo input data object derived from SAGEdata.     ~
// ~                                                                         ~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

class pedinfo_data : public APP::SAGE_Simple_Data
{
  public:

    pedinfo_data(const string& program_name, bool debug);
    ~pedinfo_data();
   
    void input(int argc, char** argv);
   
    virtual bool read_analysis();

    const vector<LSF_ptr<LSFBase> >& analyses() { return my_analyses; }

  protected:

     vector<LSF_ptr<LSFBase> >  my_analyses;
};

} // end of namespace SAGE 

#endif
