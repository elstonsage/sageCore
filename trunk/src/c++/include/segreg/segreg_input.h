#ifndef SEGREG_INPUT_H
#define SEGREG_INPUT_H

#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"
#include "segreg/model.h"

namespace SAGE
{

namespace SEGREG
{

class segreg_data : public APP::SAGE_Data
{
  public:

    segreg_data(const string& program_name, bool debug);
    ~segreg_data();
   
    void input(int argc, char** argv);
   
    virtual bool read_analysis();

    const vector<model>& analyses() { return my_analyses; }

    vector<unsigned> locus_indic_vec; // due to JA (trial)
   
    double like_cutoff; // due to JA, see comment in parser.h 

  protected:

     vector<model>  my_analyses;
};

}

} // end of namespace SAGE 

#endif
