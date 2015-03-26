#ifndef DECIPHER_INPUT_H
#define DECIPHER_INPUT_H

//============================================================================
// File:      input.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 2/12/2                                            
//                                                                          
// Notes:     defines class, decipher_data.
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"

namespace SAGE
{

namespace DECIPHER
{

using namespace APP;

class decipher_data : public SAGE_Data
{
  public:

    decipher_data(const string& program_name, bool debug);
   
    const vector<LSF_ptr<LSFBase> >& analyses() { return my_analyses; }
    size_t  true_marker_count() const;

    void input(int argc, char** argv);
    virtual bool read_analysis();
    

  private:
    vector<LSF_ptr<LSFBase> >  my_analyses;
    size_t  my_true_marker_count;
};

} 
}

#endif
