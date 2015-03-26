#ifndef LODLINK_INPUT_H
#define LODLINK_INPUT_H
//============================================================================
// File:      input.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 9/20/2                                            
//                                                                          
// Notes:     defines class, lodlink_data.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <limits>
#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "data/SAGEdata.h"

namespace SAGE
{

namespace LODLINK
{

class lodlink_data : public APP::SAGE_Data
{
  public:

    lodlink_data(const string& program_name, bool debug);
   
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
