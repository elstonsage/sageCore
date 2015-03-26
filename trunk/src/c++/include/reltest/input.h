#ifndef RELTEST_INPUT_H
#define RELTEST_INPUT_H

//==========================================================================
//  File:     input.h
//
//  Author:   Yeunjoo Song
//
//  History:  Initial implementation.                                Jul 03
//
//  Notes:
//
//  Copyright (c) 2003 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "data/SAGEdata.h"

namespace SAGE
{

namespace RELTEST
{

class reltest_data : public APP::SAGE_Data
{
public:

  reltest_data(const string& program_name, bool debug);

  ~reltest_data();

  bool input(int argc, char** argv);

  virtual bool read_analysis();
  
  const vector<LSF_ptr<LSFBase> >& get_analysis() const;

private:

  vector<LSF_ptr<LSFBase> >    my_analysis;

};

// ================
// Inline Functions
// ================

inline
const vector<LSF_ptr<LSFBase> >& reltest_data::get_analysis() const
{
  return my_analysis;
}


} // end of namespace RELTEST

} // end of namespace SAGE

#endif
