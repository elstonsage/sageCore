#ifndef __TEST_IVG_NEW_INPUT_H
#define __TEST_IVG_NEW_INPUT_H

//
//  test inheritance vector generator input routines
//
//  Author: Yeunjoo Song
//
//  History:  Initial Interface design  Mar. 2002
//
//  Copyright (c) 2002 R. C. Elston

#include <iostream>
#include <fstream>
#include "LSF/Attr.h"
#include "LSF/LSFsymbol.h"
#include "data/SAGEdata.h"
#include "error/errorstream.h"
#include "rped/genome_description.h"


namespace SAGE
{

class test_ivg_data : public APP::SAGE_Data
{
public:

  typedef RPED::genome_description::region_type region;

  test_ivg_data(const string& program_name, bool debug);

  ~test_ivg_data();

  bool input(int argc, char** argv);

  virtual bool read_analysis();
  
  const vector<LSF_ptr<LSFBase> >& analysis() const;

  const region get_region() const;

private:

  void build_genome();

  RPED::genome_description*        my_genome;

  vector<LSF_ptr<LSFBase> >  my_analysis;
};


// ================
// Inline Functions
// ================

inline
const vector<LSF_ptr<LSFBase> >& test_ivg_data::analysis() const
{ 
  return my_analysis;
}

inline
const test_ivg_data::region test_ivg_data::get_region() const
{
  return my_genome->region("REGION");
}

} // end of namespace

#endif
