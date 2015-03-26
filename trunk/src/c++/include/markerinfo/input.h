#ifndef __MARKERINFO_NEW_INPUT_H
#define __MARKERINFO_NEW_INPUT_H

//============================================================================
//  File:    input.h
//
//  Author:  Yeunjoo Song
//
//  History: Initial Interface design                                Mar. 2002
//
//  Notes:   MARKERINFO input routines
//
//  Copyright (c) 2002 R. C. Elston
//    All Rights Reserved
//============================================================================

#include <iostream>
#include <fstream>
#include "LSF/Attr.h"
#include "LSF/LSFsymbol.h"
#include "data/SAGEdata.h"
#include "error/errorstream.h"
#include "rped/genome_description.h"


namespace SAGE
{

class markerinfo_data : public SAGE::APP::SAGE_Data
{
public:

  typedef RPED::genome_description::region_type region;

  markerinfo_data(const string& program_name, bool debug);

  ~markerinfo_data();

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
const vector<LSF_ptr<LSFBase> >& markerinfo_data::analysis() const
{ 
  return my_analysis;
}

inline
const markerinfo_data::region markerinfo_data::get_region() const
{
  return my_genome->region("REGION");
}

} // end of namespace

#endif
