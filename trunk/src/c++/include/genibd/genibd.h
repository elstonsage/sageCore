#ifndef GENIBD_H
#define GENIBD_H

//============================================================================
// File:      genibd.h
//
// Author:    Yeunjoo Song
//
// History:   Initial implementation                                   Nov. 03
//
// Notes:
//
// Copyright (c) 2003 R.C. Elston
// All Rights Reserved
//============================================================================

#include "app/SAGEapp.h"
#include "genibd/params.h"

namespace SAGE
{

namespace GENIBD
{

//----------------------------------------------------------------------------
//  Class:    genibd
//                                                                          
//  Purpose:  Application class for the GENIBD program.
//                                                                          
//----------------------------------------------------------------------------
//
class genibd : public APP::SAGEapp
{
  public:

    // Constructor/destructor.
    genibd(int argc = 0, char** argv = NULL);
    ~genibd();
    
    // Run the application.
    virtual int main();

  protected:
  
    LSF_ptr<genome_description> build_genome_description(const RefMPedInfo& minfo,
                                                         LSFBase*           regions,
                                                         bool               multipoint,
                                                         double             distance,
                                                         cerrorstream&      err = sage_cerr);

    void print_analysis_table(genibd_parameters& params, ostream& o);

};

} // end of namespace GENIBD

} // end of namespace SAGE

#endif
