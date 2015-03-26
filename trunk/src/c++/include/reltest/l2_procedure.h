#ifndef  RELTEST_L2_PROCEDURE_H
#define  RELTEST_L2_PROCEDURE_H

//==========================================================================
//  File:       l2_procedure.h
//
//  Author:     Yeunjoo Song
//
//  History:    Initial Implementation.                              Sep. 01
//              Took out from analysis.                              Jul. 03
//
//  Notes:      This class defines a nonparametric estimation procedure.
//
//  Copyright (c) 2001 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include "maxfun/maxfun.h"
#include "reltest/putative_pair.h"

namespace SAGE
{

namespace RELTEST
{

class L2_error_procedure : public MaxFunction
{
  public :

    L2_error_procedure(bool is_Yj, vector<putative_pair>& p)
    : my_Yj(is_Yj), my_ptt_pairs(p) { nfe = 0; }
    
  protected :

    virtual double evaluate     (vector<double>& theta);
    virtual int    update_bounds(vector<double>& theta);

    bool                     my_Yj;
    vector<putative_pair>&   my_ptt_pairs;  //current putative pairs
};

} // end of namespace RELTEST

} // end of namespace SAGE

#endif
