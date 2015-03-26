#ifndef SEGREG_OPTIONS_H
#define SEGREG_OPTIONS_H
//============================================================================
// File:      segreg_options.h
//                                                                          
// Author:    Geoff Wedig (wedig@darwin.cwru.edu)
//                                                                          
// History:   0.1 gcw Initial Implementation                    Aug 2002
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================



namespace SAGE
{
namespace SEGREG
{

/** @name General options for SEGREG
 *  These options control the availability of many of the features in SEGREG.
 *  Eventually, they should all be set to true.
 *
 *  They are contained in the header as constants rather than externs to allow SEGREG
 *  to optimize itself based upon their value.  This means that changing the value would
 *  require a complete rebuild of all SEGREG components.  However, this is desireable, as
 *  the program won't perform unnecessary checks on constant variables.
 */
//@{

const bool  ONSET_AVAILABLE = true;
const bool  FPMM_AVAILABLE = true;
const bool  PREV_CONSTRAINT_AVAILABLE = true;
const bool  PREV_ESTIMATE_AVAILABLE = true;
const bool  TYPE_SUSCEPT_AVAILABLE = true;
const bool  INTERACTIONS_AVAILABLE = true;
const bool  TYPE_CORR_AVAILABLE = false;
const bool  EACH_PEDIGREE_AVAILABLE = false;
const bool  TYPE_PROBABILITIES_AVAILABLE = true;
const bool  PEN_FUNC_OUT_AVAILABLE = true;

//@}

}
}

#endif
