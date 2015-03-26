#ifndef PEDPARSER_H
#define PEDPARSER_H
//============================================================================
// File:      pedparse.h  
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 11/00                                        
//                                                                          
// Notes:     Declares the following classes -
//              pedinfo_parser
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <iostream>
#include <string>
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "LSF/LSF.h"
#include "rped/rped.h"

namespace SAGE
{

//----------------------------------------------------------------------------
//  Class:    pedinfo_parser
//                                                                          
//  Purpose:  Defines functions to parse pedinfo_analysis parameters.
//                                                                          
//----------------------------------------------------------------------------
//
class pedinfo_parser
{
  public:
    typedef std::vector<size_t>    trait_list;
  
    // Constructor/destructor.
    pedinfo_parser();
    
    // Gets.
    const trait_list&   traits()               const;
    const std::string&  file_name()            const;
    bool                show_each_pedigree()   const;
    bool                suppress_general()     const;
    
    void                parse(const LSFBase* analysis, cerrormultistream& errors,
                             const RPED::RefMultiPedigree* mp);
    
  private:
    void clear();
  
    // Data members.
    trait_list          my_traits;
    bool                my_show_each_pedigree;
    bool                my_suppress_general;
    std::string         my_file_name;
  
};

#include "pedinfo/pedparse.ipp"
}

#endif
