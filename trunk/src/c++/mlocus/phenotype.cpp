//============================================================================
//  File:       phenotype.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#ifdef _MSC_VER
    #include <app/SAGEconfig.h>
    #pragma hdrstop
#endif
#include "mlocus/phenotype.h"
 
namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  IMPLEMENTATION: phenotype_info
//============================================================================
//
phenotype_info::phenotype_info()
  : name(), id(NPOS)
{}


phenotype_info::phenotype_info(const string& n, uint i)
  : name(n), id(i)
{}


bool
phenotype_info::operator ==(const phenotype_info& p) const
{
    return name == p.name;
}


bool
phenotype_info::operator !=(const phenotype_info& p) const
{
    return name != p.name;
}


//============================================================================
//  IMPLEMENTATION: phenotype
//============================================================================
//
phenotype::phenotype()
  : my_info(new phenotype_info())
{}


phenotype::phenotype(const phenotype& p)
  : my_info(p.my_info)
{}


phenotype::phenotype(const string& n, uint i)
  : my_info(new phenotype_info(n, i))
{}


phenotype::~phenotype()
{}


const phenotype&
phenotype::operator =(const phenotype& p)
{
    if (&p != this  &&  p.my_info != my_info)
    {
        my_info = p.my_info;
    }
    return *this;
}

} // End namespace MLOCUS
} // End namespace SAGE

