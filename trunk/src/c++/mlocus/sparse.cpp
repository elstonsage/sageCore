//=============================================================================
//  File:       "sparse.cpp"
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
    #pragma warning (disable : 4661)
#endif

#include <memory.h>
#include <assert.h>
#include <sparse.h>

namespace SAGE   {
namespace MLOCUS {

template class sparse_matrix<double>;

void f()
{
    sparse_matrix<double>   smd;
}

} // End namespace MLOCUS
} // End namespace SAGE

