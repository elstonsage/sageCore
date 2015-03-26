//============================================================================
//  File:       phenotype.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved

#ifndef PHENOTYPE_H
#include "mlocus/phenotype.h"
#endif

namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  IMPLEMENTATION: phenotype
//============================================================================
//
inline bool
phenotype::operator ==(const phenotype& rhs) const
{
    return my_info == rhs.my_info  ||  *my_info == *rhs.my_info;
}

inline bool
phenotype::operator !=(const phenotype& rhs) const
{
    return my_info != rhs.my_info  &&  *my_info != *rhs.my_info;
}

//----------
//
inline const string&
phenotype::name() const
{
    return my_info->name;
}

inline uint
phenotype::id() const
{
    return my_info->id;
}

inline 
phenotype::phenotype(phenotype_info* info)
  : my_info(info)
{ }

//----------
//
//lint -e{1762}
inline void
phenotype::set_id(uint i)
{
    my_info->id = i;
}

inline void
phenotype::uniquify()
{
    if (!my_info.unique())
    {
        my_info = boost::shared_ptr<phenotype_info>
                     (new phenotype_info(*my_info));
    }
}

} // End namespace MLOCUS
} // End namespace SAGE

