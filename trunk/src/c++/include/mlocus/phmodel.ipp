//============================================================================
//  File:       phmodel.ipp
//
//  Author:     Bob Steagall
//
//  History:    Version 1.00.00
//
//  Notes:    
//
//  Copyright (c) 1999 R.C. Elston
//  All Rights Reserved

#ifndef PHMODEL_H
#include "mlocus/phmodel.h"
#endif

namespace SAGE   {
namespace MLOCUS {

//============================================================================
//  IMPLEMENTATION: phenotype_model
//============================================================================
//
inline uint
phenotype_model::phenotype_count() const
{
    return (my_info->phenotypes.size() - 1);
}

inline uint
phenotype_model::alias_count() const
{
    return (my_info->phenotype_names.size() - phenotype_count()) - 2;
}


//----------
//
inline phenotype_model::phenotype_iterator
phenotype_model::phenotype_begin() const
{
    return my_info->phenotypes.begin() + 1;
}

inline phenotype_model::phenotype_iterator
phenotype_model::phenotype_end() const
{
    return my_info->phenotypes.end();
}

inline const phenotype&
phenotype_model::get_phenotype(uint id) const
{
  if(id == my_info->missing_ptid)
    return get_missing_phenotype();
  else
    return my_info->phenotypes[id];
}

//----------
//
inline uint
phenotype_model::get_missing_phenotype_id() const
{
    return my_info->missing_ptid;
}

//----------
//
inline const string&
phenotype_model::name() const
{
    return my_info->name;
}

inline const string&
phenotype_model::missing_phenotype_name() const
{
    return my_info->missing_ptname;
}

/// Returns \c true if there are phenotypes which we generated
/// from genotype names, \c false otherwise
inline bool
phenotype_model::has_generated_phenotypes() const
{
  return my_info->my_has_generated_phenotypes;
}

/// Set to \c true if there are phenotypes which were created
/// from external sources (add_phenotype()), \c false otherwise.
inline bool
phenotype_model::has_external_phenotypes() const
{
  return my_info->my_has_external_phenotypes;
}

inline void
phenotype_model::set_name(const string& nm)
{
  uniquify();
  my_info->name = nm;
}

} // End namespace MLOCUS
} // End namespace SAGE

