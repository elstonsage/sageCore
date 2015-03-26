//============================================================================
// File:      rpedtest_input.ipp
//                                                                          
// Author:    
//                                                                          
// History:   3-27-01 modified to add test of genome_description class.  - djb                                                   
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2001 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

// ================
// Inline Functions
// ================

namespace SAGE {
namespace RPED {


inline SymbolTable* rptest_data::parameters() const
{ return params; }

inline AttrVal rptest_data::parameter (string s) const
{ return params->query(s); }

inline const RefMultiPedigree* rptest_data::pedigrees() const
{ return &pedigree_data; }

} // End namespace RPED
} // End namespace SAGE
