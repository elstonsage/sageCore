//============================================================================
// File:      data.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   5/30/02 created        -djb
//                                                                          
// Notes:     Inline implementation for class, data.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef MLOD_DATA_H
#include "mlod/data.h"
#endif

//============================================================================
// IMPLEMENTATION:  data
//============================================================================
//
namespace SAGE
{
namespace MLOD
{
  
inline
Data::Data(const string& program_name, bool debugging_active)
    : SAGE_Data(program_name, debugging_active)
{ 
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PARAMETER_FILE,   APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::PEDIGREE_FILE,    APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::TRAIT_MODEL_FILE, APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::LOCUS_FILE,       APP::ArgumentRuleset::ONE));
  my_cmdline_rules.add_rule(APP::ArgumentRuleset::Rule(APP::GENOME_FILE,      APP::ArgumentRuleset::ONE));
}

inline
Data::~Data()
{ }

inline const vector<AnalysisParameters>&
Data::get_analyses() const
{
  return my_analyses;
}

inline const MLOCUS::inheritance_model_map&
Data::get_marker_loci() const
{
  return my_marker_loci;
}

inline const MLOCUS::inheritance_model_map&
Data::get_trait_loci() const
{
  return my_trait_loci;
}

/*
inline const Data::region
Data::get_trait_region() const
{
  return my_trait_genome->region((size_t) 0); 
}
*/

}
}

