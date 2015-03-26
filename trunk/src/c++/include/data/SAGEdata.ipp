//============================================================================
// File:      SAGEdata.ipp
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

namespace SAGE {
namespace APP  {

// ================
// Inline Functions
// ================


inline 
void
SAGE_Data::evaluate_functions()
{
  evaluate_functions(my_pedigrees);
}

inline const ArgumentsFound&  
SAGE_Data::parsed_arguments()
{
  return  my_parsed_arguments;
}

inline 
LSFBase*
SAGE_Data::parameters()
{
  return my_params;
}

inline 
const LSFBase*
SAGE_Data::parameters() const
{
  return my_params;
}

inline
const RPED::RefMultiPedigree&
SAGE_Data::pedigrees() const
{
  return my_pedigrees;
}

inline
RPED::RefMultiPedigree&
SAGE_Data::pedigrees()
{
  //lint --e{1536} returning access to the pedigree is ok
  return my_pedigrees;
}

inline
vector<pair<string, string> >&
SAGE_Data::first_ten_ind()
{
  return my_first_ten_ind;
}

inline
const vector<pair<string, string> >&
SAGE_Data::first_ten_ind() const
{
  return my_first_ten_ind;
}

inline
Output_Streams&
SAGE_Data::get_ostreams() const
{
  return my_output;
}

inline
const MLOCUS::inheritance_model_map&
SAGE_Data::markers() const
{
  return my_pedigrees.info().markers();
}

inline
RPED::genome_description*
SAGE_Data::genome()
{
  return my_genome;
}

inline
const RPED::genome_description*
SAGE_Data::genome() const
{
  return my_genome;
}

inline
LSFBase*
SAGE_Data::regions()
{
  return my_regions;
}

inline
const LSFBase*
SAGE_Data::regions() const
{
  return my_regions;
}

inline
cerrorstream&
SAGE_Data::errors() const
{
  return my_output.errors();
}

inline
std::ostream&
SAGE_Data::info() const
{
  return my_output.info();
}

inline
std::ostream&
SAGE_Data::screen() const
{
  return my_output.screen();
}

inline
std::ostream&
SAGE_Data::messages() const
{
  return my_output.messages();
}

inline
const marker_order&
SAGE_Data::original_marker_order() const
{
  return  my_marker_order;
}

inline
const ArgumentsFound&
SAGE_Data::parsed_arguments() const
{
  return  my_parsed_arguments;
}

} // End namespace APP
} // End namespace SAGE

