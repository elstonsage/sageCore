//============================================================================
// File:      pedparse.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   Created 11/00                                       
//                                                                          
// Notes:     Inline implementation of the following classes -
//              pedinfo_parser
//                                                                          
// Copyright (c) 2000 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

inline const pedinfo_parser::trait_list&
pedinfo_parser::traits() const
{
  return my_traits;
}
    
inline const std::string&       
pedinfo_parser::file_name() const
{
  return my_file_name;
}
    
inline bool             
pedinfo_parser::show_each_pedigree() const
{
  return my_show_each_pedigree;
}

inline bool             
pedinfo_parser::suppress_general() const
{
  return my_suppress_general;
}
