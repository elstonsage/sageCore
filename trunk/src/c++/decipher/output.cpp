//============================================================================
// File:      output.cpp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   created 11/15/4                               djb
//                                                                          
// Notes:     
//                                                                          
// Copyright (c) 2004 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "decipher/output.h"

using namespace std;

namespace SAGE
{

namespace DECIPHER
{

bool  
order_by_sub_pop::operator ()(const base_em_phenotype_map* map1, 
                                 const base_em_phenotype_map* map2 )
{
  return  map1->sub_pop_name() < map2->sub_pop_name();
}

  

//============================================================================
// IMPLEMENTATION:  writer
//============================================================================
//


}
}
