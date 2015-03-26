#ifndef SEGREG_DEFINITIONS_H
#define SEGREG_DEFINITIONS_H
//============================================================================
// File:      definitions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   3/8/2 created         -djb
//                                                                          
// Notes:     For definitions common to model and sub-models.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include <vector>

namespace SAGE
{

namespace SEGREG
{

using std::vector;

enum model_class { model_A, model_D, model_FPMM, model_MLM, model_INVALID };
enum primary_type { pt_NONE, pt_CONTINUOUS, pt_BINARY, pt_ONSET };

/*
template<class E> 
inline bool  has_element(const vector<E>& vect, const E& element)
{
  vector<E>::const_iterator  begin = vect.begin();
  vector<E>::const_iterator  end   = vect.end();
  vector<E>::const_iterator  result;
  
  result = std::find_if<vector<E>::const_iterator, E>(begin, end, element);
  return result != vect.end();
}
*/
}
}

#endif

