//============================================================================
// File:      Calculator.ipp
//
// Author:    Dan Baechle
//
// History:   2/28/12 - created.                                   djb
//
// Notes:    
//
//
// Copyright (c) 2012 R.C. Elston
// All Rights Reserved
//============================================================================


//============================================================================
// IMPLEMENTATION:  Calculator       
//============================================================================
//
inline const map<string, FortranMatrix<double> >&   
Calculator::getSharedEffects() const
{
  return  my_shared_effects;
}                             
       
inline const FortranMatrix<double>&   
Calculator::getRandomEffect() const
{
  return  my_random_effect;
}

inline const vector<FPED::MemberConstPointer>&  
Calculator::getMemberLookup() const
{
  return  my_member_lookup;
}
                
inline const map<FPED::MemberConstPointer, size_t>&  
Calculator::getIndexLookup() const
{
  return  my_index_lookup;
}

