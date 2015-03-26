//=======================================================================
//  File:	MemberCovariateCalculator.h
//
//  Author:	Stephen Gross
//
//  Hisory: 12/7/9  moved inline functions from MemberCovariateCalculator.h
//                  -djb
//
//  Copyright 2002 R. C. Elston
//
//=======================================================================


//======================================================================
// IMPLEMENTATION:  MemberCovariateCalculator
//======================================================================
//
inline
MemberCovariateCalculator::~MemberCovariateCalculator()
{}

inline double 
MemberCovariateCalculator::getZScore (size_t id) const 
{ 
  return my_z_scores[id]; 
}

inline const vector<double>&
MemberCovariateCalculator::getDiffs() const
{
  return  my_diffs;
}


//======================================================================
// IMPLEMENTATION:  MccContinuousBoth
//======================================================================
//
inline
MccContinuousBoth::~MccContinuousBoth()
{}

//======================================================================
// IMPLEMENTATION:  MccContinuousDiff
//======================================================================
//
inline
MccContinuousDiff::~MccContinuousDiff()
{}

//======================================================================
// IMPLEMENTATION:  MccBinaryDiff
//======================================================================
//
inline
MccBinaryDiff::~MccBinaryDiff()
{}




