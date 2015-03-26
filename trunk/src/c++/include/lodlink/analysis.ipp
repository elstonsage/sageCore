//============================================================================
// File:      analysis.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/20/2 created        -djb
//                                                                          
// Notes:     Inline implementation of analysis.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  analysis
//============================================================================
//
inline
analysis::analysis(cerrorstream& errors, const RPED::RefMultiPedigree& mped, const instructions& instr)
      : my_errors(errors), my_original_mped(mped), my_mped(mped), my_instructions(instr),
        aborted(false)
{
  build_filtered_mped();
}

// - Deletion of tasks not necessary because boost smart_ptr's are 
//   used.
//
inline
analysis::~analysis()
{}





