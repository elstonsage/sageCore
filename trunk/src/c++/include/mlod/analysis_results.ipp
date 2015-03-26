//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Geoff Wedig
//                                                                          
// History:   5/24/02 created                     - gcw
//            5/28/02 modified                    - djb
//                                                                          
// Notes      Defines classes to hold analysis parameters and analysis data.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef ANALYSIS_RESULTS_H
#include "mlod/analysis_results.h"
#endif

namespace SAGE {
namespace MLOD {

inline
AnalysisResults::AnalysisResults(const PedigreeAnalysisSample& ped_data)
  : my_pedigree_sample(ped_data),
    my_sped_lod_scores(ped_data.get_subpedigree_count())
{ }

inline
void AnalysisResults::set_sped_lod_table(SpedIterator spiter, const SpedLodTable& table)
{
  // Determine the index of the subpedigree
  //lint -e{732} loss of sign ok
  size_t spindex = spiter - my_pedigree_sample.get_subpedigree_begin();
  
  // Set the subpedigree's table
  my_sped_lod_scores[spindex] = table;
  
  // Check to see if the summary_lod_scores table has been initialized.
  if(my_summary_lod_scores.get_point_count() == 0 &&
     my_summary_lod_scores.get_trait_count() == 0)
  {
    my_summary_lod_scores = SummaryLodTable(table.get_trait_count(), table.get_point_count());
  }
  
  // Add the table to the summary table.
  //lint -e{534} Ignored return type
  my_summary_lod_scores.add_table(table);
}

inline
const SpedLodTable& AnalysisResults::get_sped_lod_table(SpedIterator spiter) const
{
  //lint -e{732} loss of sign ok
  size_t spindex = spiter - my_pedigree_sample.get_subpedigree_begin();
  
  return my_sped_lod_scores[spindex];
}

inline
const SummaryLodTable& AnalysisResults::get_summary_lod_table() const
{
  return my_summary_lod_scores;
}


}
}

