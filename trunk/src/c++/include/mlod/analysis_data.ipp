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
// Copyright (c) 2005 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#ifndef ANALYSIS_DATA_H
#include "mlod/analysis_data.h"
#endif

namespace SAGE {
namespace MLOD {

inline
AnalysisData::AnalysisData()
{ }

inline
AnalysisData::AnalysisData(const AnalysisData& ad)
  : my_data_storage(ad.my_data_storage)
{ }

inline
AnalysisData::~AnalysisData()
{ }

inline
AnalysisData& AnalysisData::operator=(const AnalysisData& ad)
{
  if(this == &ad) return *this;
  
  my_data_storage = ad.my_data_storage;
  
  return *this;
}

inline
const AnalysisParameters& AnalysisData::get_parameters() const
{
  return my_data_storage->my_parameters;
}

inline
const PedigreeAnalysisSample& AnalysisData::get_pedigree_sample()  const
{
  return my_data_storage->my_pedigree_sample;
}

inline
const AnalysisResults&      AnalysisData::get_results()  const
{
  return my_data_storage->my_results;
}

inline
AnalysisData::AnalysisData(const AnalysisDataImplShPtr& data_storage)
  : my_data_storage(data_storage)
{ }
  
inline
AnalysisDataImpl::AnalysisDataImpl
   (const AnalysisParameters&  params,
    const RPED::MultiPedigree& source_rped,
    cerrorstream&              errors)
  : my_parameters      (params),
    my_pedigree_sample (source_rped, params),
    my_results         (my_pedigree_sample),
    my_genome          (source_rped.info(), errors)
{ 
  setup_genome();
}

}
}

