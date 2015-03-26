#ifndef LOD_SCORE_ANALYZER_H
#include "mlod/lod_score_analyzer.h"
#endif

namespace SAGE {
namespace MLOD {

inline
LodScoreAnalyzer::~LodScoreAnalyzer()
{
  my_result_target = 0;
}

inline
bool LodScoreAnalyzer::is_valid() const
{
  return my_valid;
}
  
inline
void LodScoreAnalyzer::set_result_target(AnalysisResults& results)
{
  my_result_target = &results;
}
  
inline
size_t LodScoreAnalyzer::num_loci()   const
{
  return my_parameters.get_region().locus_count();
}

inline
size_t LodScoreAnalyzer::num_points() const
{
  if(using_intervals())
  {
    return my_parameters.get_region().point_count();
  }
  else
  {
    return num_loci();
  }
}
inline
size_t LodScoreAnalyzer::num_traits() const
{
  return my_parameters.get_trait_list().size();
}

inline
size_t LodScoreAnalyzer::num_bits()   const
{
  return my_peds.get_largest_bit_count();
}

inline
bool LodScoreAnalyzer::using_intervals() const
{
  return my_parameters.get_scan_type() != AnalysisParameters::ST_MARKER;
}

} // end namespace MLOD
} // end namespace SAGE

