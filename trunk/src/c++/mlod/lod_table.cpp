#include "numerics/isnan.h"
#include "mlod/lod_table.h"

namespace SAGE {
namespace MLOD {

bool SummaryLodTable::add_table (const SpedLodTable& l)
{
  if(my_num_points != l.my_num_points || my_num_traits != l.my_num_traits)
    return false;

  for(size_t i = 0; i < my_lod_scores.size(); ++i)
  {
    const LodScoreInfoType& lsinfo = l.my_lod_scores[i];
    
    if(SAGE::isnan(lsinfo.lod_score) || SAGE::isnan(lsinfo.info_content))
      continue;

    // If we haven't had this particular element added before, we initialize the element
    if(my_lod_scores[i].table_count == 0)
    {
      my_lod_scores[i].info        = lsinfo;
      my_lod_scores[i].dimension   = l.get_dimension();
      my_lod_scores[i].table_count = 1;
    }
    else
    {
      // Otherwise, we must add in the current values and calculate the new information
      // content
      
      size_t new_dim = my_lod_scores[i].dimension + l.get_dimension();
  
      my_lod_scores[i].info.lod_score += lsinfo.lod_score;
  
      // new information content is the weighted sum of the current information content
      // and the information content of the new table.
      my_lod_scores[i].info.info_content = 
          (my_lod_scores[i].info.info_content * my_lod_scores[i].dimension + 
           lsinfo.info_content                * l.get_dimension()) 
           / new_dim;

      // Set the new dimension for future weighted sums.
      my_lod_scores[i].dimension = new_dim;

      // Add one table.  Stir.
      ++my_lod_scores[i].table_count;
    }
  }

  return true;
}

}
}

