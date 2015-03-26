#include <sstream>
#include "mlod/analysis_parameters.h"

namespace SAGE
{
namespace MLOD
{

/// \internal
/// Returns a comma delimited, parentheses enclosed list
/// of traits.
///
/// \param o     The ostream where output goes
/// \param tlist The trait list.
inline
std::string create_trait_list_string
    (const AnalysisParameters::TraitModelList& tlist)
{
  std::stringstream s;
  
  // Output the singular or plural traits as needed
  if(tlist.size() == 1)
  {
    s << "trait (" << tlist.begin()->second.name() << ")";
  } 
  else
  {
    AnalysisParameters::TraitModelList::const_iterator iter = tlist.begin(); 
  
    s << "traits: (" << iter->second.name();
    
    ++iter;
    
    for( ; iter != tlist.end(); ++iter)
    {
      s << ", " << iter->second.name();
    }
    
    s << ")";
  }
  
  return s.str();
}


void AnalysisParameters::print_analysis_table(ostream& out) const
{
  std::stringstream s;
  
  // Create a list (bulletted) for the analysis description to go into
  OUTPUT::List desc("Analysis: " + my_analysis_title);
  
  // 1.  basic analysis on regions with traits
  
  s << "LOD Score analysis on region (" << get_region().name() << ") with "
    << create_trait_list_string(get_trait_list());
    
  desc << OUTPUT::List::makeBullet(s.str());
  s.str("");
  
  // 2.  Scanning
  switch(get_scan_type())
  {
    case ST_MARKER   :
      s << "Scan only markers.";
      break;
    case ST_INTERVAL :
      s << "Scan intervals at a distance of "
        << get_scan_distance();
      break;
    case ST_BOTH :
      s << "Scan markers and intervals at a distance of "
        << get_scan_distance();
      break;
  }

  desc << OUTPUT::List::makeBullet(s.str());
  s.str("");
  
  // 3.  Max Size
  s << "Maximum 2n-f pedigree size cutoff is "
    << get_max_ped_size();

  desc << OUTPUT::List::makeBullet(s.str());
  s.str("");
  
  // 4.  Pedigree Lods
  s << "Individual Pedigree LOD scores ";
  
  switch(get_ped_output_detail_option())
  {
    case PD_NONE      : s << "will not be printed."                                << endl; break;
    case PD_MARKERS   : s << "will be printed for the marker locations."           << endl; break;
    case PD_INTERVALS : s << "will be printed for the interval locations."         << endl; break;
    case PD_ALL       : s << "will be printed for markers and interval locations." << endl; break;
  }
  desc << OUTPUT::List::makeBullet(s.str());
  s.str("");
  
  // 5. Filenames
  s << "Output will be sent to files of the form: " << get_base_output_filename()
    << ".(det,sum)" << endl;

  desc << OUTPUT::List::makeBullet(s.str());
  
  out << desc << endl;
}

}
}

