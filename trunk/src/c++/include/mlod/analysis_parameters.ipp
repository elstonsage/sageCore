#ifndef MLOD_ANALYSIS_PARAMETERS
#include "mlod/analysis_parameters.h"
#endif

namespace SAGE
{
namespace MLOD
{

//============================================================================
// IMPLEMENTATION:  analysis_parameters
//============================================================================
//
inline
AnalysisParameters::AnalysisParameters
    (size_t analysis_id)
  : my_scan_type       (ST_MARKER),
    my_ped_option      (PD_NONE),
    my_ind_option      (IS_REMOVED),
    my_distance        (2.0),
    my_max_ped_size    (MLOD_DEFAULT_MAX_PED_SIZE),
    my_region          (),
    my_trait_loci      ()
{
  string  id = long2str((long) analysis_id);

  my_analysis_title = "MLOD Analysis " + id;
  my_base_filename  = "mlod_analysis" + id;
}
inline
AnalysisParameters::AnalysisParameters
  (const AnalysisParameters& a)
  : my_analysis_title  (a.my_analysis_title),
    my_base_filename   (a.my_base_filename),
    my_scan_type       (a.my_scan_type),       
    my_ped_option      (a.my_ped_option),      
    my_ind_option      (a.my_ind_option),
    my_distance        (a.my_distance),        
    my_max_ped_size    (a.my_max_ped_size),
    my_region          (a.my_region),         
    my_trait_loci      (a.my_trait_loci)      
{ }
                        
inline AnalysisParameters&
AnalysisParameters::operator=(const AnalysisParameters& a)
{
  if(&a != this)
  {
    my_analysis_title   = a.my_analysis_title;
    my_base_filename    = a.my_base_filename;
    my_scan_type        = a.my_scan_type;       
    my_ped_option       = a.my_ped_option;      
    my_ind_option       = a.my_ind_option;      
    my_distance         = a.my_distance;        
    my_max_ped_size     = a.my_max_ped_size;
    my_region           = a.my_region;         
    my_trait_loci       = a.my_trait_loci;      
  }

  return *this;
}

inline
AnalysisParameters::~AnalysisParameters()
{ }

inline bool
AnalysisParameters::is_valid() const
{
  return !my_trait_loci.empty() && 
          my_region.valid();
}

inline const string& 
AnalysisParameters::get_title() const
{
  return my_analysis_title;
}

inline void 
AnalysisParameters::set_title(const string& s)
{
  my_analysis_title = s;
}

inline const string&
AnalysisParameters::get_base_output_filename() const
{
  return my_base_filename;
}

inline void 
AnalysisParameters::set_base_output_filename(const string& s)
{
  my_base_filename = s;
}

inline AnalysisParameters::PedOutputDetailEnum 
AnalysisParameters::get_ped_output_detail_option() const
{
  return my_ped_option;
}

inline
string AnalysisParameters::get_ped_output_detail_string() const
{
  switch(my_ped_option)
  {
    case PD_NONE      : return "none";
    case PD_MARKERS   : return "markers";
    case PD_INTERVALS : return "intervals";
    case PD_ALL       : return "all";
  }
  
  SAGE_internal_error();
  
  return "";
}

inline void 
AnalysisParameters::set_ped_output_detail_option(AnalysisParameters::PedOutputDetailEnum p)
{
  my_ped_option = p;
}

inline AnalysisParameters::IndSampleTableEnum 
AnalysisParameters::get_ind_sample_table_option() const
{
  return my_ind_option;
}

inline
string AnalysisParameters::get_ind_sample_table_string() const
{
  switch(my_ind_option)
  {
    case IS_NONE    : return "none";
    case IS_REMOVED : return "removed";
    case IS_ALL     : return "all";
  }
  
  SAGE_internal_error();
  
  return "";
}

inline void 
AnalysisParameters::set_ind_sample_table_option(AnalysisParameters::IndSampleTableEnum i)
{
  my_ind_option = i;
}

inline AnalysisParameters::ScanTypeEnum 
AnalysisParameters::get_scan_type() const
{
  return my_scan_type;
}

inline string
AnalysisParameters::get_scan_type_string() const
{
  switch(my_scan_type)
  {
    case ST_MARKER   : return "markers";
    case ST_INTERVAL : return "intervals";
    case ST_BOTH     : return "markers and intervals";
  }
  
  SAGE_internal_error();

  return "";
}

inline void 
AnalysisParameters::set_scan_type(ScanTypeEnum s)
{
  my_scan_type = s;
}

inline double 
AnalysisParameters::get_scan_distance() const
{
  return my_distance;
}

inline void
AnalysisParameters::set_scan_distance(double d)
{
  my_distance = d;
}

inline size_t 
AnalysisParameters::get_max_ped_size() const
{
  return my_max_ped_size;
}

inline void
AnalysisParameters::set_max_ped_size(size_t sz)
{
  my_max_ped_size = sz;
}

inline const AnalysisParameters::RegionType& 
AnalysisParameters::get_region() const
{
  return my_region;
}

inline void
AnalysisParameters::set_region(const RegionType& r)
{
  my_region = r;
}

inline AnalysisParameters::TraitModelList& 
AnalysisParameters::get_trait_list()
{
  //lint -e{1536}  Exposing this data member is both ok, and desired.
  return my_trait_loci;
}

inline const AnalysisParameters::TraitModelList&
AnalysisParameters::get_trait_list() const
{
  return my_trait_loci;
}

inline
size_t get_trait_model_pair_first_element(const AnalysisParameters::TraitModelPair& p)
{
  return p.first;
}

inline
bool is_same_index(const AnalysisParameters::TraitModelPair& p, size_t idx)
{
  return p.first == idx;
}

}
}
