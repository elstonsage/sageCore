#include "mlod/analysis_data.h"
#include "boost/bind.hpp"

namespace SAGE {
namespace MLOD {

void AnalysisDataImpl::setup_genome()
{
  // 1. Get the region from my_parameters
  RPED::genome_description::region_type orig_region = my_parameters.get_region();
  
  // 2. Initialize genome variables
  my_genome.set_scan_distance(my_parameters.get_scan_distance());
  my_genome.set_mapping_function(orig_region.map());
  
  // 3. Build Analysis Specific Regions
  build_analysis_specific_region();
  build_analysis_trait_region();
  
  // 4. Finalize the genome

  //lint -e{534} Ignored return
  my_genome.build();

  //lint -e{534} Ignored return
  my_genome.freeze();
  
  // 5. Connect my_parameters to the region
  my_parameters.set_region(my_genome.region("REGION"));
}

void AnalysisDataImpl::build_analysis_specific_region()
{
  // 1. Add the region to my_genome
  
  //lint -e{534} Ignored return
  my_genome.add_region("REGION");
  
  // 2. Get the region from my_parameters
  RPED::genome_description::region_type orig_region = my_parameters.get_region();

  for(size_t i = 0; i < orig_region.locus_count(); ++i)
  {
    my_genome.add_locus(orig_region.locus(i).marker_index(),
                        orig_region.locus(i).location());
  }
}

void AnalysisDataImpl::build_analysis_trait_region()
{
  //lint -e{534} Ignored return
  my_genome.add_region("TRAITS");

  // Add each trait in the AnalysisPameters to the trait genome
  //lint -e{534} Ignored return
  for_each(my_parameters.get_trait_list().begin(),
           my_parameters.get_trait_list().end(),
           boost::bind(&RPED::genome_description::add_locus, boost::ref(my_genome), 
               boost::bind(get_trait_model_pair_first_element, _1), 0.0));
}

}} // End namespaces
