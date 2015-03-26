#include "mlod/lod_score_analyzer.h"

namespace SAGE
{
namespace MLOD
{

LodScoreAnalyzer::LodScoreAnalyzer(
    const AnalysisDataImpl&               data,
    APP::Output_Streams&                  out)
  : my_data(data),
    my_parameters(data.my_parameters),
    my_peds(data.my_pedigree_sample),
    my_result_target(NULL),
    my_output(out),
    my_marker_data(out.errors()),
    my_trait_data (out.errors())
{ 
  // Initially set the my_valid flag to true.  It will be reset if anything goes wrong
  my_valid = true;
  
  // WORAROUND: Shouldn't need to const cast the genome to get the region, 
  // but since region() isn't provided const, we have to do it
  RPED::genome_description& g = const_cast<RPED::genome_description&>(data.my_genome);
  
  my_traits = g.region("TRAITS");
  
  // Perform data structure building.  Each step may make the whole thing invalid,
  // so check that for each operation, stoping if we hit a problem.
  if(is_valid()) request_marker_resources();
  if(is_valid()) build_likelihood_data_structures();
  if(is_valid()) build_temporary_likelihood_vectors();
}

void LodScoreAnalyzer::request_marker_resources()
{
  my_marker_data.request_resource(mpoint_likelihood_data::MP_COMBINED);

  if(using_intervals())
    my_marker_data.request_resource(mpoint_likelihood_data::MP_SEPARATE);

  my_trait_data.request_resource(mpoint_likelihood_data::SINGLE_POINT);
}
  
void LodScoreAnalyzer::build_likelihood_data_structures()
{
  //lint --e{713} Lots of conversion, but none of it matters
  
  // Allocate our data structures.  The marker and trait data are the biggest, so
  // do them first.
  if(!my_marker_data.build(num_loci(),   num_bits()) ||
     !my_trait_data.build (num_traits(), num_bits()) )
  {
    my_output.errors() << SAGE::priority(SAGE::critical)
                 << "Unable to allocate memory for analysis.  "
                    "Please reduce the maximum cutoff." 
                 << endl;

    my_valid = false;

    return;
  }

  if(!my_marker_data.built() || !my_trait_data.built())
  {
    my_output.errors() << SAGE::priority(SAGE::critical)
                 << "Unable to perform analysis." << endl;

    my_valid = false;
    
    return;
  }
}

void LodScoreAnalyzer::build_temporary_likelihood_vectors()
{
  // Now allocate our two temporary likelihood vectors
  temp1 = Likelihood_Vector(num_bits());
  temp2 = temp1;

  if(!temp1.is_valid() || !temp2.is_valid())
  {
    my_output.errors() << SAGE::priority(SAGE::critical)
           << "Unable to allocate memory for analysis.  "
              "Please reduce the maximum cutoff." << endl;

    my_valid = false;
    return;
  }
}

void LodScoreAnalyzer::initialize_meiosis_map(const FPED::Subpedigree& s)
{
  my_meiosis_map = meiosis_map(&s);
}

void LodScoreAnalyzer::initialize_marker_data(const FPED::Subpedigree& s)
{
  // Generate our pedigree specific marker data
  pedigree_region sped_markers = pedigree_region(s, my_parameters.get_region(), my_output.errors());

  // Initialize the likelihood vector data structure for markers.
  //lint -e{534} Ignored return
  my_marker_data.set_markers(sped_markers);
  
  my_marker_data.set_meiosis_map(my_meiosis_map);
}

void LodScoreAnalyzer::initialize_trait_data(const FPED::Subpedigree& s)
{
  // Generate our pedigree specific trait loci data
  pedigree_region sped_traits  = pedigree_region(s, my_traits,  my_output.errors());

  // Initialize the likelihood vector data structure for traits.
  //lint -e{534} Ignored return
  my_trait_data .set_markers(sped_traits);
  my_trait_data .set_meiosis_map(my_meiosis_map);
}

bool LodScoreAnalyzer::analyze_subpedigree(SpedIterator spiter)
{
  //lint --e{713} lots of sign issues, but none important
  
  if(!is_valid())
    SAGE_internal_error();

  if(!my_result_target)
    SAGE_internal_error();

  initialize_meiosis_map   (**spiter);
  initialize_marker_data   (**spiter);
  initialize_trait_data    (**spiter);

  // Verify that the initalization worked.
  if(!my_marker_data.valid() || !my_trait_data.valid())
  {
    // Print error?
    return false;
  }
  
  SpedLodTable lod_table(num_traits(), num_points(), 1 << my_meiosis_map.nonfounder_meiosis_count());

  // Generate our stored variables

  size_t sptable_index = 0;

  //cout << sptable_index << "\t";

  // Compute the lod scores at the first marker
  compute_lod_scores(lod_table, my_marker_data.multi_point_vector(0), sptable_index);

  ++sptable_index;

  for(size_t marker_idx = 0; marker_idx < num_loci() - 1; ++marker_idx)
  {
    if(using_intervals())
    {
      // Determine how many points exist between marker marker_idx and marker_idx+1
      size_t num_pts_in_intval = 
          my_parameters.get_region().locus(marker_idx).interval_point_count(1);
      
      // For each point, compute a multipoint likelihood vector and then
      // calculate the lod scores from that.
      for(size_t pt = 1; pt < num_pts_in_intval; ++pt)
      {
        compute_interval_multipoint_lvec(marker_idx, pt);
        
        //cout << sptable_index << "\t";

        compute_lod_scores(lod_table, temp1, sptable_index);
        ++sptable_index;
      }
    }

    //cout << sptable_index << "\t";

    // Compute the lod scores for the next marker.
    compute_lod_scores(lod_table,
                       my_marker_data.multi_point_vector(marker_idx+1),
                       sptable_index);

    ++sptable_index;
  }
  
  // Add the SpedLodTable into the results
  //lint -e{613} We check for null above
  my_result_target->set_sped_lod_table(spiter, lod_table);

  return true;
}

void LodScoreAnalyzer::compute_lod_scores
    (SpedLodTable&            slt,
     const Likelihood_Vector& mkr_vect,
     size_t                   sptable_index)
{
  // Compute the no linkage amount, the average likelihood for each inheritance
  // pattern.
  double no_linkage = mkr_vect.total() / mkr_vect.size();

#if 0
  cout << mkr_vect.total() << "\t"
       << mkr_vect.log_scale() << "\t"
       << mkr_vect.log_scale().get_log() << "\t"
       << mkr_vect.size() << "\t";
#endif

  if(no_linkage == 0) return;

  // If the marker has non-zero likelihoods, we calculate the multipoint likelihood
  // for each value.
  for(long i = 0; i < my_trait_data.lvector_count(); ++i)
  {
    // Use temp2 here as temp1 may *be* the mkr_vect, and we don't want to
    // corrupt it.
    temp2 = mkr_vect;

    temp2 *= my_trait_data.single_point_vector(i);

#if 0
    cout << my_trait_data.single_point_vector(i).total() << "\t"
         << my_trait_data.single_point_vector(i).log_scale() << "\t"
         << my_trait_data.single_point_vector(i).log_total() << "\t"
         << temp2.total() << "\t"
         << temp2.log_scale() << "\t"
         << temp2.log_scale().get_log() << "\t"
         << temp2.size() << "\t";
#endif

    double scaling_factor = (temp2.log_scale() / mkr_vect.log_scale()).get_double();
    double lod_score = log10((temp2.total() / no_linkage) * scaling_factor);

#if 0
    cout << no_linkage << "\t"
         << temp2.log_scale() / mkr_vect.log_scale() << "\t"
         << temp2.log_scale().get_log() - mkr_vect.log_scale().get_log() << "\t"
         << scaling_factor << "\t"
         << lod_score << "\t"
         << temp2.information() << endl;
#endif

    SpedLodTable::LodScoreInfoType lsi(lod_score, temp2.information());

    //lint -e{732} loss of sign ok
    slt.set_lod_score_info(i, sptable_index, lsi);
  }
}

void LodScoreAnalyzer::compute_interval_multipoint_lvec(size_t marker_idx, size_t pt)
{
  //lint --e{713}, --e{534}
  temp1 = my_marker_data. left_sided_vector(marker_idx    );
  temp2 = my_marker_data.right_sided_vector(marker_idx + 1);

  // Get the region for convenience
  const AnalysisParameters::RegionType& r = my_parameters.get_region();
  
  // Multiply temp1 times the rec frac to pt from marker marker_idx, and
  // multiply temp2 times the rec frac to pt from marker marker_idx+1, and
  temp1(&my_meiosis_map, r.locus(marker_idx).point_theta(pt));
  temp2(&my_meiosis_map, r.locus(marker_idx+1).point_theta(pt - r.locus(marker_idx).interval_point_count(1)));

  // temp1 is now the left hand vector, and temp2 the right.  Multiply them together toget
  // the multipoint vector at pt.
  temp1 *= temp2;
}

}
}

