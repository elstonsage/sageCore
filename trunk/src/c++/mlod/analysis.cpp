#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "app/SAGEapp.h"
#include "mlod/analysis.h"

namespace SAGE
{
namespace MLOD
{

Analyzer::Analyzer(const RPED::RefMultiPedigree&   rped,
                   const RPED::genome_description& genome,
                   APP::Output_Streams&            output)
  : my_peds(rped),
    my_genome(genome),
    my_out(output)
{ }

AnalysisData Analyzer::run_analysis(const AnalysisParameters& params) const
{
  // Initialize our analysis outputs
  open_output_files(params);
  
  // Print our analysis headers
  print_analysis_headers(params);
  my_out.messages() << std::endl;
  
  DataShPtr data_impl = process_analysis(params);
  
  close_output_files();
  
  return AnalysisData(data_impl);
}

inline void
  output_sped_header(ostream& o, const string& sp)
{
  o << "  Generating LOD scores for " << sp << std::endl;
}

inline void output_sped_footer(ostream& o)
{
  o << "  =================================================================="
    << std::endl << std::endl;
}

inline boost::shared_ptr<std::ofstream>
open_output_file(const AnalysisParameters& params, const std::string& ext)
{
  string fname = params.get_base_output_filename() + "." + ext;
  
  boost::shared_ptr<std::ofstream> f(new std::ofstream(fname.c_str()));
  
  (*f) << APP::SAGEapp::getReleaseString() << std::endl;
  
  return f;
}

void Analyzer::open_output_files(const AnalysisParameters& params) const
{
  my_sum = open_output_file(params, "sum");
  my_det = open_output_file(params, "det");
}  

void Analyzer::close_output_files() const
{
  my_sum = boost::shared_ptr<std::ofstream>();
  my_det = boost::shared_ptr<std::ofstream>();
}

void Analyzer::print_analysis_headers(const AnalysisParameters& params) const
{
  my_out.messages() << std::endl << "Running Analysis: " << params.get_title() << std::endl
      << "==============================================================" << std::endl
      << std::endl;

  params.print_analysis_table(my_out.messages());
  params.print_analysis_table(*my_sum);
  params.print_analysis_table(*my_det);
  
  my_out.messages() << std::endl;
  *my_sum           << std::endl;
  *my_det           << std::endl;
}

Analyzer::DataShPtr
    Analyzer::generate_data_impl(const AnalysisParameters& params) const
{
  my_out.messages() << "Generating Valid Pedigree Data Sample..." << flush;

  DataShPtr data_impl(new AnalysisDataImpl(params, my_peds, my_out.errors()));

  my_out.messages() << "...Done." << std::endl << std::endl;

  // Print basic summary information
  data_impl->my_pedigree_sample.print_summary_statistics_table(my_out.messages());
  data_impl->my_pedigree_sample.print_summary_statistics_table(*my_sum);
  data_impl->my_pedigree_sample.print_summary_statistics_table(*my_det);

  // print detailed information if applicable
  if(params.get_ind_sample_table_option() != AnalysisParameters::IS_NONE)
  {
    bool detailed = params.get_ind_sample_table_option() == AnalysisParameters::IS_ALL;
    
    data_impl->my_pedigree_sample.print_member_table(my_out.info(), detailed);
    data_impl->my_pedigree_sample.print_member_table(*my_det,       detailed);
  }

  my_out.errors() << flush;
  my_out.messages() << flush;

  // IF data set empty, stop
  if(data_impl->my_pedigree_sample.get_pedigree_count() == 0)
  {
    my_out.errors() << priority(error);
    
    std::ostringstream s;
    
    s << "No valid peigree data for analysis "
      << params.get_title() << "." << std::endl;
      
    my_out.errors() << s.str() << flush;
    *my_sum         << s.str();
    *my_det         << s.str();
    
    // print error and return an empty analysis data
    return boost::shared_ptr<AnalysisDataImpl>();
  }
  
  return data_impl;
}

Analyzer::DataShPtr Analyzer::process_analysis(const AnalysisParameters& params2) const
{
  // Create our data store and analysis sample.
  DataShPtr data_impl = generate_data_impl(params2);
  
  if(!data_impl.get())
    return data_impl;
    
  // Create Analysis object and set it up.  If this fails, we just return 
  my_out.messages() << "Initializing Analysis..................." << flush;
  
  // Create the lod analysis object
  LodScoreAnalyzer analyzer(*data_impl, my_out);
  
  if(!analyzer.is_valid())
  {
    // print error and return an empty analysis data
    my_out.errors() << priority(error)
                    << "Initialization Error.  Aborting analysis." << std::endl;
    *my_sum         << "Initialization Error.  Aborting analysis." << std::endl;
    *my_det         << "Initialization Error.  Aborting analysis." << std::endl;
        
    return data_impl;
  }

  analyzer.set_result_target(data_impl->my_results);

  // Is this needed?
  if(!analyzer.is_valid())
  {
    // print error and return an empty analysis data
    my_out.errors() << priority(error)
                    << "Initialization Error.  Aborting analysis." << std::endl;
    *my_sum         << "Initialization Error.  Aborting analysis." << std::endl;
    *my_det         << "Initialization Error.  Aborting analysis." << std::endl;

    return data_impl;
  }
  
  my_out.messages() << "...Done." << std::endl;

  // Finally begin running our lod scores
  my_out.messages() << "Calculating LOD Scores.................." << std::endl << std::endl;
  
  LodTableFormatter formatter(my_peds.info(), data_impl->my_parameters);

  // Run the analysis for each subpedigree
  for(PedigreeAnalysisSample::SpedIterator spiter = data_impl->my_pedigree_sample.get_subpedigree_begin();
      spiter != data_impl->my_pedigree_sample.get_subpedigree_end(); ++spiter)
  {
    string sp_name;
    
    if((*spiter)->pedigree()->subpedigree_count() == 1)
      sp_name = "pedigree " + (*spiter)->pedigree()->name();
    else
    {
      std::ostringstream temp;
      
      temp << "constituent pedigree " << (*spiter)->index() << " of pedigree "
           << (*spiter)->pedigree()->name();
           
      sp_name = temp.str();
    }
    
    output_sped_header(my_out.messages(), sp_name);

    analyzer.analyze_subpedigree(spiter);
    
    if(data_impl->my_parameters.get_ped_output_detail_option() != AnalysisParameters::PD_NONE)
    {
      OUTPUT::Table t = formatter.formatTable
          (data_impl->my_results.get_sped_lod_table(spiter),
           sp_name,
           data_impl->my_parameters.get_ped_output_detail_option() != AnalysisParameters::PD_MARKERS);
           
      (*my_det) << t << flush;
    }
    
    output_sped_footer(my_out.messages());
  }
  
  my_out.messages() << "LOD Score Calculation complete." << std::endl << std::endl;

  // Output final lod tables
  
  OUTPUT::Table t = formatter.formatTable
      (data_impl->my_results.get_summary_lod_table(), "Analysis", true);
      
  *my_sum << t;
  *my_det << t;
  
  return data_impl;
}

}
}
