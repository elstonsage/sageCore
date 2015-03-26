#ifndef FREQ_RESULTS
#define FREQ_RESULTS

#include "maxfunapi/maxfunapi.h"
#include "output/Output.h"
#include "freq/Configuration.h"
#include "freq/Sample.h"

namespace SAGE {
namespace FREQ {

struct AlleleFrequencies
{
  AlleleFrequencies() { founder_freqs.clear(); all_freqs.clear(); }
  
  AlleleFrequencies(const AlleleFrequencies & other) : founder_freqs(other.founder_freqs), all_freqs(other.all_freqs) { }

  AlleleFrequencies& operator= (const AlleleFrequencies & other) { if(this != &other) { founder_freqs = other.founder_freqs; all_freqs = other.all_freqs; } return *this; }  
  
  void dump() const
  {
    std::cout << "ALLELE FREQUENCIES DUMP" << std::endl;
    
    for(size_t i = 0; i < founder_freqs.size(); ++i)
    {
      std::cout << i << " founders only: " << std::setprecision(10) << founder_freqs[i] << " everyone: " << all_freqs[i] << std::endl;
    }
      
    std::cout << std::endl;
  }

  /// Frequencies for each allele for only the founders
  std::vector<double> founder_freqs;

  /// Frequencies for each allele for the entire dataset
  std::vector<double> all_freqs;
};

struct MarkerAnalysis
{
  MarkerAnalysis() { }
  
  MarkerAnalysis(const MarkerAnalysis & other) :
    marker_id                 (other.marker_id),
    allele_freqs              (other.allele_freqs),
    maxfun_results            (other.maxfun_results),
    maxfun_results_inbreeding (other.maxfun_results_inbreeding)
    { }

  MarkerAnalysis& operator=(const MarkerAnalysis & other)
  {
    if(this != &other)
    {
      marker_id                 = other.marker_id;
      allele_freqs              = other.allele_freqs;
      maxfun_results            = other.maxfun_results;
      maxfun_results_inbreeding = other.maxfun_results_inbreeding;
    }

    return *this;
  }

  /// Marker id analyzed.
  size_t marker_id;
  
  /// The initial allele frequencies.
  AlleleFrequencies allele_freqs;
    
  /// The maximization results (without inbreeding coefficient).
  MAXFUN::Results maxfun_results;
  
  /// The maximization results (with inbreeding coefficient).
  MAXFUN::Results maxfun_results_inbreeding;
};

typedef std::vector<MarkerAnalysis> MarkerAnalyses;

struct Results
{
  Results() { }
  
  Results(const Results & other) : config(other.config), analyses(other.analyses) { }
  
  Results& operator=(const Results & other) { if(this != &other) { config = other.config; analyses = other.analyses; } return *this; }
  
  /// Configuration used for analysis
  Configuration config;

  /// Analyses of each requested marker
  MarkerAnalyses analyses;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
