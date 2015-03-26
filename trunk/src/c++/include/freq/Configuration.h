#ifndef FREQ_CONFIGURATION
#define FREQ_CONFIGURATION

#include <string>
#include <iostream>

namespace SAGE {
namespace FREQ {

struct Configuration
{
  /// A vector of marker id's.
  typedef std::vector<size_t> MarkerVector;

  /// Constructor.
  Configuration() 
  { 
    ofilename           = "freq"; 
    markers             . clear(); 
    use_founder_weight  = false; 
    founder_weight      = 0.0; 
    skip_mle            = false; 
    estimate_inbreeding = false; 
    maxfun_debug        = false;
  }

  /// Copy constructor.
  Configuration(const Configuration & other) :
    ofilename           (other.ofilename),
    use_founder_weight  (other.use_founder_weight),
    founder_weight      (other.founder_weight),
    skip_mle            (other.skip_mle),
    estimate_inbreeding (other.estimate_inbreeding),
    maxfun_debug        (other.maxfun_debug),
    markers             (other.markers)
  { }

  /// operator=
  Configuration& operator=(const Configuration & other)
  {
    if(this != &other)
    {
      ofilename           = other.ofilename;
      use_founder_weight  = other.use_founder_weight;
      founder_weight      = other.founder_weight;
      skip_mle            = other.skip_mle;
      estimate_inbreeding = other.estimate_inbreeding;
      maxfun_debug        = other.maxfun_debug;
      markers             = other.markers;
    }

    return *this;
  }

  /// Adds the given id to the marker list if it isn't already in there.
  void addMarker(size_t id)
  {
    for(size_t i = 0; i < markers.size(); ++i)
      if(id == markers[i])
        return;
        
    markers.push_back(id);
  }

  /// Returns a const reference to the vector of markers.
  const MarkerVector & getMarkers() const { return markers; }

  /// Filename for output.
  std::string ofilename;

  /// Whether or not to use a founder weight.
  bool use_founder_weight;

  /// The founder weight to use.
  double founder_weight;

  /// Whether or not to skip maximum likelihood estimation.
  bool skip_mle;

  /// Whether or not to skip the estimation of the inbreeding coefficient.
  bool estimate_inbreeding;
  
  /// Whether or not to enable maxfun debugging output
  bool maxfun_debug;

  ///
  /// Dumps contents to screen (debugging use only)
  void dump() const
  {
    std::cout << "CONFIG DUMP"                                    << std::endl
              << "  ofilename          = " << ofilename           << std::endl
              << "  use founder weight = " << use_founder_weight  << std::endl
              << "  founder weight     = " << founder_weight      << std::endl
              << "  skip mle           = " << skip_mle            << std::endl
              << "  estimate inb.      = " << estimate_inbreeding << std::endl
              << "  markers            = ";    
          
              
    for(size_t i = 0; i < markers.size(); ++i)
      std::cout << markers[i] << " ";
      
    std::cout << std::endl << std::endl;
  }

private:

  /// List of markers to analyze.
  MarkerVector markers;
};

} // End namespace FREQ
} // End namespace SAGE

#endif
