#ifndef GROUP_AGGREGATE_PROCESS_H
#define GROUP_AGGREGATE_PROCESS_H

#include <string>
#include <map>
#include <vector>
#include <list>
#include <set>
#include <cctype>
#include <cstdlib>
#include <limits>
#include <cassert>
#include <setjmp.h>
#include <signal.h>
#include <unistd.h>
#include <iomanip>
#include <algorithm>
#include <utility>
#include <math.h>

#include "util/StringUtils.h"
#include "util/RegexUtils.h"
#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"
#include "error/internal_error.h"
#include "numerics/cephes.h"
#include "numerics/constants.h"
#include "numerics/sinfo.h"
#include "numerics/histogram.h"
#include "numerics/mt.h"
#include "numerics/functions.h"
#include "util/AutoTrace.h"
#include "func/FunctionParser.h"

namespace SAGE {
namespace FUNC {

/// Takes care of generating trait values for 'mean_adj', 'var_adj', and 'z_score'

/// Please note that if the function is not successful, it will report errors to the user
/// (via the errorstream passed to the function). Remember to check the return value
/// of createAdjustedTrait() !

/// This class has been designed to take advantage of function block's 'VERBOSE' parameter.
/// If 'VERBOSE' is included in the function block, there will be a fair amount of useful
/// debugging information provided in the .inf file.

/// The simplest operation of this class is for a non-class-based adjustment. If the
/// user only specified one argument to the adjustment function (see user manual), that
/// means that the trait adjustment is to be based on the summary data for the trait.

/// In this case, the class will simply populate its 'my_summary_stats' with valid
/// values from the rmp, then perform the requested adjustment on the basis of those stats.

/// If the user wants class-based adjustment, but with no minimum bin size ...
/// For each individual i, AdjustedTraitCreator will determine i's classification. The individual's
/// original trait value will be added to a SampleInfo JUST FOR THE CLASS IN QUESTION.
///
/// Later on, when it comes time to actually adjust the values, the stats used as the basis
/// for that adjustment will be the class-based stats, NOT the summary stats.
///
class AdjustedTraitCreator
{
public:
  struct TraitStats
  {
    double  binned_mean;
    double  binned_stdev;
    SampleInfo  sinfo;
  };

  typedef std::map<size_t, TraitStats> ClassBasedStats;
  typedef std::vector<size_t> TraitVector;

  /// \returns \c true if successful, \c false otherwise
  static bool createAdjustedTrait(RPED::RefMultiPedigree& rmp, cerrorstream errors, const FunctionParser::TraitData& data);
  
private:
  enum FunctionType { MEAN_ADJ, VAR_ADJ, Z_SCORE };

  /// \brief Indicates whether the trait to be created should be adjusted on the basis of a single trait, or on two traits
  /// AdjustmentType will be determined by parsing the function expression. If there is only one argument,
  /// the expression shall be interpreted as a SELF_ADJUSTMENT. If there are 3 or more, the expression shall
  /// be interpreted as ALTERNATE_ADJUSTMENT.
  enum AdjustmentType
  {
    /// In the case of non-class-based adjustment, a variable is adjusted by the overall descriptive statistic
    /// in question.
    ///
    /// That is, given \f$ theta \f$ (the new trait to create), \f$ X \f$ (the descriptive statistic about some
    /// trait 'x', and \f$ f(x_i | X) \f$ (some function operating on the individual's value for x), then:
    ///
    /// \f$ theta_i = f(x_i | X) \f$
    /// 
    /// For mean adjustment, for example:
    ///
    /// \f$ theta_i = x_i - X \f$
    ///
    NON_CLASS_BASED, 

    /// In the case of class-based adjustment, a variable is adjusted by the descriptive statistic
    /// regarding a subset of the data in question (on the basis of the individual's classification).
    ///
    /// That is, given \f$ theta \f$ (the new trait to create), \f$ j \f$ (the classification of the
    /// individual), \f$ X_j \f$ (the descriptive statistic about the class-based subset of some
    /// trait 'x', and \f$ f(x_i | X_j) \f$ (some function operating on the individual's value for x), then:
    ///
    /// \f$ theta_i = f(x_i | X_j) \f$
    /// 
    /// For mean adjustment, for example:
    ///
    /// \f$ theta_i = x_i - X_j \f$
    ///
    CLASS_BASED
  };

  AdjustedTraitCreator(RPED::RefMultiPedigree& rmp, cerrorstream errors);
    
  /// Creates an X-adjusted trait (where x is mean_adj, var_adj, or z_score). 
  bool internal_createAdjustedTrait(const FunctionParser::TraitData& data);
  
  /// Calculates the necessary statistics for the summary information, and for trait-based information.
  int calculateStats(const FunctionParser::TraitData& data);

  /// Having calculated the necessary statistics, creates the new trait and populates the individual
  /// trait values.
  void createIndividualNewTraitValues(const FunctionParser::TraitData& par_data);

  ///
  /// Populates the summary SampleInfo and the class-based TraitStats's with values.
  void populateSampleInfos(const FunctionParser::TraitData& data);

  /// Applies the binning algorithm (if required) to the class-based TraitStats's.
  bool performBinning(const FunctionParser::TraitData& data);

  /// Parses the expression into tokens (populating class trait nums, base trait num, bin size).
  /// \retval 0 Everything's ok!
  /// \retval 1 Malformed string
  /// \retval 2 At least one of the traits referred to doesn't exist
  /// \retval 3 The bin size argument is misspecified
  int parseExpressionString(const FunctionParser::TraitData& data);

  /// Figures out the means / variances / etc. for each class, as well as for the trait altogether.
  /// \retval 0 Everything calculated ok
  /// \retval 1 Not enough values to flesh out the bin size
  int calculateTraitStats();

  /// Adds the trait entry to the RefMPedInfo's and the RefPedInfo's.
  size_t addTrait(const FunctionParser::TraitData& par_data);


  RPED::RefMultiPedigree&  my_rmp;
  cerrorstream  my_errors;
  
  SampleInfo  my_summary_stats;
  ClassBasedStats  my_class_stats;
  
  size_t  my_reference_trait;  
  TraitVector  my_class_traits;
  
  size_t  my_bin_size;
  FunctionType  my_function_type;
  AdjustmentType  my_adj_type;
};

} // End namespace FUNC
} // End namespace SAGE

#endif
