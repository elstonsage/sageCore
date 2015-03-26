#ifndef MAXIMIZER_H
#define MAXIMIZER_H

#include "maxfun/maxfun.h"
#include "maxfunapi/Function.h"
#include "maxfunapi/SequenceCfg.h"
#include "maxfunapi/DebugCfg.h"
#include "maxfunapi/ParameterMgr.h"
#include "maxfunapi/APIMaxFunction.h"
#include "maxfunapi/Results.h"

namespace SAGE   {
namespace MAXFUN {

/**
 * \internal
 * \class Maximizer
 * \brief The maximizing class.
 */
class Maximizer
{
  public:
    //===================================================================
    // Public utility functions:
    //===================================================================

    ///
    /// This function performs a maximization of the function \c func according to configuration options \c config.
    /// \param mgr The ParameterMgr object that manages the parameters.
    /// \param func The MaxFunction object that returns the function value.
    /// \param config The SequenceCfg object that describes how to carry out this maximization.
    /// \param debug A DebugCfg object that describes the user's debugging preferences.
    /// \returns A Results object containing information about the maximization requested.
    static Results Maximize(
      const SequenceCfg  & config,  
            ParameterMgr & mgr, 
            MaxFunction  & func, 
      const DebugCfg     & debug);
      
    static double  calculateScoreTestStatistic(Function& func, const DebugCfg& debug);

    ///
    /// Performs a default maximization of the given Function.
    /// \param config The SequenceCfg object that describes how to carry out this maximization.
    /// \param func A MAXFUN::Function to maximize.
    /// \param debug A DebugCfg object that describes the user's debugging preferences.
    /// \returns A Results object containing information about the maximization requested.
    static Results maximize(
            Function    & func,
      const SequenceCfg & config = SequenceCfg (SequenceCfg::DEFAULT_MAXIMIZATION), 
      const DebugCfg    & debug  = DebugCfg    (DebugCfg::NO_DEBUG_INFO));

    ///
    /// This functions performs a single function evaluation of the given function and returns the results.
    /// Please note that it is entirely possible for this function to return INF or QNAN!
    /// \param mgr The ParameterMgr object that manages the parameters
    /// \param func The MaxFunction obejct that returns the function value
    /// \param use_initial_ests If \c true, uses the initial estimates from the ParameterMgr. If \c false, uses the 
    /// \b current estimates from the ParameterMgr.
    /// \returns A pair where the first elements is the function value, and the second element is the integer exit code (see Maxfun manual)
    static pair<double, int> evaluateOnce(
      ParameterMgr & mgr,
      MaxFunction  & func,
      bool           use_initial_ests = true);

    ///
    /// This functions performs a single function evaluation of the given function and returns the results.
    /// Please note that it is entirely possible for this function to return INF or QNAN!
    /// \param func The MAXFUN::Function object that returns the function value
    /// \param use_initial_ests If \c true, uses the initial estimates from the ParameterMgr. If \c false, uses the 
    /// \b current estimates from the ParameterMgr.
    /// \returns A pair where the first elements is the function value, and the second element is the integer exit code (see Maxfun manual)
    static pair<double, int> evaluateOnce(
      Function  & func,
      bool        use_initial_ests = true);

  private:
    //===================================================================
    /// Private utility functions:
    //===================================================================

    ///
    /// Adds the parameters from the ParameterMgr to the maxfun instance.
    /// \param mgr The ParameterMgr
    /// \param maxfun The maxfun instance
    /// \param use_initial_ests If \c true, uses the initial estimates from the 
    /// ParameterMgr. If \c false, uses the \b current estimates from the ParameterMgr.
    static void AddParameters(
      ParameterMgr & mgr,
      Maxfun       & maxfun,
      bool           use_initial_ests = true);

    ///
    /// Adds the give parameter to the maxfun instance.
    /// \param param The Parameter to add
    /// \param maxfun The maxfun instance
    /// \param use_initial_est If \c true, uses the initial estimate from the
    /// Parameter. If \c false, uses the \b current estimate from the Parameter.
    static void AddParameter(
      Parameter     & param,   
      Maxfun        & maxfun,
      bool            use_initial_est = true);

    static void setParamIndices(ParameterMgr& mgr, const Maxfun_Data& data, Results& results);

    static Results MaximizeFunction(
            ParameterMgr   & mgr,
      const SequenceCfg    & config,
            Maxfun         & maxfun,
            APIMaxFunction & api_maxfunction,
      const DebugCfg       & debug);
};

/// @name Maximization functions
//@{

///
/// This function performs a maximization of the function \c func according to configuration options \c config.
/// \param mgr The ParameterMgr object that manages the parameters.
/// \param func The MaxFunction object that returns the function value.
/// \param config The SequenceCfg object that describes how to carry out this maximization.
/// \param debug A DebugCfg object that describes the user's debugging preferences.
/// \returns A Results object containing information about the maximization requested.
Results maximize(
        ParameterMgr   & mgr, 
        MaxFunction    & func,
  const SequenceCfg    & config, 
  const DebugCfg       & debug);

///
/// This functions performs a single function evaluation of the given function and returns the results.
/// Please note that it is entirely possible for this function to return INF or QNAN!
/// \param mgr The ParameterMgr object that manages the parameters
/// \param func The MaxFunction obejct that returns the function value
/// \param use_initial_ests If \c true, uses the initial estimates from the ParameterMgr. If \c false, uses the
/// \b current estimates from the ParameterMgr.
/// \returns A pair where the first elements is the function value, and the second element is the integer exit code (see Maxfun manual)
pair<double, int> evaluateOnce(ParameterMgr & mgr, MaxFunction  & func, bool use_initial_ests = true);
                                   

///
/// This function performs a default maximization of the function \c func.
///
/// Please note that invoking \code MaximizeDefault(my_mgr, my_func) \endcode is equivalent to invoking
/// \code Maximize(MAXFUN::SequenceCfg(MAXFUN::max_default), my_mgr, my_func, MAXFUN::DebugCfg(MAXFUN::NO_DEBUG_INFO)) \endcode
Results maximizeDefault(ParameterMgr & mgr, MaxFunction & func);

//@}

} // End namespace MAXFUN
} // End namespace SAGE

#endif
