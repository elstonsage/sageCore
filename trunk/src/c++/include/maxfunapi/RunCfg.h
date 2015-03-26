#ifndef RUNCFG_H
#define RUNCFG_H

#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/DebugCfg.h"

namespace SAGE   {
namespace MAXFUN {

/** \class RunCfg
 *  \brief Describes a \b single run of Maxfun.
 *
 * \par Purpose
 *
 * RunCfg allows you to configure a single run of Maxfun.
 *
 * \par Detailed explanation
 *
 * It is important to remember that a complete maximization (invoking Maximize()) represents
 * a \b sequence of individual maximizations. RunCfg represents a configuration for a single
 * individual maximization.
 *
 * \par Getting started
 *
 * Generally speaking, you will add RunCfg's through your SequenceCfg object. SequenceCfg
 * has two functions (SequenceCfg::AddRunCfg() and SequenceCfg::getLatestRunCfg()) that
 * allow you to add and modify single RunCfg's.
 *
 * With SequenceCfg::AddRunCfg(), you can add another RunCfg to the sequence. This function
 * requires two arguments (maximization method and number of iterations). Please consult
 * SequenceCfg::AddRunCfg for more information.
 *
 * With SequenceCfg::getLatestRunCfg(), you can get a reference to the most recently added
 * RunCfg. With this reference, you can then directly change maximization configuration features.
 * Please consult the documentation on each feature for more specific information.
 */
class RunCfg
{
public:

  /// Identifies which type of maximization method to use.
  /// Please note that this is different from MAXFUN::maxtype, which indicates an
  /// entire preset maximization template, not simply a maximization method.
  enum MaximizationMethodEnum 
  {
    COMPLETE_DIRECT       = 1, /*!< \hideinitializer Complete direct search  */
    DIRECT_WITHOUT        = 2, /*!< \hideinitializer Direct search without \f$2^{N_i}\f$-trial search  */
    NEWTON_WITHOUT_SECOND = 3, /*!< \hideinitializer Newton-Raphson method, without repeated calculation of second derivatives  */
    NEWTON_WITH_SECOND    = 4, /*!< \hideinitializer Newton-Raphson method with second derivatives recalculated for each iteration  */
    VAR_METRIC_IDENTITY   = 5, /*!< \hideinitializer Variable metric method, with initial B = identity  */
    VAR_METRIC_ESTIMATE   = 6  /*!< \hideinitializer Variable metric method, with initial B = -H calculated at initial parameter estimates.  */
  };

  /// Variance-covariance matrix calculation option
  /// Note that the variance-covariance matrix is available \b only in case of convergence (including that
  /// presumed in the user-specified zero-iteration case).
  enum VarCovOptionEnum 
  {
    // Unspecifed is used internally only. If the user does not explicitly set the var_cov option, then
    // a default will be used. The UNSPECIFIED value is there so the program can check to see if the option
    // was set at all.
    VAR_COV_UNSPECIFIED = -1,
    NO_MATRIX = 0, /*!< \hideinitializer Not computed beyond what is necessary for parameter estimation. */
    INITIAL   = 1, /*!< \hideinitializer Computed for final estimates unless one is available that is no more than one iteration old. */
    FINAL     = 2  /*!< \hideinitializer Computed for final estimates (using iterative method). */
  };

  /// Second derivative 
  /// approximation calculation option
  enum SecDerivOptionEnum
  {
    // Unspecifed is used internally only. If the user does not explicitly set the second_deriv option, then
    // a default will be used. The UNSPECIFIED value is there so the program can check to see if the option
    // was set at all.
    SEC_DERIV_UNSPECIFIED = -1,
    SINGLE    = 0, /*!< \hideinitializer Use a single approximation */
    ITERATIVE = 1  /*!< \hideinitializer Estimate through an iterative process */
  };

  /// Control feature 
  /// for a RunCfg.
  enum MaxfunControlOptionEnum
  {
    UNCONDITIONAL           = 0, /*!< \hideinitializer Do unconditionally. */
    PREVIOUS_NONCONVERGENCE = 1, /*!< \hideinitializer Do only if previous procedure \b didn't reach convergence. */
    PREVIOUS_CONVERGENCE    = 2, /*!< \hideinitializer Do only if previous procedure \b did reach converge. */
  };

  friend class SequenceCfg;
  friend class ParameterMgr;
  friend class Maximizer;

  /// @name Constructors & operators
  //@{

	///
	/// Default constructor. Please use the MAXFUN::ParameterMgr interface for creating RunCfg's.
	RunCfg(MaximizationMethodEnum method, int max_iterations);

	///
	/// Copy constructor.
	RunCfg(const RunCfg &);

	///
	/// Assignment operator.
	RunCfg& operator= (const RunCfg &);
  //@}

  /// @name Required information
  //@{

	///
	/// Optimization method
	MaximizationMethodEnum method; // METHOD

	///
	/// Maximum number of iterations; if 0, function value and derivatives will be computed for initial estimates, 
	/// as well as the variance-covariance matrix if var_cov > 0. Please note that if you set max_iterations to
	/// zero, you should set max_method to MAXFUN::direct_without.
	///
	/// Please note that when max_iterations is input into Maxfun, it is first multiplied by the number of
	/// estimated parameters.
	int max_iterations; // MAXIT

  //@}

  /// @name Epsilon information
  //@{

	///
	/// Convergence criterion 1 (\f$\varepsilon_{C1}\f$): maximum relative change in parameter estimates during 
	/// last iteration; default is \f$10^{-3}\f$
	double epsilon1; // EPSC1

	///
	/// Convergence criterion 2 (\f$\varepsilon_{C2}\f$): maximum normalized gradient at last iteration (for 
	/// Newton-Raphson or variable metric methods); default is \f$10^{-15}\f$
	double epsilon2; // EPSC2

	///
	/// Standard \f$\varepsilon_D\f$ for testing and sometimes initializing stepsize factors (\f$\delta\f$'s); 
	/// default is \f$\varepsilon_{C1}\f$ if it is greater than \f$10^{-9}\f$, otherwise \f$10^{-3}\f$.
	double epsilon_delta; // EPSD

	///
	/// Upper bound \f$\varepsilon_T\f$ on truncation error in gradient computation for variable metric methods; 
	/// default is 10 times the square root of \f$\iota\f$.
	double epsilon_trunc; // EPST

  //@}

  /// @name Extra information
  //@{

	///
	/// Variance-covariance matrix calculation option
	///
	/// Default = \c none
	VarCovOptionEnum var_cov; // IXVC

	///
	/// Second derivative approximation calculation option
	///
	/// Default = \c single
	SecDerivOptionEnum second_deriv; // IHIT

	///
	/// Stepsize factor \f$\iota\f$ used in computing first and possibly second derivatives; default is 10 times 
	/// the square root of PRECIS(\f$2^{-45}\f$)
	double stepsize; // YOTA

	///
	/// Each RunCfg can be subjected to a control option. Information on control options
	/// is available in the MAXFUN::maxfun_control_option documentation.
	///
	/// Default = \c unconditional
	MaxfunControlOptionEnum control_option;
  //@}

  protected:

	void  copy(const RunCfg&);

	void TransferContentsToMaxfun(Maxfun& maxfun) const;

	// Dumps the contents of this RunCfg to debug.os;
	void Dump(const DebugCfg& debug) const;

  private:

	RunCfg();
};

} // End namespace MAXFUN
} // End namespace SAGE

#endif
