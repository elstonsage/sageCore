#ifndef RESULTS_H
#define RESULTS_H

#include "maxfun/maxfun.h"
#include "maxfunapi/CovarianceMatrix.h"
#include "maxfunapi/DebugCfg.h"
#include "maxfunapi/ParameterMgr.h"
#include "maxfunapi/OutputFormatter.h"

namespace SAGE   {
namespace MAXFUN {

/** \class Results
 *  \brief This class is used to keep track of maximization results.
 *
 *
 * \par Purpose
 *
 * Results stores the entire set of maximization results, including parameter
 * information, maximization status, and the variance-covariance matrix.
 *
 * \par How to debug your results
 *
 * You can use the checkSelf() function to see if there were any problems with
 * your maximization. If you invoke checkSelf() with no arguments, then checkSelf()
 * will simply return 0 if there were no errors detected, and 1 if were were errors.
 * Please see the documentation on checkSelf() for more information.
 *
 * \par Basic information
 *
 * The basic information includes the maximization's name, initial validity, final
 * function value, convergence status, exit flag, and iterations completed.
 *
 * \par Variance-covariance matrix
 *
 * If the maximization was configured to calculate a variance-covariance matrix (a default
 * option for MAXFUN::max_default), then Maximize() will have attempted to calculate
 * the matrix. You can generally use the SAGE::FormatMatrix() function to report the matrix,
 * but you can query the Results object directly for this information as well.
 *
 * First, make sure the matrix is available by checking the getCovMatrixStatus() flag.
 *
 * Then, use getParameterMgr().getCovarianceMatrix() to get access to the matrix's values.
 *
 * \par Additional information
 *
 * There are a variety of other informative functions available. Please consult the documentation
 * for those functions for more information.
 */
class Results
{
  friend class Maximizer;
  friend class OutputFormatter;
  MEM_FRIEND(SAGE::MAXFUN::Results);

  public:

  /// @name Constructors and operators:
  //@{

    ///
    /// Constructor #1
    Results();
    
    ///
    /// Constructor #2
    /// \param data The Maxfun_Data instance whose results this object will represent
    /// \param debug The DebugCfg instance that describes which debugging features to
    /// enable when importing the data from the Maxfun_Data instance.
    Results(const Maxfun_Data & data, const DebugCfg & debug);


    ///
    /// Copy constructor
    Results(const Results &);

    ///
    /// Destructor
    ~Results();

    ///
    /// Assignment operator
    Results& operator= (const Results &);

  //@}

  /// @name Debugging functions
  //@{
    ///
    /// Performs a check of the maximization results.
    /// \param debug Contains information about which debugging output to produce.
    /// \retval 0 No error condition detected.
    /// \retval 1 Error condition detected. It is recommended that you run 
    /// \code checkSelf(MAXFUN::DebugCfg(MAXFUN::DebugCfg::COMPLETE)); \endcode to examine the specific 
    /// error condition.
    int checkSelf(const DebugCfg & debug, const Maxfun_Data & data) const;

  //@}

  /// @name Basic information
  //@{

         double score; // due to JA
         vector<std::string> flat_dir; // due to JA
    ///
    /// Returns the name associated with this sequence.
    ///
    /// The name is assigned by the MAXFUN::SequenceCfg object that is used in invoking maximize().
    const std::string & getSequenceName() const;

    ///
    /// Sets the name associated with this sequence.
    /// \param sequence_name The name for this sequence
    /// \retval 0 Value was set successfully
    /// \retval 1 Value was not set successfully
    int setSequenceName(const std::string & sequence_name);

    ///
    /// Returns the name for this sequence's function.
    ///
    /// The name is assigned by the MAXFUN::SequenceCfg object that is used in invoking maximize().
    const std::string & getFunctionName() const;

    ///
    /// Sets the name associated with this function.
    /// \param function_name The name for this function
    /// \retval 0 Value was set successfully
    /// \retval 1 Value was not set successfully
    int setFunctionName(const std::string & function_name);

    ///
    /// This the RunCfg corresponding to this Results object was skipped (because of
    /// the RunCfg::control_option), then WasSkipped will return true.
    /// \retval true This maximization was skipped.
    /// \retval false This maximization was \b not skipped.
    bool getWasSkipped() const;


    bool compscorestat(); // due to JA, to facilitate score test implementation

    // \param was_skipped Indicates whether or not this analysis was skipped.
    // \retval 0 Value was set successfully.
    // \retval 1 Value was not set successfully.
    int setWasSkipped(bool was_skipped);

    ///
    /// Identifies whether or not this maximization had a valid initial
    /// function evaluation.
    /// \retval true This maximization \b did have a valid initial value, and maximization continued.
    /// \retval false This maximization did \b not have a valid initial value, so maximization was aborted.
    bool getValidInitialValue() const;

    // \param valid Indicates whether or not this analysis was valid.
    // \retval 0 Value was set successfully.
    // \retval 1 Value was not set successfully.
    int setValidInitialValue(bool valid);

    ///
    /// Identifies whether or not this maximization converged.
    /// \retval true This maximization \b did converge.
    /// \retval false This maximization did \b not converge.
    bool getConverged() const;

    ///
    /// Returns the final exit flag describing conditions of termination of MAXFUN:
    ///
    /// \retval 0 No errors found, but zero iterations requested
    ///
    /// \retval 1 Convergence by criterion 1
    ///
    /// \retval 2 Convergence by criterion 2
    ///
    /// \retval 3 Convergence by criterion 3
    ///
    /// \retval 4 Reached maximum number of iterations without convergence
    ///
    /// \retval 5 Accumulation of round-off errors or boundary problem
    /// prevents further progress (in variable metric methods)
    ///
    /// \retval 6 Computed search direction not upward (in variable metric methods)
    ///
    /// \retval 7 All significant digits lost in optimal conditioning during matrix update
    /// (in variable metric methods)
    ///
    /// \retval 8 Gradient could not be computed (because too close to a bound or other constraint)
    ///
    /// \retval 9 Variance-covariance matrix could not be computed (in Newton-Raphson methods); either
    /// second partial derivatives could not be computed (because too close to a bound or other 
    /// constraint) or matrix of second partials could not be inverted.
    ///
    /// \retval 10 All independent parameters converged to bounds.
    ///
    /// \retval 11 Too many parameters to estimate using derivatives (N(i) > NPV)
    ///
    /// \retval 12 Initial parameter estimates not in domain of the function to be maximized.
    ///
    /// \retval 13 Control input invalid.
    ///
    int getExitFlag() const;

    ///
    /// Returns the final function value.
    double getFinalFunctionValue() const;

    ///
    /// If for any reason you need to adjust the final function value, you can use
    /// setFinalFunctionValue() to do so.
    /// \param val The new final function value.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setFinalFunctionValue(double val);

    ///
    /// Returns the number of iterations completed.
    int getIterations() const;

  //@}

  /// @name Variance-covariance matrix information
  //@{

    ///
    /// Returns the variance-covariance matrix status.
    /// \retval 0 No problem
    /// \retval 1 Round-off error in Hessian matrix
    /// \retval 2 Variance-covariance matrix has negative values on the diagonal (round-off error may also be a problem)
    /// \retval 3 Hessian matrix cannot be inverted (round-off error may also be a problem)
    /// \retval 4 Hessian matrix could not be computed (too close to a bounse or other constraint)
    /// \retval 5 Not attempted
    int getCovMatrixStatus() const;

    ///
    /// Returns a non-const reference to the MAXFUN::CovarianceMatrix.
    CovarianceMatrix & getCovarianceMatrix();

    ///
    /// Returns a const reference to the MAXFUN::CovarianceMatrix.
    const CovarianceMatrix & getCovarianceMatrix() const;

  //@}

  /// @name Parameter information
  //@{

    ///
    /// Once the maximization has been performed and the Results object generated,
    /// the ParameterMgr object is copied into the Results object. If you want to
    /// perform any post-maximization modifications to the final parameter estimates, you
    /// should use the getParameterMgr() function to do so.
    ParameterMgr & getParameterMgr();

    ///
    /// Returns the ParameterMgr object.
    const ParameterMgr & getParameterMgr() const;

    ///
    /// Returns the number of independent parameters that have converged to a bound.
    /// (Those for which FinalType = 5-8)
    int getNumOfConvergedIndependentParams() const;
  
    ///
    /// Returns the number of dependent parameters.
    /// (Those for which FinalType = 3)
    int getNumOfDependentParams() const;

    ///
    /// Returns the number of estimated parameters.
    /// (Those for which FinalType = 1-3, 5-10)
    int getNumOfEstimatedParams() const;

    ///
    /// Returns the number of independent parameters.
    /// (Those for which FinalType = 1-2, 5-10)
    int getNumOfIndependentParams() const;

    ///
    /// Returns the number of varying independent parameters.
    /// (Those for which FinalType = 1-2)
    int getNumOfVaryingIndependentParams() const;

  //@}

  /// @name Additional information
  //@{

    ///
    /// Returns the maximum absolute difference between function values considered during last direct 
    /// search iteration and function value at the start of the iteration.
    double getMaximumDifference() const;

    ///
    /// Returns the maximum (over i) magnitude of change \f$|\Delta \Theta_i|\f$ in last iteration.
    double getERM() const;

    ///
    /// Returns the change in function value \f$f(\vec{\Theta}) - f(\vec{\Theta^0})\f$ in last iteration 
    double getChangeInFunctionValue() const;

    ///
    /// Returns the value of function for previous parameter estimates \f$f(\vec{\Theta^0})\f$.
    double getPreviousFunctionValue() const;

    ///
    /// Returns the norm of gradient \f$ \vec{g}\prime\vec{g} \f$
    double getGradientNorm() const;

    ///
    /// Returns the normalized gradient \f$ \vec{p}\prime\vec{g} \f$; 
    /// for Newton-Raphson or variable metric methods only.
    double getNewtonGradientNorm() const;

    ///
    /// Returns the first derivative appromixation indicator.
    /// \retval 1 Forward difference used
    /// \retval 2 Central difference used
    int getFirstDerivApproxIndic() const;

    ///
    /// Returns the indicator of age of gradient vector.
    /// \retval -1 No valid gradient vector
    /// \retval 0 Have gradient \f$ \vec{g}(\vec{\Theta}) \f$ corresponding to final parameter estimates.
    int getAgeOfGradientVector() const;

    ///
    /// Returns the gradient status.
    /// \retval 0 No problem
    /// \retval 1 Gradient could not be computed (too close to bound or other constraint)
    /// \retval 2 Gradient could not be attempted
    int getGradientStatusFlag() const;

    ///
    /// Returns the termination at constraint status.
    /// \retval 0 Implied boundary not involved in termination.
    /// \retval 1 Terminating at or near bound of dependent parameter or boundary implied by other constraint.
    int getTerminationAtConstraint() const;

    ///
    /// Returns the age of variance-covariance matrix status.
    /// \retval -1 No valid variance-covariance matrix
    /// \retval >=0 Number of iterations since V computed
    int getAgeOfCovMatrix() const;

    ///
    /// Returns number of times \f$2^{N_i}\f$-trial search improved estimates (for direct search methods only)
    /// \retval 0 - 2 As indicated above
    /// \retval 3 The search made a significant improvement twice, but was not tried a third time
    int getNumOfTrialSearch() const;

    ///
    /// Returns the step size of change in established direction.
    /// For Newton-Raphson or variable metric methods.
    double getStepSizeChange() const;

    int    getAugvStatus()   const; 
  //@}


//    const Maxfun_Data & getMaxfunData() const;

  protected:


    // Copy operation
    void copy (const Results &);

    // Update operation
    void inputMaxfunData(const Maxfun_Data &, const DebugCfg & debug);

    // Update operation #2
    void inputParameterMgr(const ParameterMgr &);


  private:
    ParameterMgr     my_ParameterMgr;
    CovarianceMatrix my_CovarianceMatrix;

    std::string my_sequence_name;
    std::string my_function_name;

    bool   my_ValidInitialValue;
    bool   my_WasSkipped;
    int    my_ExitFlag;                             // LFL
    double my_FinalFunctionValue;                   // F
    int    my_Iterations;                           // IT
    int    my_CovMatrixStatus;                      // IVFL
    int    my_NumOfBoundConvergedIndependentParams; // NB
    int    my_NumOfDependentParams;                 // ND
    int    my_NumOfEstimatedParams;                 // NE
    int    my_NumOfIndependentParams;               // NI
    int    my_NumOfVaryingIndependentParams;        // NV
    int    my_NumOfTerms;                           // NT
    double my_MaximumDifference;                    // DIXMAX
    double my_ERM;                                  // ERM
    double my_ChangeInFunctionValue;                // FCH
    double my_PreviousFunctionValue;                // FPR
    double my_GradientNorm;                         // GTG
    double my_NewtonGradientNorm;                   // PTG
    int    my_FirstDerivApproxIndic;                // IDIF
    int    my_AgeOfGradientVector;                  // IGAGE
    int    my_GradientStatusFlag;                   // IGFL
    int    my_TerminationAtConstraint;              // IMPBND
    int    my_AgeOfCovMatrix;                       // IVAGE
    int    my_NumOfTrialSearch;                     // NSURF2
    double my_StepSizeChange;                       // TSTEP

    int    my_AugvStatus;                            // IRETURN (in augv)
};

}} // End namespace

MEM_COUNT_BEGIN(SAGE::MAXFUN::Results)
{
  return get_mem(t.my_ParameterMgr) +
         get_mem(t.my_CovarianceMatrix) +
         get_mem(t.my_sequence_name) +
         get_mem(t.my_function_name) +
         get_mem(t.flat_dir) +
         sizeof(SAGE::MAXFUN::Results);
}
MEM_COUNT_END

#include "maxfunapi/Results.ipp"

#endif
