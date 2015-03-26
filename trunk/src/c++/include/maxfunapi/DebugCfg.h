#ifndef DEBUGCFG_H
#define DEBUGCFG_H

#include <iostream>
#include <fstream>
#include "error/internal_error.h"
#include "app/aparser.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"

namespace SAGE   {
namespace MAXFUN {

/** \class DebugCfg
 *  \brief Identifies which runtime debugging options to enable/disable.
 *
 * \par Purpose
 *
 * The DebugCfg class allows you to specify a variety of debugging options for your
 * maximization. Debugging output is produced at runtime, during maximization.
 *
 * \par Getting started
 *
 * The first step in making use of this class is creating your own instance of it:
 *
 * \code
 * DebugCfg dbg(SAGE::MAXFUN::DebugCfg::NO_DEBUG_INFO);
 * \endcode
 *
 * The parameter passed to the constructor is a MAXFUN::DebugCfg::DebugLevelEnum, and it specifies a predefined
 * debugging template. 
 *
 * \par More options
 * 
 * Once you have created an instance of DebugCfg, however, you can customize
 * various features. Documentation for the different categories of debugging output follows.
 *
 * \par Basic configuration
 *
 * If setReportBasicConfig() is enabled (set to true, that is), then the following information will be reported
 * for each RunCfg, prior to maximization:
 * 
 * \c o sequence template 
 *
 * \c o SequenceCfg::KeepBestRun (if enabled)
 *
 * \par Non-parametric input & output
 *
 * If setReportNonParametricInput() is enabled (set to true, that is), then for each RunCfg the following information
 * will be reported prior to maximization:
 *
 * \c o Maximization method (MAXFUN::RunCfg::MaximizationMethodEnum)
 *
 * \c o Control option (MAXFUN::RunCfg::MaxfunControlOptionEnum)
 *
 * \c o Maximum iterations
 *
 * \c o Epsilon 1
 *
 * \c o Epsilon 2
 *
 * \c o Epsilon delta
 *
 * \c o Epislon trunc
 *
 * \c o Variance-covariance option (MAXFUN::RunCfg::VarCovOptionEnum)
 *
 * \c o Second derivative option (MAXFUN::RunCfg::SecDerivOptionEnum)
 *
 * \c o Stepsize
 *
 * If setReportNonParametricOutput() is enabled (set to true, that is), then for each RunCfg the following information
 * will be reported after maximization:
 *
 * \c o Exit status
 *
 * \c o Final function value
 *
 * \c o Number of iterations performed
 *
 * \c o Number of converged independent parameters
 *
 * \c o Number of dependent parameters
 *
 * \c o Number of estimated parameters
 *
 * \c o Number of independent parameters
 *
 * \c o Number of varying independent parameters
 *
 * \par Parametric input
 *
 * If setReportParametricInput() is enabled (set to true, that is), then the following information
 * will be reported for every parameter, prior to maximization:
 *
 * \c o Group name and parameter name
 *
 * \c o Initial type (MAXFUN::Parameter::ParamTypeEnum)
 *
 * \c o Initial estimate
 *
 * If setReportParametricOutput() is enabled (set to true, that is), then the following information
 * will be reported for every parameter, after maximization:
 *
 * \c o Group name and parameter name
 *
 * \c o Final type (MAXFUN::Parameter::ParamTypeEnum)
 *
 * \c o Final estimate
 *
 * \c o Standard error (if available)
 *
 * \c o Derivative (if available)
 *
 * \c o Covariance information (if available)
 *
 * \par Per-evaluation output
 *
 * If setReportEvaluations() is set to some positive, non-zero value i, then the following information
 * will be produced for every i'th function evaluation:
 *
 * \c o Current parameter estimates
 *
 * \c o Current function evaluation
 *
 * \par Final results
 *
 * It is possible to have the final results reprinted. Please see the documentation on setReportFinalResults()
 * for more information.
 *
 * \par Passing the debug object to SAGE::MAXFUN::Maximize()
 *
 * Remember to pass your DebugCfg object to Maximize()!
 * \code
 * MAXFUN::Maximize(my_cfg, my_info, my_func, my_dbg);
 * \endcode
 */
class DebugCfg
{
  public:

  /// \brief Identifies how much maximization debug output should be produced.
  enum DebugLevelEnum 
  {
    NO_DEBUG_INFO     = 0, /*!< \hideinitializer No debugging output. 
                            */
    BASIC             = 1, /*!< \hideinitializer Information for the final results reported by maximize().
                            * This setting is equivalent to the following code:
                            * \code
                            * setReportBasicConfig  (true);
                            * setReportFinalResults (true);
                            * \endcode
                            */
    PER_RUN           = 2, /*!< \hideinitializer Information for each MAXFUN::RunCfg. 
                            * This setting is equivalent to the following code:
                            * \code
                            * setReportBasicConfig         (true);
                            * setReportNonParametricInput  (true);
                            * setReportParametricInput     (true);
                            * setReportNonParametricOutput (true);
                            * setReportParametricOutput    (true);
                            * setReportFinalResults        (true);
                            * \endcode
                            */
    COMPLETE          = 3  /*!< \hideinitializer Complete information for the entire maximization process.
                            * This setting is equivalent to the following code:
                            * \code
                            * setReportBasicConfig         (true);
                            * setReportNonParametricInput  (true);
                            * setReportParametricInput     (true);
                            * setReportNonParametricOutput (true);
                            * setReportParametricOutput    (true);
                            * setReportFinalResults        (true);
                            * setReportEvaluations         (1);
                            * \endcode
                            */
  };

  friend class SequenceCfg;
  friend class APIMaxFunction;
  friend class Maximizer;
  friend class Parameter;
  friend class RunCfg;
  friend class Results;

  /// @name Constructor / destructor / operators
  //@{

    ///
    /// \param debug_level A preset template that enables/disables a variety of debug options. Documentation available at MAXFUN::DebugCfg::DebugLevelEnum.
    ///
    /// Please note that \b all debugging options are \b disabled by default.
    ///
    /// A MAXFUN::DebugCfg object stores quite a bit of information about which
    /// aspects of maximization should be reported. When you construct the object,
    /// you can use the \c debug_level parameter to indicate a preset debug level.
    /// Information on these debugging presets are available at MAXFUN::DebugCfg::DebugLevelEnum.
    ///
    /// After constructing the object, you can choose to enable/disable any particular
    /// debugging feature. These features are documented below.
    explicit DebugCfg(DebugLevelEnum debug_level = NO_DEBUG_INFO);
    
    ///
    /// Copy constructor.
    /// \param other The object to copy     
    DebugCfg (const DebugCfg & other);
    
    ///
    /// Destructor.
    ~DebugCfg();

    ///
    /// Assignment operator.
    /// \param other The object to copy
    DebugCfg& operator= (const DebugCfg & other);

  //@}

  /// @name Setting the debugging template & output stream
  //@{

    ///
    /// Sets the maximization type.
    /// \param debug_level The debug level requested
    void setType(DebugLevelEnum debug_level);

    ///
    /// Returns the maximization type.
    DebugLevelEnum getType() const;

    ///
    /// Indicates the type of output that this object will produce.
    /// \retval 0 Output will be directed to a user-passed output stream
    /// \retval 1 Output will be directed to a file
    int getOutputType();

    ///
    /// Sets the output stream to which debugging messages will be sent.
    /// \param output_stream The output stream
    void setDebugOutput(ostream & output_stream);
    
    ///
    /// Sets the output file to which debugging messages will be sent.
    ///
    /// Please note that invoking this function \b necessarily creates the indicated file. If there
    /// is to be no debugging output, make sure to check hasAnyReportedOutput() before you invoke
    /// this function.
    ///
    /// \param ofilename The name of the output file
    /// \param append If \c true, debugging messages will be appended to the indicated file. 
    /// If \c false, the file will be cleared before debugging messages are added.
    void setDebugOutput(const string & ofilename, bool append);

    ///
    /// Returns the output stream for directing output.
    ostream & getOutputStream() const;

  //@}

  /// @name Basic configuration
  //@{

    ///
    /// Indicates whether or not this object is configured to produce any output at all.
    /// \retval true There is at least one reported component
    /// \retval false There will be \b not reported output
    bool hasAnyReportedOutput() const;
    
    ///
    /// Reports basic information about the MAXFUN::SequenceCfg object, including sequence template and KeepBestRun.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportBasicConfig(bool report);

    ///
    /// Reports basic information about the MAXFUN::SequenceCfg object, including sequence template and KeepBestRun.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportBasicConfig() const;

    ///
    /// Reports maxfun's non-parametric input from each MAXFUN::RunCfg.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportNonParametricInput(bool report);

    ///
    /// Reports maxfun's non-parametric input from each MAXFUN::RunCfg.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportNonParametricInput() const;

    ///
    /// Reports maxfun non-parametric output from each MAXFUN::RunCfg.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportNonParametricOutput(bool report);

    ///
    /// Reports maxfun non-parametric output from each MAXFUN::RunCfg.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportNonParametricOutput() const;

    ///
    /// Reports maxfun parametric input from each MAXFUN::RunCfg.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportParametricInput(bool report);

    ///
    /// Reports maxfun parametric input from each MAXFUN::RunCfg.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportParametricInput() const;

    ///
    /// Reports maxfun parametric output from each MAXFUN::RunCfg.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportParametricOutput(bool report);

    ///
    /// Reports maxfun parametric output from each MAXFUN::RunCfg.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportParametricOutput() const;

  //@}

  /// @name Per-evaluation output
  //@{

    ///
    /// If setReportEvaluations() is set to some non-zero, postive integer i, then
    /// maximize() will report current parameter estimates and function value for every i-th function evaluation.
    /// \param evaluations Indicates how often to report the function evaluation.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportEvaluations(int evaluations);

    ///
    /// If ReportEvaluations() is set to some non-zero, postive integer i, then
    /// Maximize() will report current parameter estimates and function value for every i-th function evaluation.
    /// \returns x Where the function evaluation will be reported for every x evaluations.
    int getReportEvaluations() const;

  //@}

  /// @name Final results
  //@{

    ///
    /// Reports the final results that maximize() returns.
    /// \param report Indicates whether or not to include this information.
    /// \retval 0 Value was set successfully.
    /// \retval 1 Value was \b not set successfully.
    int setReportFinalResults(bool report);

    ///
    /// Reports the final results that Maximize() returns.
    /// \retval true Option is enabled.
    /// \retval false Option is disabled.
    bool getReportFinalResults() const;

    bool iter_end; // due to JA for dump at end of iteration only

  //@}

  protected:

    // Reports each time a RunCfg is begun.
    // Equivalent to getReportNonParametricInput  () || getReportParametricInput  () ||
    //               getReportNonParametricOutput () || getReportParametricOutput () ||
    //               getReportEvaluations               () 
    bool reportEachRunCfg() const;

    // Equivalent to:
    // getReportNonParametricOutput() || getReportParametricOutput()
    bool reportMaxfunOutput() const;

    // Copy operation (used in copy constructor and operator=)
    void copy(const DebugCfg &);

    // Resets the object.
    void reset();

  private:

    mutable ostream * my_ostream_ptr;

    boost::shared_ptr<ofstream> my_ofstream_shptr;

    DebugLevelEnum my_debug_level;

    bool my_ReportBasicConfig;
    bool my_ReportNonParametricInput;
    bool my_ReportParametricInput;
    bool my_ReportNonParametricOutput;
    bool my_ReportParametricOutput;
    bool my_ReportFinalResults;

    int  my_ReportEvaluations;
};

/// @name Parameter file parsing
//@{

/** \brief Parses a MAXFUN block from a parameter file.
  *
  * It order to faciliate the specification of Maxfun debug output options on a per-analysis
  * basis, this function is made available for parsing a standard Maxfun sub-block in a parameter
  * file.
  *
  * Please see the MAXFUN sub-block parsing section of the main page for more details.
  *
  * \code
  * maxfun
  * {
  *   level=[no_debug_info|basic|per_run|complete]
  * }
  * \endcode
  *
  * \param debug_cfg The DebugCfg object whose options will be set by this function.
  * \param param The LSFBase pointer that points to the MAXUFN block to be parsed.
  * \param errors The errorstream to which error messages will be directed.
  * \retval true MAXFUN Block parsed successfully
  * \retval false MAXFUN Block \b not parsed successfully
  */
bool parseDebugParams(
         DebugCfg     & debug_cfg,
   const LSFBase      * param,
         cerrorstream & errors = sage_cerr);

//@}

//=========================================================
//  INLINE FUNCTIONS
//=========================================================

inline ostream & DebugCfg::getOutputStream() const { if(!my_ostream_ptr) my_ostream_ptr = &cout; return *my_ostream_ptr; }

inline DebugCfg::DebugLevelEnum DebugCfg::getType() const { return my_debug_level; }

}} // End namespace

#endif
