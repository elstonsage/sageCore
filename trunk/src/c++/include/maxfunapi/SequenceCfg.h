#ifndef SEQUENCECFG_H
#define SEQUENCECFG_H

#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/DebugCfg.h"
#include "maxfunapi/RunCfg.h"

namespace SAGE   {
namespace MAXFUN {

/** \class SequenceCfg
 *  \brief Represents a specific configuration for a maximization.
 *
 * \par Purpose
 *
 * This class is used to keep track of maximization options.
 *
 * \par Maximization
 *
 * A complete maximization actually consists of a sequence of individual maximizations, each
 * identified by a single MAXFUN::RunCfg object. 
 *
 * The simplest way to maximize: Invoke MaximizeDefault(), passing it simply the a MAXFUN::MaxfunInfo
 * object and a MaxFunction object.
 *
 * Alternately, you can invoke Maximize() and pass it your own customized MAXFUN::SequenceCfg and
 * MAXFUN::DebugCfg objects:
 *
 * \code
 * SequenceCfg cfg(MAXFUN::max_default);
 * DebugCfg dbg(MAXFUN::debug_none);
 * Maximize(my_info, my_func, cfg, dbg);
 * \endcode
 *
 * The argument in the SequenceCfg() constructor indicates what type of maximization you want. There
 * are a number of predefined maximization templates available; in addition, you can specify that your
 * SequenceCfg instance will represent a user-customized maximization. More information on these templates
 * is available at MAXFUN::SequenceTemplateEnum.
 *
 * If, however, you wish to construct your own maximization, you can do so with the addRunCfg()
 * and getLatestRunCfg() functions. Each call to addRunCfg() adds one more individual
 * maximization to the sequence. getLatestRunCfg() can be used to tweak specific aspects of
 * the RunCfg at the end of the current list.
 *
 * For instance, let's say you wanted to start with a direct search method, with 4 iterations, and
 * epsilon1 = 2. The following code would accomplish this:
 *
 * \code
 *  my_SequenceCfg.addRunCfg(direct_search, 4);
 *  my_SequenceCfg.getLatestRunCfg().epsilon1 = 2;
 * \endcode
 *
 * If you want to tweak aspects of a maximization template, you can do that too. First, remember
 * to create a MAXFUN::SequenceCfg object with some predefined template:
 *
 * \code
 *  SequenceCfg my_SequenceCfg(MAXFUN::max_default);
 * \endcode
 *
 * Now you can edit individual maximizations with the getRunCfg() function. Recall that the
 * maximization template \c max_default creates a sequence of MAXFUN::RunCfg objects. With the
 * getRunCfg() function, you can fetch individual maximizations created by the template and
 * then edit them. For instance, if you wanted epsilon1 changed to 1e-2 for the first individual
 * maximization, you could use the following code:
 *
 * \code
 *  my_SequenceCfg.getRunCfg(0).epsilon1 = 1e-2;
 * \endcode
 *
 * Please note that there are a number of options available for each MAXFUN::RunCfg. You can consult
 * the MAXFUN::RunCfg documentation for information on these features.
 *
 */
class SequenceCfg
{
public:

	///
	/// Identifies a preset maximization template.
	///
	/// Since a MAXFUN::SequenceCfg represents a sequence of MAXFUN::RunCfg's, there are a variety
	/// of ways to configure that sequence. The standard method is \c DEFAULT_MAXIMIZATION, which is detailed below.
	///
	/// Alternately, you can construct your own sequence of MAXFUN::RunCfg's. Please see the MAXFUN::SequenceCfg
	/// documentation for information on how to do this.
	enum SequenceTemplateEnum 
	{ 
	  USER_DEFINED = 0,
	  /*!< \hideinitializer 
	   * User-defined maximization technique.
	   */

	  DEFAULT_MAXIMIZATION = 1 
	  /*!< \hideinitializer 
	   * The DEFAULT_MAXIMIZATION method is equivalent to the following user-defined method:
	   * \code
	   *
	   * addRunCfg(RunCfg::DIRECT_WITHOUT, 1);
	   * getLatestRunCfg().epsilon1 = 1e-3;
	   * getLatestRunCfg().epsilon2 = 1e-12;
	   *
	   * addRunCfg(RunCfg::VAR_METRIC_IDENTITY, 20);
	   * getLatestRunCfg().epsilon1       = 1e-3;
	   * getLatestRunCfg().epsilon2       = 1e-12;
	   * getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;
	   * getLatestRunCfg().var_cov        = RunCfg::FINAL;
	   *
	   * addRunCfg(RunCfg::DIRECT_WITHOUT, 50);
	   * getLatestRunCfg().epsilon1       = 1e-4;
	   * getLatestRunCfg().epsilon2       = 1e-12;
	   * getLatestRunCfg().control_option = RunCfg::PREVIOUS_NONCONVERGENCE;
	   *
	   * addRunCfg(RunCfg::VAR_METRIC_ESTIMATE, 20);
	   * getLatestRunCfg().epsilon1 = 1e-12;
	   * getLatestRunCfg().epsilon2 = 1e-12;
	   * getLatestRunCfg().var_cov  = RunCfg::FINAL;
	   *
	   * \endcode
	   */
	};

	friend class Maximizer;

	/// @name Constructor
	//@{
		///
		/// This constructor will initialize the configuration to a template.
		/// \param sequence_template Indicates which template to use.
		/// \param sequence_name The name given to this sequence.
		/// \param function_name The name given to this function (ie: "likelihood")
		///
		/// Please consult the documentation on SequenceCfg::SequenceTemplateEnum for information on
		/// configuration templates.
		SequenceCfg(SequenceTemplateEnum sequence_template = DEFAULT_MAXIMIZATION, 
		            string sequence_name = "Default analysis", 
		            string function_name = "function value");

		///
		/// Copy constructor
		SequenceCfg(const SequenceCfg &);

		///
		/// Assignment operator
		SequenceCfg & operator= (const SequenceCfg &);
	//@}

	/// @name Optional maximization settings
	//@{

		///
		/// Returns the KeepBestRun status. Please see SetKeepBestRun() for more information.
		bool getKeepBestRun() const;

		///
		/// By default, Maximize() will return the results for the last maximization
		/// step. If, however, you want to force Maximize() to
		/// return the last stage maximization step that actually converged, you can set set 
		/// setKeepBestRun() to \c true:
		/// \code setKeepBestRun(true); \endcode
		/// \param keep Indicates whether or not to keep the best run.
		/// \retval 0 Value was set successfully.
		/// \retval 1 Value was not set successfully.
		int setKeepBestRun(bool keep);

	//@}

	/// @name Informative functions
	//@{
	
	        ///
	        /// Assigns the sequence name.
	        void setSequenceName(const std::string & name);
	        
		///
		/// Returns the name assigned to this sequence (in the constructor).
		const std::string & getSequenceName() const;

		///
		/// Assigns the function name.
		void setFunctionName(const std::string & name);
		
		///
		/// Returns the name assigned to the function associated with this analysis (in the constructor).
		const std::string & getFunctionName() const;

		///
		/// Returns the configuration template number for this object.
		SequenceTemplateEnum getSequenceTemplate();

	//@}

	/// @name User-defined maximization functions...
	//@{
		///
		/// Use addRunCfg() to add a run maxfun.
		/// \param method Maximization method (see MAXFUN::RunCfg::method for information)
		/// \param max_iterations (see MAXFUN::RunCfg::max_iterations for information)
		/// \retval 0 Success
		/// \retval 1 Failure
		int addRunCfg(RunCfg::MaximizationMethodEnum method, int max_iterations);

		///
		/// adds another MAXFUN::RunCfg to the end of the RunCfg sequence that is identical
		/// to the most recently added MAXFUN::RunCfg.
		int duplicateLatestRunCfg();

		///
		/// Returns a reference to the most recently added RunCfg.
		RunCfg & getLatestRunCfg();

		///
		/// Returns a reference to i-th RunCfg.
		/// \param i Index number of requested RunCfg. Please note that indices start at 0, not at 1.
		RunCfg & getRunCfg(int i);

	//@}

	protected:

		// Copys the contents of one SequenceCfg instance to another (used in copy constructor
		// and operator=).
		void copy(const SequenceCfg &);

		// Dumps information about this SequenceCfg object to debug.os.
		void dump(const DebugCfg & debug) const;

		// Sets up the SequenceCfg object according to a preset template.
		void setSequenceTemplate(SequenceTemplateEnum sequence_template);

		// Returns the vector of RunCfg's.
		const vector<RunCfg> & getRunCfgs() const;

		// Returns the vector of RunCfg's.
		vector<RunCfg> & getRunCfgs();

	private:
		string               my_sequence_name;
                string               my_function_name;
		SequenceTemplateEnum my_sequence_template;
		vector<RunCfg>       my_RunCfgs;
		bool                 my_KeepBestRun;
};

}} // End namespace

#endif
