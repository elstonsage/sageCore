#ifndef PARAMETERMGR_H
#define PARAMETERMGR_H

#include <list>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <utility>
#include <algorithm>
#include <sstream>
#include <exception>
#include "numerics/cephes.h"
#include "numerics/functions.h"
#include "util/get_mem.h"
#include "maxfun/maxfun.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/DebugCfg.h"
#include "maxfunapi/Parameter.h"
#include "maxfunapi/ParameterInput.h"
#include "maxfunapi/Submodel.h"

namespace SAGE   {
namespace MAXFUN {

typedef boost::shared_ptr<const ParamCalculator> ParamCalculatorShCstPtr;
typedef boost::shared_ptr<ParamCalculator>       ParamCalculatorShPtr;

/// \brief Calculates a dependent parameter
/// For dependent parameters, it may be useful to calculate the parameter with
/// a functor (rather than place the calculation code directly in your program).
///
/// With ParamCalculator, you can! 
class ParamCalculator
{
public:
  virtual ~ParamCalculator() { }
  virtual double calculateParam (const ParameterMgr * mgr) const = 0;
};

/// Adds together all the parameters whose indices are given in this class' constructor.
class AdditiveParamCalculator : public ParamCalculator
{
public:
  AdditiveParamCalculator(const vector<int> & ids);
  virtual double calculateParam (const ParameterMgr * mgr) const;
private:
  vector<int> my_ids;
};

//==============================================================================
//                 Transformer stuff 
//===============================================================================

/** \class Transformer
 *
 * Parameter transformation is a very neat feature of the Maxfun API!
 *
 * Let's say you've got two parameters, where one parameter is dependent on the other.
 * For instance, "stdev" is estimated, while "var" (where var = stdev * stdev) is simply
 * a dependent parameter. The standard way of accomplishing the "var = stdev*stdev" portion
 * would be to include the calculation somewhere in your update_bounds() dependent parameter
 * calculation portion.
 *
 * Well, that's fine and dandy. It will work perfectly well. But what if you wanted to be able
 * to "plug in" a standard independent-to-dependent transformation \b without adding it
 * to your dependent parameter update code? The solution is to use a Transformer. A Transformer
 * is a simple chunk of code that will convert the estimated parameter ("maximization scale")
 * to the dependent parameter ("reported scale").
 *
 * Consider the following code snippet:
 *
 * \code
 * class SpecialTransformer : public Transformer
 * {
 * public:
 *   virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const { return sqrt(val); }
 *   virtual double transformToReportedScale     (const ParameterMgr * mgr, double val) const { return val * val; }
 * };
 *
 * ...
 *    
 * // Add the maximization-scale parameter
 * my_mgr.addParameter("global", "stdev", MAXFUN::Parameter::INDEPENDENT_FUNCTIONAL, 0.0, -MAXFUN::MF_INFINITY, MAXFUN::MF_INFINITY);
 * 
 * // Add the reported-scale parameter
 * my_mgr.addParameter("reported global", "var", MAXFUN::Parameter::DEPENDENT, 0.0, 0.0, MAXFUN::MF_INFINITY);
 *
 * // Create your transformer!
 * MAXFUN::TransformerShCstPtr tr(new MAXFUN::SpecialTransformer());
 *              
 * // Associate the transformer with the indicated parameters
 * my_mgr.addTransformer(tr, "global", "stdev", "reported global", "var");
 * \endcode
 */
class Transformer
{
  public:

    // Required to make compiler happy.
    virtual ~Transformer() { }

    Transformer            () {};
    Transformer            (const Transformer &);
    Transformer& operator= (const Transformer &);

    virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const = 0;
    virtual double transformToReportedScale    (const ParameterMgr * mgr, double val) const = 0;

    virtual double getLowerBound() const { return -MF_INFINITY; }
    virtual double getUpperBound() const { return MF_INFINITY; }
};

class SimpleTransformer : public Transformer
{
public:
  virtual double transformToMaximizationScale (const ParameterMgr * mgr, double val) const { return val; }
  virtual double transformToReportedScale     (const ParameterMgr * mgr, double val) const { return val; }
};

typedef boost::shared_ptr<const Transformer> TransformerShCstPtr;
typedef boost::shared_ptr<Transformer> TransformerShPtr;

struct TransformerInfo
{
  TransformerInfo();
  TransformerInfo(const TransformerInfo & other);
  TransformerInfo& operator=(const TransformerInfo & other);

  int reported_param_id;
  int maximization_param_id;
  TransformerShCstPtr transformer;
};

inline TransformerInfo::TransformerInfo() {}

inline TransformerInfo::TransformerInfo(const TransformerInfo & other)
{ reported_param_id = other.reported_param_id; maximization_param_id = other.maximization_param_id; transformer = other.transformer; }

inline TransformerInfo& TransformerInfo::operator=(const TransformerInfo & other)
{ reported_param_id = other.reported_param_id; maximization_param_id = other.maximization_param_id; transformer = other.transformer; return *this; }

//===============================================================================
//                 NewSubmodel stuff 
//===============================================================================

template<class SM>
struct SMType
{
  typedef SM sm_type;
};

typedef boost::shared_ptr<NewSubmodel> NewSubmodelShPtr;
typedef boost::shared_ptr<const NewSubmodel> NewSubmodelShCstPtr;

class NewSubmodel
{
  friend class ParameterMgr;

public:

  /// @name Constructor / operators
  //@{

    inline NewSubmodel(const string & name, cerrorstream & errors = sage_cerr);
    inline NewSubmodel(const NewSubmodel & other);
    inline NewSubmodel& operator=(const NewSubmodel & other);
    virtual inline ~NewSubmodel() { }
    
  //@}
  
  /// @name Required virtual interface
  //@{
  
    ///
    /// Updates any relevant internal state from ParameterMgr.
    ///
    /// Note: You are \b required to define this function in a derived class.
    /// \retval 0 Internal state was updated succesfully.
    /// \retval 1 internal state was \b not updated successfully.
    virtual int update() = 0;

    virtual NewSubmodelShPtr clone() = 0;

  //@}
  
  
  /// @name Optional virtual interface
  //@{
  
    ///
    /// If for any reason you elected \b not to resize the my_parameters vector
    /// in your initial setup functions, you can do so within the finalizeConfiguration()
    /// function. This function will be called immediately prior to actually \b adding 
    /// the contents of my_parameters to Maxfun.
    /// \retval 0 Configuration was finalized successfully.
    /// \retval 1 Configuration was \b not finalized successfully.
    virtual inline int finalizeConfiguration();
                                                                
    /**
      *  If you want to make use of any of MAXFUN::Parameter's advanced features
      *  (output formatting, p-value selection, etc.), you should do so within the confines
      *  of the setAdvancedParameterOptions() function. Within this function, you can use
      *  the getParameterMgr()->getParameter() function to access a parameter's full interface.
      *  For instance, if you want to exclude a parameter's p-value from the final output, you
      *  could use the following code:
      *
      *  \code  
      *  void mysubmodel::setAdvancedParameterOptions()
      *  {
      *    getParameterMgr()->getParameter("submode1", "X1").IncludepValue() = false;
      *  }
      *  \endcode
      */
    virtual inline void setAdvancedParameterOptions();
                                                                                                                                              
  //@}
                                                                                                                                                
  /// @name Basic info
  //@{

    ///
    /// Sets the name of the NewSubmodel.
    inline void setName(const string & name);

    ///
    /// Returns the name of the submodel.
    inline const string & getName() const;
    
    ///
    /// Returns the parameter mgr to which this object is bound.
    inline const ParameterMgr & getParameterMgr() const;

    ///
    /// Returns whether or not this object has been finalized (that is, have its
    /// parameters been added to a ParameterMgr for maximization).
    inline bool isFinalized() const;

  //@}

  /// @name Runtime maximization functions
  //@{

    ///
    /// Within run-time maximization, you should fetch parameter values by using
    /// the getParam() function.
    /// \param local_id The index (within my_params) of the parameter you want.
    /// \returns A double reference to the parameter's current estimate.
    inline double & getParam(int local_id);

    ///
    /// Const version of getParam().
    inline double getParam(int local_id) const;

  //@}

protected:

    ///
    /// my_params is a vector of MAXFUN::ParameterInput's, designed to help
    /// you to set up initial options for your parameters.
    vector<ParameterInput> my_parameters;

    /// \internal
    ///
    /// Stream for reporting errors.
    mutable cerrorstream  my_errors;
                

private:
    string my_name;

    bool my_finalized;

    /// \internal
    /// The ParameterMgr that "owns" this submodel.
    ParameterMgr * my_mgr;

    /// \internal
    /// Points this submodel to the given ParameterMgr.
    inline void setParameterMgr(ParameterMgr * mgr);

    /// \internal
    /// Adds the necessary parameters to ParameterMgr.
    /// \retval 0 Parameters added successfully.
    /// \retval !=0 Error code (unsuccessful)
    inline int addParametersToMgr();

    /// \internal
    /// Sets the finalized flag.
    inline void setFinalized(bool finalized);
};


/** \class ParameterMgr
 *  \brief Keeps track of a set of Maxfun parameters.
 *
 * \par Purpose
 *
 * ParameterMgr acts as a parameter manager, allowing you to add & setup parameters as well as fetch & set their
 * values during maximization.
 *
 * \par Parameters
 *
 * Every parameter that you want to use \b must be added first to the ParameterMgr object. You can do this with the
 * AddParameter() function.
 *
 * \par Sub-models
 *
 * Every submodel that you want to use \b must be added first to the ParameterMgr object. You can do this with the
 * AddSubModel() function.
 *
 * \par Parameter groups
 *
 * Every parameter belongs to one and only group. These groups are identified by name.
 *
 * Although a parameter group is automatically added whenever you invoke AddParameter(), you can also explicitly
 * add empty groups with the addGroup() function. 
 *
 * It is possible to iterate across the parameters within a single group. The ParamBegin() and ParamEnd()
 * iterators are available for this purpose. Each ParamBegin() or ParamEnd() function returns either a
 * ParameterIterator or a ParameterConstIterator. This iterator can be dereferenced to return a
 * MAXFUN::Parameter reference.
 *
 * When you invoke ParamBegin() or ParamEnd(), you can either pass it no arguments (in which case it will
 * iterate across \b all parameters), or the name of a group (in which case it will iterate across the parameters
 * in that group).
 *
 * \par Dependent parameters
 *
 * If you are using any dependent parameters in your maximization, you will have to calculate the dependent's
 * value in your update_bounds() function. 
 *
 * \code
 * int my_func::update_bounds(vector<double> & params)
 * {
 *   my_ParameterMgr("group1", "dependent_X") = 2.0 * my_ParameterMgr("group1", "independent_Y");
 * }
 * \endcode
 *
 */
class ParameterMgr
{
public:

  friend class APIMaxFunction;
  friend class Function;
  friend class OutputFormatter;
  friend class Results;
  friend class Submodel;
  friend class Maximizer;
  MEM_FRIEND(SAGE::MAXFUN::ParameterMgr);

  /// @name Constructors & operators
  //@{
    ///
    /// Default constructor.
    ParameterMgr();

    ///
    /// Copy constructor.
    ///
    /// WARNING:  Cannot be used if either ParameterMgr object
    /// involved is in use (see isInUse()).  If this is
    /// attempted, the program will exit.
    /// \param other The ParameterMgr whose contents will be copied into this object.
    ParameterMgr(const ParameterMgr & other);

    ///
    /// Assignment operator.
    ///
    /// WARNING:  Cannot be used if either ParameterMgr object
    /// involved is in use (see isInUse()).  If this is
    /// attempted, the program will exit.
    /// \param other The ParameterMgr whose contents will be copied into this object.
    ParameterMgr& operator= (const ParameterMgr & other);

  //@}

  /// @name Destructor
  //@{

    ///
    /// Destructor.
    ~ParameterMgr();

  //@}

  /// @name Initial setup functions
  //@{

    ///
    /// Reset() clears all parameter and group entries from the ParameterMgr instance.
    ///
    /// WARNING:  Cannot be used if the ParameterMgr is in use
    /// (see isInUse()).  If this is attempted, the program will
    /// exit.
    void reset();

    ///
    /** Adds a parameter to the vector of Maxfun parameters.
     *
     * \param group_name The name of the group to which the parameter belongs. Note: 'ALL' is a reserved 
     * group name.
     * \param param_name The name of the parameter.
     * \param initial_type The initial type for the parameter.
     * \param initial_estimate The initial estimate for the parameter.
     * \param lower_bound The lower bound for estimating the parameter.
     * \param upper_bound The upper bound for estimating the parameter.
     *
     * Please note that if \c group_name == \c param_name, FormatEstimates() will interpret this
     * to mean that it should \b not include a group_name header for this parameter.
     *
     * \b Example: \code AddParameter("global", "beta", MAXFUN::Parameter::INDEPENDENT, 1.0, MAXFUN::-MF_INFINITY, MAXFUN::MF_INFINITY); \endcode
     *
     * \retval -1 Failed to add parameter
     * \retval >=0 The id number of the parameter.
     */
    int addParameter(
      string                   group_name,
      string                   param_name, 
      Parameter::ParamTypeEnum initial_type     =  MAXFUN::Parameter::INDEPENDENT, 
      double                   initial_estimate =  0.0, 
      double                   lower_bound      = -MAXFUN::MF_INFINITY, 
      double                   upper_bound      =  MAXFUN::MF_INFINITY,
      double                   init_stepsize    =  0.1);
      
    Parameter&  addParameterAlt(
      string                   group_name,
      string                   param_name, 
      Parameter::ParamTypeEnum initial_type     =  MAXFUN::Parameter::INDEPENDENT, 
      double                   initial_estimate =  0.0, 
      double                   lower_bound      = -MAXFUN::MF_INFINITY, 
      double                   upper_bound      =  MAXFUN::MF_INFINITY,
      double                   init_stepsize    =  0.1);      

    ///
    /// Alternately, you can use the version of AddParameter() that takes a MAXFUN::ParameterInput 
    /// reference. This structure is designed to help you organize only the required information for
    /// adding a parameter. It is used primarily in the submodels; please note, however, that the
    /// submodel's particular use of MAXFUN:ParameterInput's obviates the need to explicitly invoke
    /// addParameter().
    /// \param param The ParameterInput to be added.
    /// \retval -1 Failed to add parameter
    /// \retval >=0 The id number of the parameter.
    int addParameter(ParameterInput& param);
    Parameter&  addParameterAlt(ParameterInput& param);    
                
    ///
    /// Adds a ParamCalculator.
    /// \param idx The index of the parameter to be calculated
    /// \param calculator The ParamCalculator which will calculate this param's value.
    int addParamCalculator(
      int idx,
      ParamCalculatorShPtr calculator);

    ///
    /// Adds the indicated transformer to the named parameter.
    int addTransformer(
      TransformerShCstPtr transformer, 
      const string& maximization_group_name, 
      const string& maximization_param_name,
      const string& reported_group_name,
      const string& reported_param_name);
                
    ///
    /// A ParameterMgr class is 'in use' when:
    /// - It has any active submodel links (hasLinkedSubmodels).
    bool isInUse() const;

    ///
    /// Adds an empty group of parameters.
    /// \param group_name The name of the group to be added.
    void addGroup(string group_name) const;

    ///
    /// Prints information about the current configuration to stdout.
    void dumpConfiguration() const;

  //@}

  /// @name Submodels
  //@{

    ///
    /// Adds a submodel to the ParameterMgr object (and consequently to the parameter vector).
    /// Please note that for every submodel you use, you \b must explicitly add that submodel
    /// to your ParameterMgr object prior to function evaluation and maximization.
    ///
    /// Please note that there is a very specific order of events necessary for using a submodel.
    /// To use a submodel, you must do the following:
    ///
    /// \c 1 Invoke AddSubModel() immediately prior to invoking Maximize().
    ///
    /// \c 2 Invoke Maximize()
    ///
    /// \c 3 Invoke removeSubmodels() immediately thereafter.
    ///
    /// \b Example:
    ///
    /// \code
    /// my_maxfun_info.addSubModel(&my_submodel);
    /// Results my_results = maximizeDefault(my_maxfun_info, my_max_function);
    /// my_maxfun_info.removeSubmodels();  
    /// \endcode
    ///
    /// The submodel being added is \b assumed to be both valid and not in use by
    /// any ParameterMgr.  These conditions are checked, and the program exits
    /// if they are not met. 
    /// \retval 0  Successfully added sub_mod
    /// \retval 1+ Submodel specific error codes 
    int addSubModel(Submodel* sub_mod);

    ///
    /// Disconnects from the SubModels.  Removes all
    /// Submodels from internal storage and unlinks the
    /// submodels.  This should \b only be done after the use of
    /// the Submodel to ParameterMgr connection is finished
    /// (such as when maximization is totally complete), as it
    /// can't be re-established.
    void removeSubmodels();

    ///
    /// Returns true if there are any submodels currently linked
    /// to this ParameterMgr object.
    bool hasLinkedSubmodels() const;

    ///
    /// Creates a new submodel according to the type given.
    /// \param t The submodel type (use SMType<foo_type>) when invoking.
    /// \returns A non-const shared pointer to the newly created submodel.
    template<class SMType> boost::shared_ptr<typename SMType::sm_type> createSubmodel(SMType t);

    ///
    /// Returns a non-const shared pointer to the named submodel.
    /// Please note: If the submodel is not found, it will cause an internal error.
    /// \param name The name of the requested submodel
    NewSubmodelShPtr getSubmodel(const string& name);

    ///
    /// Returns a const shared pointer to the named submodel.
    /// Please note: If the submodel is not found, it will cause an internal error.
    /// \param name The name of the requested submodel
    NewSubmodelShCstPtr getSubmodel(const string& name) const;

    ///
    /// Returns a const reference to the named, typed submodel.
    /// Please note: This function static_cast's the return type for you. If the
    /// name of the requested submodel returns an object that cannot be downcast
    /// to the requested type, you'll have a runtime error.
    /// Please note: If the submodel is not found, it will cause an internal error.
    /// \param name The name of the submodel
    /// \param t The type descriptor of the submodel
    template<class SM> const typename SM::sm_type& getSubmodel(const string & name, SM t) const;

    ///
    /// Returns a non-const reference to the named, typed submodel.
    /// Please note: This function static_cast's the return type for you. If the
    /// name of the requested submodel returns an object that cannot be downcast
    /// to the requested type, you'll have a runtime error.
    /// Please note: If the submodel is not found, it will cause an internal error.
    /// \param name The name of the submodel
    /// \param t The type descriptor of the submodel
    template<class SM> typename SM::sm_type& getSubmodel(const string & name, SM t);

  //@}

  /// @name Basic informative functions
  //@{

    ///
    /// Indicates whether or not the named parameter has been added to this object.
    /// \param group_name The name of the group to which the parameter should belong.
    /// \param param_name The name of the parameter.
    /// \retval true The parameter \b does exist.
    /// \retval false The parameter does \b not exist.
    bool doesParamExist(string group_name, string param_name) const;

    ///
    /// Returns the ID number of the parameter named by \c param_name.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param param_name The name of the parameter.
    int getParamID(string group_name, string param_name) const;

    ///
    /// Returns the ID number of the parameter whose group_id = \c group_id.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param group_id The group ID number of the parameter.
    int getParamID(string group_name, int group_id) const;

    ///
    /// Returns the total number of parameters.
    int getParamCount() const;

    ///
    /// Returns the number of parameters in the named group.
    /// \param group_name The name of the group whose number of parameters will be returned.
    int getParamCount(string group_name) const;

    ///
    /// Returns the total number of parameters that will be estimated.
    int getEstimatedParamCount() const;
                
    ///
    /// Returns the number of parameter groups.
    int getGroupCount() const;

  //@}

  /// @name Parameter accessor/mutators
  //@{

    ///
    /// Returns the current estimate of the given paramid. This function is present
    /// in case you need to use boost::lambda::bind to access parameter estimates.
    double getEst(int param_id) const;

    ///
    /// Returns a reference to the named parameter's current estimate.
    /// \param group_name The name of the group to which this parameter belongs.
    /// \param param_name The name of parameter whose estimate will be returned.
    double& operator() (string group_name, string param_name);

    ///
    /// Returns the value of to the named parameter's current estimate.
    /// \param group_name The name of the group to which this parameter belongs.
    /// \param param_name The name of parameter whose estimate will be returned.
    double operator() (string group_name, string param_name) const;

    ///
    /// Returns a reference to the parameter's current estimate.
    /// \param param_id The id number of the parameter whose estimate will be returned.
    double& operator() (int param_id);

    ///
    /// Returns a the value of the parameter's current estimate.
    /// \param param_id The id number of the parameter whose estimate will be returned.
    double operator() (int param_id) const;

    ///
    /// Returns a reference to the named Parameter.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param param_name The name of the parameter.
    Parameter& getParameter (string group_name, string param_name);

    ///
    /// Returns a const reference to the named Parameter.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param param_name The name of the parameter.
    const Parameter& getParameter (string group_name, string param_name) const;

    ///
    /// Returns a reference to the indicated Parameter.
    /// \param param_id The ID number of the parameter.
    Parameter& getParameter (int param_id);

    ///
    /// Returns a const reference to the indicated Parameter.
    /// \param param_id The ID number of the parameter.
    const Parameter& getParameter (int param_id) const;

    ///
    /// Returns a reference to the indicated Parameter.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param group_id The group ID number of the parameter.
    Parameter& getParameter(string group_name, int group_id);

    ///
    /// Returns a const reference to the indicated Parameter.
    /// \param group_name The name of the group to which the parameter belongs.
    /// \param group_id The group ID number of the parameter.
    const Parameter& getParameter(string group_name, int group_id) const;

  //@}

  /// @name Parameter iterators
  //@{

    ///
    /// Returns a non-const begin iterator for the list of all parameters.
    ParameterIterator getParamBegin();

    ///
    /// Returns a const begin iterator for the list of all parameters.
    ParameterConstIterator getParamBegin() const;

    ///
    /// Returns a non-const end iterator for the list of all parameters.
    ParameterIterator getParamEnd();

    ///
    /// Returns a const end iterator for the list of all parameters.
    ParameterConstIterator getParamEnd() const;

    ///
    /// Returns a non-const begin iterator for the list of all parameters in the named group.
    /// \param group_name The name of the group whose begin iterator will be returned.
    ParameterIterator getParamBegin(string group_name);

    ///
    /// Returns a const begin iterator for the list of all parameters in the named group.
    /// \param group_name The name of the group whose begin iterator will be returned.
    ParameterConstIterator getParamBegin(string group_name) const;

    ///
    /// Returns a non-const end iterator for the list of all parameters in the named group.
    /// \param group_name The name of the group whose end iterator will be returned.
    ParameterIterator getParamEnd(string group_name);

    ///
    /// Returns a const end iterator for the list of all parameters in the named group.
    /// \param group_name The name of the group whose end iterator will be returned.
    ParameterConstIterator getParamEnd (string group_name) const;
    
    // Added 8-29-7. djb
    //
    void  dumpParameterLookupTable() const;    

  //@}

  // Added JA July -09 
  // Are all parameters fixed ?
     bool allfixed();

  protected:

    ///
    /// Gets this ParameterMgr's submodels to fully add themselves to the mgr.
                int finalizeSubmodels();

    // Copy function (for assignment operator and copy constructor).
                //
    // WARNING:  Cannot be used if either ParameterMgr object
    // involved is in use (see isInUse()).  If this is
    // attempted, the program will exit.
    void copy(const ParameterMgr &);

    // Copy function (to copy an In-Use ParameterMgr into a not-in-use one.
    //
    // This function is only used in the Maximize() function. It copies all contents
    // of a ParameterMgr *except* linkages to submodels.
    //
    // This function is used according the following flowchart:
    //
    // 1. Invoke addSubModel() on ParameterMgr
    // 2. Invoke Maximize()
    // 2.1. Maximize() maximizes the function.
    // 2.2. Maximize constructs a Results object, into which the ParameterMgr is copied. It is at this 
    //      point that copyInUseParameterMgr is invoked.
    // 3. Invoke removeSubModels() on ParameterMgr.
    void copyInUseParameterMgr(const ParameterMgr &);

    // update (for Maxfun synchronization)
    int update(vector<double> & params);

                // Transforms the internal/external parameters
                void transformParameters();

                // Immediately prior to maximization, this function sticks in the correct initial estimates for internal/external parameters:
                void calculateInternalInitialEstimates();

    // update (for Maxfun synchronization of dependent parameters)
    int updateDependents(vector<double> & params);

    // dump() will print a list of the parameters and their current
    // estimates to stdout.
    void dump(DebugCfg & debug) const;

    // Returns a list of all the group names.
    vector<string> getGroupNames() const;

  public:
    // Returns the group names in the order in which they were added.
    vector<string> getOrderedGroupNames() const;

  protected:
    // Returns whether or not a group exists.
    bool doesGroupExist(string group_name) const;

    // Returns const & non-const param_group's.
    param_group       & getGroup(string group_name);
    const param_group & getGroup(string group_name) const;
    
    // Removes a submodel from the list
    void removeSubmodel(Submodel* sm);



  private:

  // Master list of parameters:
    param_vector params;

  // Lookup table for parameter names:
    std::map<pair<string, string>, int> param_lookup_table;

  // Sub-lists of parameters:
    mutable std::map<string, param_group> groups;

  // Group ordering:
    mutable vector<string> group_ordering;
    
  // List of Submodels
    list<Submodel*> my_submodels;

  // List of NewSubmodels
     vector<NewSubmodelShPtr> my_new_submodels;

  // Parameter transformation
     vector<TransformerInfo> my_transformer_infos;

  // Parameter calculators:
     std::map<int, ParamCalculatorShPtr> my_param_calcs;
};

//================================================================
//  NewSubmodel INLINES
//================================================================

inline NewSubmodel::NewSubmodel(const string & name, cerrorstream & errors) 
{ 
  my_finalized  = false;
  my_errors     = errors;
  my_name       = name;
  my_parameters . clear(); 
}

inline NewSubmodel::NewSubmodel(const NewSubmodel & other)
{
  my_finalized  = other.my_finalized;
  my_name       = other.my_name;
  my_errors     = other.my_errors;
  my_parameters = other.my_parameters;
}

inline NewSubmodel& 
NewSubmodel::operator=(const NewSubmodel & other)
{
  my_finalized  = other.my_finalized;
  my_name       = other.my_name;
  my_errors     = other.my_errors;
  my_parameters = other.my_parameters;

  return *this;
}

inline void NewSubmodel::setName(const string & name) { my_name = name; }

inline const string & NewSubmodel::getName() const { return my_name; }

inline bool NewSubmodel::isFinalized() const { return my_finalized; }

inline void NewSubmodel::setFinalized(bool finalized) { my_finalized = finalized; }

inline double & NewSubmodel::getParam(int local_id)       { return my_mgr->operator()(my_parameters[local_id].index); }
inline double   NewSubmodel::getParam(int local_id) const { return my_mgr->operator()(my_parameters[local_id].index); }

inline void NewSubmodel::setParameterMgr(ParameterMgr * mgr) { my_mgr = mgr; }

inline int
NewSubmodel::addParametersToMgr()
{
  for(size_t i = 0; i < my_parameters.size(); i++)
    my_mgr->addParameter(my_parameters[i]);
        
  return 0;
}

inline int 
NewSubmodel::finalizeConfiguration()
{
  return 0;
}

inline void 
NewSubmodel::setAdvancedParameterOptions()
{
}
                                                                                                
inline const ParameterMgr & 
NewSubmodel::getParameterMgr() const
{
  return *my_mgr;
}


class FooSubmodel : public NewSubmodel
{
  friend class ParameterMgr;
  
public:

  FooSubmodel() : NewSubmodel("Foo") 
  {
    my_parameters.push_back(ParameterInput("Foo", "tempy", Parameter::FIXED, 1.0, -1.0, 1.0));
  }
  
  FooSubmodel(const FooSubmodel & other) : NewSubmodel(other) {}

  virtual int update() { return 0; }
    
  virtual NewSubmodelShPtr clone() { return NewSubmodelShPtr(new FooSubmodel(*this)); }
        
};

//=======================================================================
//  hasLinkedSubmodels()
//=======================================================================

inline bool
ParameterMgr::hasLinkedSubmodels() const
{
  return !my_submodels.empty();
}

//=======================================================================
//  isInUse()
//=======================================================================
inline bool
ParameterMgr::isInUse() const
{
  return hasLinkedSubmodels();
}

//=======================================================================
//  getParamID(...) #1
//=======================================================================
inline int
ParameterMgr::getParamID(string group_name, string param_name) const
{
  std::map<pair<string, string>, int>::const_iterator i = param_lookup_table.find(make_pair(group_name, param_name));

  if(i == param_lookup_table.end())
  {
    cout << "Error: Parameter '" << param_name << "' in group '" << group_name << "' does not exist." << endl;
    exit(1);
  }

  return i->second;
}

//=======================================================================
//  getParamID(...) #2
//=======================================================================
inline int
ParameterMgr::getParamID(string group_name, int group_id) const
{
  if(!doesGroupExist(group_name))
  {
    cout << "Error: Group " << group_name << " does not exist!" << endl;
    exit(1);
  }

  const param_group & g = getGroup(group_name);

  if((size_t)group_id >= g.size())
  {
    cout << "Error: Group id " << group_id << " too large for group " << group_name << endl;
    exit(1);
  }

  return g[group_id];
}

//=====================================================================
//  createSubmodel(...)
//=====================================================================
template<class SMType> 
inline boost::shared_ptr<typename SMType::sm_type> 
ParameterMgr::createSubmodel(SMType t) 
{
  boost::shared_ptr<typename SMType::sm_type> sm_sh_ptr = boost::shared_ptr<typename SMType::sm_type>(new typename SMType::sm_type);

  sm_sh_ptr->setParameterMgr(this);

  my_new_submodels.push_back(sm_sh_ptr);

  return sm_sh_ptr;
}

//=====================================================================
//  getSubmodel(...) TEMPLATIZED CONST
//=====================================================================
template<class SM> 
inline const typename SM::sm_type& 
ParameterMgr::getSubmodel(const string & name, SM t) const
{
  return *(static_cast<const typename SM::sm_type *> (getSubmodel(name).get()));
}
                
//=====================================================================
//  getSubmodel(...) TEMPLATIZED NON-CONST
//=====================================================================
template<class SM> 
inline typename SM::sm_type& 
ParameterMgr::getSubmodel(const string & name, SM t)
{
  return *(static_cast<typename SM::sm_type *> (getSubmodel(name).get()));
}

}} // End namespace

MEM_COUNT_BEGIN(SAGE::MAXFUN::ParameterMgr)
{
  size_t x = 0;
  
  for(size_t i = 0; i < t.params.size(); ++i)
    x += get_mem(t.params[i]);
    
  for(std::map<pair<string, string>, int>::const_iterator i = t.param_lookup_table.begin(); i != t.param_lookup_table.end(); ++i)
  {
    x += get_mem(i->first.first) + get_mem(i->first.second) + get_mem(i->second);
  }
 
  for(std::map<string, SAGE::MAXFUN::param_group>::const_iterator i = t.groups.begin(); i != t.groups.end(); ++i)
  {
    x += get_mem(i->first) + get_mem(i->second);
  }

  x += get_mem(t.group_ordering);

  x += sizeof(SAGE::MAXFUN::ParameterMgr);

  return x;
}
MEM_COUNT_END

#endif
