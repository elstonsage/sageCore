#ifndef SUBMODEL_H
#define SUBMODEL_H

#include "error/errorstream.h"
#include "maxfun/sub_model.h"
#include "maxfunapi/Datatypes.h"
#include "maxfunapi/ParameterInput.h"

namespace SAGE   {
namespace MAXFUN {

/** \class Submodel
 *  \brief Provides encapsulated functionality for a specialized group of parameters.
 *
 * \par Introduction
 *
 * A Submodel helps you organize a related group of parameters, as well as 
 * functions that connect those parameters, into a single organization unit. 
 *
 * For instance, let's say you want to be able to estimate transformation parameters
 * in your program. You could implement your own transformation functions, and add 
 * both power and shift parameters directly to ParameterMgr.
 *
 * Alternately, you could use one of the various already-written Submodel's, including:
 * 
 * MAXFUN::TransformationSubmodel
 *
 * Finally, you could create a transformation submodel (derived from Submodel)
 * that creates and manages the power and shift parameters. In addition, this submodel
 * would be able to provide functionality for transformation.
 *
 * \par Using an extant submodel
 *
 * Please see the detailed description for MAXFUN::ParameterMgr for more information on
 * using extant submodels.
 *
 * \par Creating your own Submodel
 *
 * First, create your own derived class from Submodel, and put in the required
 * virtual interface:
 *
 * \code
 * class mysubmodel : public Submodel
 * {
 *   public:
 *     mysubmodel();
 *
 *   protected:
 *     virtual int update();
 * };
 * \endcode
 *
 * Next, use the my_parameters vector (inherited from MAXFUN::Submodel) to set up whatever
 * initial states/options you need for your parameters. This can be done in your constructor, or in
 * any other setup helper functions you create. In the following example, you can see how to set up
 * a parameter in your submodel's constructor. Please note that the my_parameters vector is a vector
 * of MAXFUN::ParameterInput objects; remember to consult the documentation for this object.
 * 
 * [mysubmodel.cpp]
 * \code
 * #include "maxfun/ParameterMgr.h"
 * #include "maxfun/MaxfunDatatypes.h"
 * #include "mysubmodel.h"
 *
 * mysubmodel::mysubmodel()
 * {
 *   MAXFUN::ParameterInput newparam("mysubmodel", "y", MAXFUN::INDEPENDENT_FUNCTIONAL, 0.5, 0.0, 1.0);
 *   my_parameters.push_back(newparam);

 *   MAXFUN::ParameterInput newparam2("mysubmodel", "z", MAXFUN::DEPENDENT, 0.5, 0.0, 1.0);
 *   my_parameters.push_back(newparam2);
 * };
 * \endcode
 *
 * Now, define the virtual interface. Note that you can use the getParam() function to fetch a reference
 * to the parameter's current estimate. Also, note that the index number identifying which parameter you want
 * corresponds to the order in which you added that parameter to the my_parameters vector.
 *
 * \code
 * int mysubmodel::update()
 * {
 *   getParam(1) = 1 - getParam(0);
 *   return 0;
 * }
 * \endcode
 *
 * Now, add whatever specialized functions are necessary for your submodel.
 *
 * \par A complete example
 *
 * Let's say you wanted a submodel that stores two proportions (\c x1 and \c x2 ) whose sum must equal 1.0.
 * In addition, this submodel will have a function that reports the product of the two proportions.
 *
 * To achieve this, you will need to:
 *
 * \c 1 Create your own derived class from Submodel.
 *
 * \c 2 Add functions that can report the current estimates for each parameter.
 *
 * \c 3 Define the update() function such that it makes sure the sum of the two parameters is 1.0.
 *
 * \c 4 Define a function that reports the product of the two proportions.
 *
 * To achieve all this, let's look at a complete example:
 *
 * [mysubmodel.h]
 * \code
 * #include "maxfun/Submodel.h"
 *
 * class mysubmodel : public MAXFUN::Submodel
 * {
 *   public:
 *
 *     mysubmodel();
 *
 *     double getX1() const;
 *     double getX2() const;
 *
 *     double getProduct() const;
 *
 *   protected:
 *     virtual int update();
 * };
 * \endcode
 *
 * [mysubmodel.cpp]
 * \code
 * #include "maxfun/ParameterMgr.h"
 * #include "maxfun/MaxfunDatatypes.h"
 * #include "mysubmodel.h"
 *
 * mysubmodel::mysubmodel()
 * {
 *   MAXFUN::ParameterInput paramx1("submod", "x1", MAXUFN::INDEPENDENT_FUNCTIONAL, 0.5, 0.0, 1.0);
 *   my_parameters.push_back(paramx1);
 *
 *   MAXFUN::ParameterInput paramx2("submod", "x2", MAXUFN::DEPENDENT, 0.5, 0.0, 1.0);
 *   my_parameters.push_back(paramx2);
 * }
 *
 * double mysubmodel::getX1() const 
 * { 
 *   return getParam(0); 
 * }
 *
 * double mysubmodel::getX2() const 
 * { 
 *   return getParam(1); 
 * }
 *
 * double mysubmodel::getProduct() const 
 * { 
 *   return getParam(0) * getParam(1); 
 * }
 *
 * int mysubmodel::update()
 * {
 *   getParam(1) = 1.0 - getParam(0);
 * }
 * \endcode
 *
 * \par More features
 *
 * The standard method for setting up your parameters is to resize and alter the my_parameters
 * vector dynamically. As the user invokes different setup functions, the submodel modifies the
 * my_parameters vector immediately.
 *
 * If, however, you don't want to continuously modify the my_parameters vector dynamically, you
 * can instead do so in the finalizeConfiguration() function. If you implement this virtual function,
 * it is expected that you will resize and setup the my_parameters vector in preparation for maximization.
 *
 * Also, if you want to use any of MAXFUN::Parameter's advanced setup features, you can do so in
 * the setAdvancedParameterOptions() function. Please see the function's documentation for more information.
 */
class Submodel
{
  public:

    Submodel(cerrorstream& errors = sage_cerr);
    Submodel(const Submodel& other);
    Submodel&  operator=(const Submodel& other);
  
    // - Virtual destructor is required due to virtual functions.
    //
    virtual ~Submodel();

    // isLinked indicates that the submodel is currently connected to a
    // ParameterMgr object.
    bool isLinked() const;

  protected:

    friend class ParameterMgr;

  /// @name Virtual interface
  //@{

    ///
    /// Updates any relevant internal state from ParameterMgr.
    ///
    /// Note: You are \b required to define this function in a derived class.
    /// \retval 0 Internal state was updated succesfully.
    /// \retval 1 internal state was \b not updated successfully.
    virtual int update() = 0;

    ///
    /// If for any reason you elected \b not to resize the my_parameters vector
    /// in your initial setup functions, you can do so within the finalizeConfiguration()
    /// function. This function will be called immediately prior to actually \b adding
    /// the contents of my_parameters to Maxfun.
    /// \retval 0 Configuration was finalized successfully.
    /// \retval 1 Configuration was \b not finalized successfully.
    virtual int finalizeConfiguration();

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
     *    getParameterMgr()->getParameter("submod1", "X1").IncludepValue() = false;
     *  }
     *  \endcode
     */
    virtual void setAdvancedParameterOptions();

  //@}

  /// @name Runtime maximization functions
  //@{

    ///
    /// Within run-time maximization, you should fetch parameter values by using
    /// the getParam() function.
    /// \param local_id The index (within my_params) of the parameter you want.
    /// \returns A double reference to the parameter's current estimate.
    double & getParam(int local_id);

    ///
    /// Const version of getParam().
    double getParam(int local_id) const;

  //@}

    ///
    /// my_params is a vector of MAXFUN::ParameterInput's, designed to help
    /// you to set up initial options for your parameters.
    vector<ParameterInput> my_parameters;

    /// \internal
    ///
    /// Adds the necessary parameters to ParameterMgr.
    ///
    /// Note: You are \b required to define this function in a derived class.
    /// \retval 0 Parameters were added successfully.
    /// \retval 1 Parameters were \b not added successfully.
    int addParametersToParameterMgr();

    /// \internal
    ///
    /// Sets the my_info pointer to the indicated ParameterMgr
    /// \param mi The ParameterMgr to which this Submodel should be linked
    void linkToParameterMgr(ParameterMgr* mi);

    /// \internal
    ///
    /// Sets the my_info pointer to NULL, "unlinking" it
    void unlinkFromParameterMgr();
    
    /// \internal
    ///
    /// Returns a non-const pointer to the ParameterMgr to which this object is linked.
    /// Please note that if this object has not actually been linked to a ParameterMgr,
    /// this function returns NULL.
    ParameterMgr* getParameterMgr();

    /// \internal
    ///
    /// Returns a const pointer to the ParameterMgr to which this object is linked.
    /// Please note that if this object has not actually been linked to a ParameterMgr,
    /// this function returns NULL.
    const ParameterMgr * getParameterMgr() const;

    /// \internal
    ///
    /// Stream for reporting errors.
    mutable cerrorstream  my_errors;

  private:

    /// \internal
    ///
    /// The info class currently using the submodel.  Is NULL when not in
    /// use.  This pointer is private so derived classes cannot (should not)
    /// modify it.
    ParameterMgr* my_info;
};

} // end MAXFUN namespace
} // end SAGE namespace

#include "maxfunapi/Submodel.ipp"

#endif


