#ifndef FUNC_PYTHON_INTERFACE_H
#define FUNC_PYTHON_INTERFACE_H

#include <string>
#include <map>
#include <vector>
#include <cctype>
#include <limits>
#include <cassert>
#include <sstream>

#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"

#include "Python.h"
#include "compile.h"
#include "eval.h"
#include "func/FunctionParser.h"

#include <setjmp.h>
#include <signal.h>
#include <unistd.h>

//#define PY_ERROR_MSG
#define DISABLE_TIMER 0
//#define NO_TIMER
#define TIMED_OUT 113

extern jmp_buf context;

namespace SAGE {
namespace FUNC {

/// \brief Provides interface to embedded python interpreter.
class PythonInterface
{
public:

  /// Describes any exception to do a problem with the underlying python
  class PythonException : public std::exception {};

  /// @name Constructor/destructor
  //@{

    ///
    /// Constructor.
    /// \param errors Errorstream for this object
    PythonInterface(cerrorstream errors);
          
    ///
    /// Destructor.
    ~PythonInterface();
          
  //@}
        
    /// ???
    static void timeout(int return_value);
      
    ///
    /// Returns a vector of all the symbolic names used in the expression.
    ///
    /// If there is any problem encountered, returns an empty vector.
    /// \param expr The expression in which to search for names
    /// \throws A PythonException is it encounters any underlying Python
    /// problems during execution.
    std::vector<std::string> getNameList(const std::string & expr);

    ///
    /// Examines the python environment to see if the name is present in it.
    bool isNameInEnvironment(const std::string & name) const;

    ///
    /// ???
    bool compile(std::string expr);                 

  /// @name Adding variables to the globals dictionary
  //@{
  
    ///
    /// Adds a variable named 'name', whose value is calculated on the basis of the
    /// given Python expression.
    /// \param name The name of the new variable
    /// \param expr The text of the Python code which will calculate the value
    /// \param time_limit A time limit for executing the Python code
    /// \returns \c true if successful, \c false otherwise
    bool calculateValueInEnvironment(const std::string& name, const std::string& expr, unsigned int time_limit);

    ///
    /// Adds the given name as a variable with value 'None' to the globals hashtable.
    /// \param name The name of the variable
    /// \returns \c true if successful, \c false otherwise.
    bool addNoneToEnvironment(const std::string & name) const;

    /// 
    /// Adds the given name/value pair to the globals hashtable (where value is a double).
    /// \param name The name of the pair
    /// \param value The value of the pair
    /// \returns \c true if successful, \c false otherwise.
    bool addDoubleToEnvironment(const std::string & name, double value) const;

    ///
    /// Adds the given name/value pair to the globals hashtable (where value is a string).
    /// \param name The name of the pair
    /// \param value The value of the pair
    /// \returns \c true if successful, \c false otherwise.
    bool addStringToEnvironment(const std::string & name, const std::string & value) const;

    ///
    /// Adds an allele list to the Python environment. This is equivalent to the following Python code:
    /// \code
    /// name = []
    /// name.append(allele1)
    /// name.append(allele2)
    /// \endcode
    /// \param name The name of the allele list
    /// \param allele1 The name of the first allele
    /// \param allele2 The name of the second allele2
    /// \returns \c true if successful, \c false otherwise
    bool addAlleleListToEnvironment(const std::string & name, const std::string & allele1, const std::string & allele2) const;

  //@}
  
    ///
    /// Executes the given string as a python statement.
    /// \returns \c true if successful, \c false otherwise
    bool runPythonString(const std::string & source) const;

        ///
        /// ???
        double run(const string& expr, unsigned int time_limit);
        
private:

        ///
        /// ???
        std::string py_error_msg() const;

  /// @name Adding functions to the Python environment
  //@{

    ///
    /// Adds the 'extract_names' function to the Python environment.
    ///
    /// Given a tree representing ???, returns a hash table where every
    /// key is a named element from the tree.
    void define_extract_names() const;

    ///
    /// ???
    void define_marker_models() const;

  //@}

    ///
    /// Create a python dictionary of names from an expression.
    /// The dictionary is named 'names' and is added to the Python environment.
    ///
    /// Note: Assumes 'extract_names' function has been defined in the
    /// Python global namespace (this is accomplished via define_extract_names() ).
    /// \returns \true if successful, \c false otherwise
    bool extract_names(const std::string& expr);
        
        ///
        /// ???
        void empty_python_path() const;

        ///
        /// ???
        void restrict_python() const;

        ///
        /// ???
        void unrestrict_python() const;

        ///
        /// ???
        void exit_on_timeout(const std::string& expr);
        
        // =========================================
        // =========== DATA MEMBERS ================
        // =========================================

        /// Errorstream for this object
        cerrorstream my_errors;

        /// Dictionary of names in Python main module.
        PyObject* my_globals;           

        /// Python byte code.
        PyObject* my_obj_code;          
    };

} // End namespace FUNC
} // End namespace SAGE

#endif
