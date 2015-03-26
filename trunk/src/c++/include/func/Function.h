#ifndef FUNC_H
#define FUNC_H

#include <string>
#include <map>
#include <vector>
#include <cctype>
#include <limits>
#include <cassert>
#include <setjmp.h>
#include <signal.h>
#include <unistd.h>

#include "LSF/LSFfile.h"
#include "rped/rped.h"
#include "error/errorstream.h"
#include "error/errormanip.h"

#include "Python.h"
#include "compile.h"
#include "eval.h"

#include "func/FunctionParser.h"
#include "func/PythonInterface.h"

//#define PY_ERROR_MSG
#define DISABLE_TIMER 0
//#define NO_TIMER
#define TIMED_OUT 113

extern jmp_buf context;

namespace SAGE {
namespace FUNC {

/// \brief Creates a new variable from an arbitrary function (possibly based on existing variables)
///                                                                          
/// The function expression is represented by a string.
/// This class embeds the Python language interpreter (via the python
/// class) to parse and evaluate the expression.
///
/// CAN'T HANDLE CASES WHERE A VARIABLE NAME IN THE FUNCTION 
/// EXPRESSION DOES NOT CONFORM TO PYTHON NAME REQUIREMENTS, IE, 
/// CONTAINS ONLY ALPHANUMERIC CHARACTERS OR UNDERSCORE 
/// - DOES NOT START W. A DIGIT. 
///
/// NAME OF THE CREATED VARIABLE MAY NOT BE THAT OF AN EXISTING
/// VARIABLE.
class Function
{
  public:
  
    ///
    /// Creates the requested trait.
    /// \param errors The error stream to use for messages.
    /// \param mp The multi pedigree in which to create the trait.
    /// \param par_data The expression information.
    static void create_trait(
      cerrorstream errors,
      RPED::MultiPedigree & mp,
      const FunctionParser::TraitData & par_data);

    ///
    /// Creates the requested trait.
    /// \param errors The error stream to use for messages.
    /// \param mp The multi pedigree in which to create the trait.
    /// \param trait_name The name of the new trait to create
    /// \param par_data The python expression as a string (eg: "a + b * 5.0")
    ///
    /// NOTE: This function creates a QUANTITATIVE trait.
    static void create_trait(
      cerrorstream errors,
      RPED::MultiPedigree & mp,
      const std::string & trait_name,
      const std::string & function_code);

  private:

    typedef std::map<std::string, size_t> input_variable_map;
  
    friend class class_aggregate_process;
    
    /// Constructor/destructor.
    //@{
    
      ///
      /// Constructor.
      /// \param errors Errorstream for this object
      Function(cerrorstream errors);
    
    //@}

    ///
    /// Create a new variable for the multipedigree whose values are a function
    /// of existing multipedigree variables as specified by expr.
    /// \param mp Multipedigree whose data will be analyzed
    /// \param par_data Object describing the how to calculate the new variable
    void create(RPED::RefMultiPedigree& mp, const FunctionParser::TraitData & par_data); 
  
    ///
    /// Extract names from the expression and remember those of existing
    /// traits or markers.
    /// \param mp Multipedigree whose data will be analyzed
    /// \param expr The python expression.
    /// \param calculator The calculator object.
    /// \returns \c True if there are any duplicated variables; false otherwise.
    bool find_variables_used(
      const RPED::RefMultiPedigree& mp, 
      const std::string& expr, 
      PythonInterface& calculator);

    ///
    /// Returns \c true if the trait is valid, \c false otherwise.
    bool  valid_trait  (const RPED::RefMultiPedigree& mp, const string& name) const;

    bool no_duplicates() const;

    bool set_constants(PythonInterface& calculator, 
                        const FunctionParser::TraitData::constant_vector& constants,
                        unsigned int time_limit);

    ///
    /// Sets the required member-specific variables (traits, markers, etc.) in the PythonInterface
    /// immediately prior to invoking the Python function.
    /// \param interface The PythonInterface that will hold the values
    /// \param member The RPED::Member in question
    /// \returns \c true if successful, \c false otherwise
    bool setMemberSpecificVariables(
          const PythonInterface & py_interface, 
          const RPED::Member & member) 
          const;

    ///
    /// Lookup trait values for a given pedigree member for each trait   
    /// in expression to be evaluated.  Set trait values in Python
    /// environment.  If any trait has a missing value, return false.
    bool setMemberSpecificTraits(const PythonInterface& calculator, const RPED::Member & member) const;


    bool setMemberSpecificStrings(const PythonInterface& calculator, const RPED::Member & member) const;


    bool  setMemberSpecificMarkers(const PythonInterface& calculator, const RPED::Member & member) const;


    void  rollback(RPED::RefMultiPedigree& mp, 
                   RPED::RefMultiPedigree::pedigree_iterator& l_iter);

    /// Error handling.
    enum error_msg  { existing_name_msg, bad_syntax_msg, eval_failure_msg, 
                      internal_error_msg, bad_const_expr_msg, missing_info_msg };
    
    void  write_error_msg(error_msg msg,
                          const std::string& name,
                          const std::string& expr );    
                                        
    void  write_error_msg(const std::string& message);
      
    /// Data members.
    cerrorstream        my_errors;
    
    input_variable_map  my_traits_used;   /// Traits used in expression.
                                          /// <trait name, trait index>
    input_variable_map  my_strings_used;  /// Strings used in expression.
                                          /// <string name, string index>
    input_variable_map  my_markers_used;  /// Markers used in expression.
                                          /// <marker name, trait index>
};

} // End namespace FUNC
} // End namespace SAGE

#endif
